#ifndef SWBAM_APP_ORDERED_BODY_H
#define SWBAM_APP_ORDERED_BODY_H

#include "bounded-cache.h"
#include <limits>
#include <vector>

namespace swbam_usage {

// The prefix is borrowed and must remain unchanged after Activate(). Only
// compressed bytes belong here; the caller owns record decisions and keys.
class OrderedBody : public swbam::RankBodySink {
public:
    OrderedBody(swbam::RankBodySink *ordered_prefix, size_t metadata_limit)
        : prefix_(ordered_prefix), metadata_limit_(metadata_limit) {}

    int Activate() {
        if (!prefix_ || sealed_ || failed_) return Error(EINVAL, "invalid activation");
        if (!active_) {
            prefix_size_ = prefix_->size();
            total_size_ = prefix_size_;
            active_ = true;
        }
        return 0;
    }

    int AppendOrdered(uint64_t key, const void *data, size_t bytes) {
        if (!prefix_ || sealed_ || failed_ || (!data && bytes))
            return Error(EINVAL, "invalid append or sealed body");
        if (!bytes) return 0;
        if (have_prefix_key_ && key <= last_prefix_key_)
            return Error(EINVAL, "key precedes the committed prefix");
        if (bytes > UINT64_MAX - size()) return Error(EOVERFLOW, "body size overflow");
        if (!active_) {
            if (prefix_->Append(data, bytes)) { failed_ = true; return -1; }
            last_prefix_key_ = key;
            have_prefix_key_ = true;
            return 0;
        }
        if (append_serial_ == UINT64_MAX) return Error(EOVERFLOW, "key count overflow");
        const uint64_t serial = append_serial_++;
        const char *src = static_cast<const char *>(data);
        while (bytes) {
            if (ReserveOne(&extents_)) { failed_ = true; return -1; }
            if (segments_.empty() || segments_.back()->file.size() == SegmentBytes()) {
                if (ReserveOne(&segments_) || !CanAllocate(sizeof(Segment))) {
                    failed_ = true;
                    return -1;
                }
                try { segments_.emplace_back(new Segment); }
                catch (...) { failed_ = true; return Error(ENOMEM, "segment allocation failed"); }
                ++live_segments_;
            }
            Segment &segment = *segments_.back();
            const size_t count = (size_t)std::min<uint64_t>(
                bytes, SegmentBytes() - segment.file.size());
            const uint64_t offset = segment.file.size();
            if (segment.file.Append(src, count)) {
                const int error = errno ? errno : EIO;
                failed_ = true;
                return Error(error, "compressed temporary write failed");
            }
            // logical_offset holds the append id until sealing, to distinguish
            // a split append from two calls accidentally using the same key.
            Extent extent = {key, offset, count, serial, segments_.size() - 1};
            extents_.push_back(extent);
            total_size_ += count;
            src += count;
            bytes -= count;
        }
        return 0;
    }

    int Append(const void *, size_t) {
        return Error(EINVAL, "AppendOrdered requires an input block key");
    }

    int Flush() {
        if (!prefix_ || failed_) return Error(EINVAL, "cannot seal failed body");
        if (sealed_) return 0;
        if (prefix_->Flush()) { failed_ = true; return -1; }
        if (!active_) {
            prefix_size_ = total_size_ = prefix_->size();
        } else {
            std::sort(extents_.begin(), extents_.end(),
                [](const Extent &a, const Extent &b) {
                    if (a.key != b.key) return a.key < b.key;
                    if (a.segment != b.segment) return a.segment < b.segment;
                    return a.offset < b.offset;
                });
            for (size_t i = 1; i < extents_.size(); ++i) {
                if (extents_[i - 1].key == extents_[i].key &&
                    extents_[i - 1].logical_offset != extents_[i].logical_offset) {
                    failed_ = true;
                    return Error(EINVAL, "duplicate compressed block key");
                }
            }
            size_t out = 0;
            for (size_t i = 0; i < extents_.size(); ++i) {
                const Extent extent = extents_[i];
                if (out && extents_[out - 1].segment == extent.segment &&
                    extents_[out - 1].offset + extents_[out - 1].bytes == extent.offset) {
                    extents_[out - 1].bytes += extent.bytes;
                } else {
                    extents_[out++] = extent;
                }
            }
            extents_.resize(out);
            uint64_t logical = prefix_size_;
            for (auto &extent : extents_) {
                extent.logical_offset = logical;
                logical += extent.bytes;
                ++segments_[extent.segment]->remaining;
            }
            if (logical != total_size_) {
                failed_ = true;
                return Error(EIO, "compressed extent size mismatch");
            }
        }
        sealed_ = true;
        return 0;
    }

    uint64_t size() const {
        return active_ || sealed_ ? total_size_ : (prefix_ ? prefix_->size() : 0);
    }

    bool is_memory() const {
        if (!prefix_) return false;
        return (prefix_size_ == 0 && active_ ? true : prefix_->is_memory()) &&
               (extents_.empty() || IsRam());
    }

    int ReadAt(uint64_t offset, void *data, size_t bytes) const {
        if (!sealed_ || failed_ || offset > total_size_ ||
            bytes > total_size_ - offset || (!data && bytes) ||
            (bytes && offset < released_))
            return Error(EINVAL, "invalid read, unsealed body, or consumed range");
        if (!bytes) return 0;
        char *dst = static_cast<char *>(data);
        if (offset < prefix_size_) {
            const size_t count = (size_t)std::min<uint64_t>(bytes, prefix_size_ - offset);
            if (prefix_->ReadAt(offset, dst, count)) return -1;
            offset += count;
            dst += count;
            bytes -= count;
        }
        if (!bytes) return 0;
        auto it = std::upper_bound(extents_.begin(), extents_.end(), offset,
            [](uint64_t pos, const Extent &extent) { return pos < extent.logical_offset; });
        if (it == extents_.begin()) return Error(EIO, "missing compressed extent");
        --it;
        while (bytes) {
            if (it == extents_.end() || offset < it->logical_offset ||
                offset - it->logical_offset >= it->bytes || !segments_[it->segment])
                return Error(EIO, "missing compressed extent");
            const uint64_t within = offset - it->logical_offset;
            const size_t count = (size_t)std::min<uint64_t>(bytes, it->bytes - within);
            if (segments_[it->segment]->file.ReadAt(it->offset + within, dst, count)) return -1;
            offset += count;
            dst += count;
            bytes -= count;
            ++it;
        }
        return 0;
    }

    int ReleaseBefore(uint64_t offset) const {
        if (!sealed_ || failed_ || offset > total_size_)
            return Error(EINVAL, "invalid release or unsealed body");
        if (offset <= released_) return 0;
        if (released_ < prefix_size_ &&
            prefix_->ReleaseBefore(std::min(offset, prefix_size_))) return -1;
        // Logical order can revisit an earlier physical segment. Reclaim a
        // whole segment only once all of its logical extents have been consumed.
        while (consumed_extents_ < extents_.size()) {
            const Extent &extent = extents_[consumed_extents_];
            if (extent.logical_offset + extent.bytes > offset) break;
            std::unique_ptr<Segment> &segment = segments_[extent.segment];
            if (!segment || !segment->remaining) return Error(EIO, "invalid segment references");
            if (--segment->remaining == 0) {
                segment.reset();
                --live_segments_;
            }
            ++consumed_extents_;
        }
        released_ = offset;
        return 0;
    }

    // Owned index/segment allocations only; borrowed prefix and stored payload
    // belong to their existing cache/storage budgets. Constant-time on hot paths.
    size_t memory_bytes() const {
        return extents_.capacity() * sizeof(Extent) +
               segments_.capacity() * sizeof(std::unique_ptr<Segment>) +
               live_segments_ * sizeof(Segment);
    }

private:
    struct Extent {
        uint64_t key, offset, bytes, logical_offset;
        size_t segment;
    };
    struct Segment {
        CachedTemporary file;
        size_t remaining;
        Segment() : file(0), remaining(0) {}
    };

    static uint64_t SegmentBytes() { return 64ull << 20; }
    static int Error(int error, const char *message) {
        fprintf(stderr, "ERROR: ordered body: %s (%s).\n", message, strerror(error));
        errno = error;
        return -1;
    }
    bool CanAllocate(size_t bytes) const {
        const size_t used = memory_bytes();
        if (used <= metadata_limit_ && bytes <= metadata_limit_ - used) return true;
        fprintf(stderr, "ERROR: ordered body metadata budget exceeded: "
                "used=%zu requested=%zu limit=%zu.\n", used, bytes, metadata_limit_);
        errno = ENOBUFS;
        return false;
    }
    template <typename T>
    int ReserveOne(std::vector<T> *items) {
        if (items->size() < items->capacity()) return 0;
        const size_t used = memory_bytes();
        const size_t available = used <= metadata_limit_ ? metadata_limit_ - used : 0;
        const size_t maximum = std::min(items->max_size(), available / sizeof(T));
        if (items->size() >= maximum) {
            Error(ENOBUFS, "metadata budget cannot grow extent/segment index");
            fprintf(stderr, "  ordered body metadata used=%zu limit=%zu.\n", used, metadata_limit_);
            return -1;
        }
        const size_t capacity = items->capacity();
        const size_t next = capacity > maximum / 2 ? maximum :
            std::min(maximum, std::max<size_t>(16, capacity * 2));
        // Include the old allocation while reserve() creates its replacement.
        if (!CanAllocate(next * sizeof(T))) return -1;
        try { items->reserve(next); }
        catch (...) { return Error(ENOMEM, "metadata allocation failed"); }
        return 0;
    }

    OrderedBody(const OrderedBody &);
    OrderedBody &operator=(const OrderedBody &);
    swbam::RankBodySink *prefix_;
    size_t metadata_limit_;
    uint64_t prefix_size_ = 0, total_size_ = 0, last_prefix_key_ = 0, append_serial_ = 0;
    bool have_prefix_key_ = false, active_ = false, sealed_ = false, failed_ = false;
    std::vector<Extent> extents_;
    mutable std::vector<std::unique_ptr<Segment> > segments_;
    mutable size_t live_segments_ = 0, consumed_extents_ = 0;
    mutable uint64_t released_ = 0;
};

}
#endif
