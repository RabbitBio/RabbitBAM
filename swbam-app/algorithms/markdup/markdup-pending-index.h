#ifndef SWBAM_MARKDUP_PENDING_INDEX_H
#define SWBAM_MARKDUP_PENDING_INDEX_H

#include <cstddef>
#include <cstdint>
#include <map>

namespace swbam_markdup {

// Non-owning ranges. Owners must erase a batch before destroying it.
template <typename Batch>
class PendingIndex {
public:
    bool Insert(Batch *batch) {
        if (batch->ordinal_base > batch->ordinal_end) return false;
        if (batch->ordinal_base == batch->ordinal_end) return true;
        auto next = batches_.lower_bound(batch->ordinal_base);
        if (next != batches_.end() && next->first < batch->ordinal_end)
            return false;
        if (next != batches_.begin()) {
            auto previous = next;
            --previous;
            if (previous->second->ordinal_end > batch->ordinal_base) return false;
        }
        batches_.emplace_hint(next, batch->ordinal_base, batch);
        return true;
    }

    Batch *Find(uint64_t ordinal) {
        if (cached_ && ordinal >= cached_->ordinal_base && ordinal < cached_->ordinal_end)
            return cached_;
        auto next = batches_.upper_bound(ordinal);
        if (next == batches_.begin()) return nullptr;
        --next;
        if (ordinal >= next->second->ordinal_end) return nullptr;
        cached_ = next->second;
        return cached_;
    }

    void Erase(const Batch *batch) {
        if (cached_ == batch) cached_ = nullptr;
        auto found = batches_.find(batch->ordinal_base);
        if (found != batches_.end() && found->second == batch) batches_.erase(found);
    }

    void Clear() { cached_ = nullptr; batches_.clear(); }
    size_t size() const { return batches_.size(); }
    static size_t EntryBytes() {
        return sizeof(typename std::map<uint64_t, Batch *>::value_type) + 4 * sizeof(void *);
    }
    size_t memory_bytes() const { return size() * EntryBytes(); }

private:
    std::map<uint64_t, Batch *> batches_;
    Batch *cached_ = nullptr;
};

} // namespace swbam_markdup
#endif
