import unittest

from markdup_halo_metrics import summarize


def fixture(partial=False):
    lines = ["OptimizedMarkdupMPI finished. ranks=2 blocks=12 records=12",
             "  duplicates pair=2 single=1 marked=3"]
    for rank in range(2):
        lines.extend([
            f"[swbam-processing rank={rank}] storage=ram seconds=1.2 status=0",
            f"[md-halo rank={rank}] enabled=1 partial={int(partial)} "
            f"conservative={int(partial)} probe_blocks=64 final_keep={0 if partial else 2} "
            "final_dup=1 pending=2 marker_pending=3 candidate_bytes=10 result_bytes=8",
            f"[md-comm rank={rank}] halo_enabled=1 conservative={int(partial)} "
            "main_candidates=5 main_markers=2 main_candidate_bytes=100 main_result_bytes=16",
            f"[rank {rank}] markdup_stats mpi_candidate_bytes=110 result_bytes=24",
        ])
    lines.append("Complete the 4.3 OptimizedMarkdupMPI cost 1.1")
    return "\n".join(lines)


class HaloMetricsTest(unittest.TestCase):
    def test_complete_and_partial(self):
        for partial in (False, True):
            result = summarize(fixture(partial), 2, 12, 3)
            self.assertEqual(result["halo"], "partial" if partial else "enabled")
            self.assertEqual(result["main_loop"]["record_candidate_entries"], 6)
            self.assertEqual(result["main_loop"]["marker_entries"], 4)
            self.assertEqual(result["main_loop"]["candidate_payload_bytes"], 200)

    def test_partial_cannot_claim_keep(self):
        with self.assertRaisesRegex(ValueError, "without final KEEP"):
            summarize(fixture(True).replace("final_keep=0", "final_keep=1"), 2, 12, 3)

    def test_inconsistent_partial_ranks(self):
        with self.assertRaisesRegex(ValueError, "partial halo status"):
            summarize(fixture(True).replace("partial=1", "partial=0", 1), 2, 12, 3)

    def test_missing_rank_or_payload_mismatch(self):
        for text in (fixture().replace("[md-comm rank=1]", "[missing rank=1]"),
                     fixture().replace("main_candidate_bytes=100", "main_candidate_bytes=101", 1)):
            with self.assertRaises(ValueError):
                summarize(text, 2, 12, 3)


if __name__ == "__main__":
    unittest.main()
