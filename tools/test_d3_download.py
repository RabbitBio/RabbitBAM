#!/usr/bin/env python3
"""Exercise download recovery without network or large input files."""

import hashlib
from pathlib import Path
import subprocess
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("prepare_nonhuman_d3.sh").resolve()


class DownloadRecovery(unittest.TestCase):
    def run_case(self, *, existing=None, replacement=b"valid fastq", size=None, download_fails=False):
        with tempfile.TemporaryDirectory(prefix="d3-download-test-") as temporary:
            root = Path(temporary)
            (root / "quarantine").mkdir()
            fixture = root / "fixture"
            fixture.write_bytes(replacement)
            output = root / "r1.gz"
            if existing is not None:
                output.write_bytes(existing)
            expected = hashlib.md5(b"valid fastq").hexdigest()
            download = "return 1" if download_fails else '[[ -f "$2" ]] || cp -- "$ROOT/fixture" "$2"'
            result = subprocess.run(
                ["bash", "-c", '''
source "$1"
ROOT=$2
D3_RETRY_WAIT=0
download() { DOWNLOADS=$((DOWNLOADS + 1)); ''' + download + '''; }
DOWNLOADS=0
if ensure_fastq unused "$ROOT/r1.gz" "$3" "$4"; then result=0; else result=$?; fi
printf 'attempts=%s\n' "$DOWNLOADS"
exit "$result"
''', "test", str(SCRIPT), str(root), str(size if size is not None else len(b"valid fastq")), expected],
                text=True, capture_output=True, timeout=10)
            quarantined = list((root / "quarantine").glob("*/r1.gz"))
            return result, [p.read_bytes() for p in quarantined], output.read_bytes() if output.exists() else None

    def test_valid_file_is_preserved(self):
        result, bad, output = self.run_case(existing=b"valid fastq")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("attempts=1", result.stdout)
        self.assertEqual(bad, [])
        self.assertEqual(output, b"valid fastq")

    def test_same_size_corruption_is_quarantined_and_retried(self):
        result, bad, output = self.run_case(existing=b"bad!! fastq")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("attempts=2", result.stdout)
        self.assertEqual(bad, [b"bad!! fastq"])
        self.assertEqual(output, b"valid fastq")

    def test_wrong_size_is_rejected_even_if_md5_matches(self):
        result, bad, output = self.run_case(size=1)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("attempts=3", result.stdout)
        self.assertEqual(len(bad), 3)
        self.assertIsNone(output)

    def test_failed_download_does_not_promote_partial_file(self):
        result, bad, output = self.run_case(download_fails=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("attempts=3", result.stdout)
        self.assertEqual(bad, [])
        self.assertIsNone(output)

    def test_aspera_keeps_aria_partial_separate_and_promotes_only_success(self):
        with tempfile.TemporaryDirectory(prefix="d3-aspera-test-") as temporary:
            root = Path(temporary)
            (root / "key").touch()
            (root / "fixture").write_bytes(b"valid fastq")
            (root / "r1.gz.part").write_bytes(b"old aria progress")
            for success in (False, True):
                mock = ('cp -- "$ROOT/fixture" "${@: -1}"' if success else
                        'printf partial > "${@: -1}"; return 1')
                result = subprocess.run(
                    ["bash", "-c", '''
source "$1"
ROOT=$2
D3_DOWNLOAD_METHOD=aspera
D3_ASPERA_KEY="$ROOT/key"
timeout() { ''' + mock + '''; }
download "$FASTQ_URL/r1.gz" "$ROOT/r1.gz"
''', "test", str(SCRIPT), str(root)], text=True, capture_output=True, timeout=10)
                self.assertEqual((root / "r1.gz.part").read_bytes(), b"old aria progress")
                self.assertEqual(result.returncode == 0, success, result.stderr)
                self.assertEqual((root / "r1.gz").exists(), success)
            self.assertEqual((root / "r1.gz").read_bytes(), b"valid fastq")
            self.assertFalse((root / "r1.gz.aspera.part").exists())


if __name__ == "__main__":
    unittest.main()
