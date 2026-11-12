#!/usr/bin/env python3
"""Small independent checks of canonical records, multiplicity and DUP identity."""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
HEADER = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000\n'
ROW = 'readA\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tAS:i:4\n'


class DigestTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.directory.name) / 'digest'
        subprocess.run(['c++', '-std=c++11', '-O2', '-Wall', '-Wextra', '-Werror',
                        str(ROOT/'tools/sam_correctness_digest.cpp'), '-lcrypto', '-o', str(cls.binary)], check=True)
        cls.portable = Path(cls.directory.name) / 'digest-portable'
        subprocess.run(['c++', '-std=c++11', '-O2', '-Wall', '-Wextra', '-Werror',
                        '-DSWBAM_DIGEST_PORTABLE_SHA256', str(ROOT/'tools/sam_correctness_digest.cpp'),
                        '-o', str(cls.portable)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.directory.cleanup()

    def digest(self, text):
        run = subprocess.run([str(self.binary)], input=text, universal_newlines=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return json.loads(run.stdout)

    def test_canonical_fields_and_header(self):
        actual = self.digest(HEADER+'@PG\tID:ignore\n'+ROW)
        canonical = ROW.replace('NM:i:0\tAS:i:4', 'AS:i:4\tNM:i:0').encode()
        self.assertEqual(actual['ordered_sha256'], hashlib.sha256(canonical).hexdigest())
        self.assertEqual(actual['header_without_pg_sha256'], hashlib.sha256(HEADER.encode()).hexdigest())
        self.assertEqual(actual['records'], 1)

    def test_order_and_aux_order(self):
        other = ROW.replace('readA', 'readB').replace('\t10\t', '\t20\t')
        a = self.digest(HEADER+ROW+other)
        b = self.digest(HEADER+other+ROW.replace('NM:i:0\tAS:i:4', 'AS:i:4\tNM:i:0'))
        self.assertEqual(a['record_multiset_sha256_sum'], b['record_multiset_sha256_sum'])
        self.assertNotEqual(a['ordered_sha256'], b['ordered_sha256'])
        self.assertEqual(b['coordinate_errors'], 1)
        self.assertNotEqual(a['record_multiset_sha256_sum'], self.digest(HEADER+ROW+ROW)['record_multiset_sha256_sum'])

    def test_duplicate_identity_and_neutral_flags(self):
        marked = ROW.replace('readA\t0\t', 'readA\t1024\t')
        plain, dup = self.digest(HEADER+ROW), self.digest(HEADER+marked)
        self.assertEqual(plain['without_dup_multiset_sha256_sum'], dup['without_dup_multiset_sha256_sum'])
        self.assertNotEqual(plain['record_multiset_sha256_sum'], dup['record_multiset_sha256_sum'])
        self.assertEqual(dup['duplicates'], 1)
        self.assertNotEqual(dup['duplicate_members_sha256_sum'], self.digest(HEADER+marked.replace('readA','readB'))['duplicate_members_sha256_sum'])

    def test_fail_closed_on_truncated_fields(self):
        run = subprocess.run([str(self.binary)], input=HEADER+'read\t0\tchr1\n', universal_newlines=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertNotEqual(run.returncode, 0)
        self.assertEqual(run.stdout, '')

    def test_portable_sha_matches_openssl(self):
        for text in (HEADER, HEADER+ROW, HEADER+ROW*1000,
                     HEADER+ROW.replace('readA\t0\t', 'readA\t1024\t'),
                     HEADER+ROW.replace('readA', 'long_qname_'+'x'*200)):
            actual = subprocess.run([str(self.portable)], input=text, universal_newlines=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            self.assertEqual(json.loads(actual.stdout), self.digest(text))

    def test_ranges_require_bam_input(self):
        run = subprocess.run([str(self.portable),'--range-end','100'],input=HEADER,
                             universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        self.assertNotEqual(run.returncode,0)
        self.assertEqual(run.stdout,'')


if __name__ == '__main__':
    unittest.main()
