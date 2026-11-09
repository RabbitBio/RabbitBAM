#!/usr/bin/env python3
import unittest
from validate_same_input import plan


class PlanTest(unittest.TestCase):
    def test_scope_and_shared_inputs(self):
        for case in ('D1','D2','D3'):
            samples=plan(case,'/bin/app','/initial.bam','/common','/out')
            self.assertEqual(len(samples),8 if case=='D3' else 7)
            by_label={item['label']:item for item in samples}
            for app,source in [('fixmate','collate'),('sort','fixmate'),('markdup','sort')]:
                self.assertEqual(by_label['standalone_'+app]['source'],'/common/'+source+'.bam')
            self.assertEqual(by_label['shared_collate']['source'],'/initial.bam')
            self.assertEqual(by_label['workflow_fixmate']['source'],by_label['shared_collate']['target'])
            self.assertEqual(by_label['workflow_sort']['source'],by_label['workflow_fixmate']['target'])
            self.assertEqual(by_label['workflow_markdup']['source'],by_label['workflow_sort']['target'])
            self.assertTrue(all('--no-dump' not in s['argv'] for s in samples))
            if case=='D3':
                self.assertEqual(by_label['integrated_sam2bam']['app'],'sam2bam')
                self.assertEqual(by_label['integrated_sam2bam']['capacity'],'40G')

    def test_d2_external_and_latest_stream_mode(self):
        for item in plan('D2','/latest','/initial','/common','/out'):
            self.assertEqual(item['argv'][0],'/latest')
            if item['app'] in ('sort','collate'):
                self.assertEqual(item['check'],'external')
                self.assertEqual(item['argv'][item['argv'].index('-m')+1],'8G')
                self.assertEqual(item['capacity'],'64G')
            if item['app']=='markdup':
                self.assertIn('--stream',item['argv'])
                self.assertEqual(item['argv'][item['argv'].index('-l')+1],'300')


if __name__=='__main__': unittest.main()
