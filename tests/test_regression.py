import csv, unittest
from pathlib import Path
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[1]
class TestPackagedRegression(unittest.TestCase):
    def rows(self,name):
        with (ROOT/'data'/name).open(newline='',encoding='utf-8-sig') as f:
            return list(csv.DictReader(f))
    def test_counts_and_ids(self):
        hits=self.rows('sancdb_hits.csv'); acc=self.rows('sancdb_hits_accessible.csv')
        self.assertEqual(len(hits),267); self.assertEqual(len(acc),185)
        self.assertEqual(len({r['id'] for r in hits}),267)
        self.assertEqual(len({r['id'] for r in acc}),185)
    def test_processed_accessible_is_subset(self):
        hits={r['id'] for r in self.rows('sancdb_hits.csv')}
        acc={r['id'] for r in self.rows('sancdb_hits_accessible.csv')}
        self.assertTrue(acc <= hits)
    def test_derived_hits_sdf(self):
        ids=[]
        for m in Chem.SDMolSupplier(str(ROOT/'data'/'sancdb_hits.sdf')):
            if m is not None: ids.append(m.GetProp('_Name'))
        self.assertEqual(len(ids),267)
        self.assertEqual(len(set(ids)),267)
if __name__=='__main__': unittest.main()
