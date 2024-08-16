import unittest
from rdkit import Chem
from . import rdShapeAlign
import pathlib

datadir = pathlib.Path(__file__).parent.parent / 'test_data'


class TestCase(unittest.TestCase):

    def setUp(self):
        suppl = Chem.SDMolSupplier(datadir / 'test1.sdf')
        self.ref = suppl[0]
        self.probe = suppl[1]

    def test1_Defaults(self):
        tpl = rdShapeAlign.AlignMol(self.ref, self.probe)
        self.assertAlmostEqual(tpl[0], 0.773, places=3)
        self.assertAlmostEqual(tpl[1], 0.303, places=3)

    def test2_NoColor(self):
        tpl = rdShapeAlign.AlignMol(self.ref, self.probe, useColors=False)
        self.assertAlmostEqual(tpl[0], 0.773, places=3)
        self.assertAlmostEqual(tpl[1], 0.0, places=3)

    def test3_FromShape(self):
        shp = rdShapeAlign.PrepareConformer(self.ref)
        self.assertTrue(type(shp) == rdShapeAlign.ShapeInput)
        tpl = rdShapeAlign.AlignMol(shp, self.probe)
        self.assertAlmostEqual(tpl[0], 0.773, places=3)
        self.assertAlmostEqual(tpl[1], 0.303, places=3)


if __name__ == '__main__':
    unittest.main()
