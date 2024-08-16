import unittest
from rdkit import Chem
from . import rdShapeAlign
import pathlib

datadir = pathlib.Path(__file__).parent.parent / 'test_data'


class TestCase(unittest.TestCase):

    def test1_Defaults(self):
        suppl = Chem.SDMolSupplier(datadir / 'test1.sdf')
        ref = suppl[0]
        probe = suppl[1]
        tpl = rdShapeAlign.AlignMol(ref, probe)
        self.assertAlmostEqual(tpl[0], 0.773, places=3)
        self.assertAlmostEqual(tpl[1], 0.303, places=3)

    def test2_NoColor(self):
        suppl = Chem.SDMolSupplier(datadir / 'test1.sdf')
        ref = suppl[0]
        probe = suppl[1]
        tpl = rdShapeAlign.AlignMol(ref, probe, useColors=False)
        self.assertAlmostEqual(tpl[0], 0.773, places=3)
        self.assertAlmostEqual(tpl[1], 0.0, places=3)

    def test3_FromShape(self):
        suppl = Chem.SDMolSupplier(datadir / 'test1.sdf')
        ref = suppl[0]
        probe = suppl[1]
        shp = rdShapeAlign.PrepareConformer(ref)
        self.assertTrue(type(shp) == rdShapeAlign.ShapeInput)
        tpl = rdShapeAlign.AlignMol(shp, probe)
        self.assertAlmostEqual(tpl[0], 0.773, places=3)
        self.assertAlmostEqual(tpl[1], 0.303, places=3)


if __name__ == '__main__':
    unittest.main()
