import math
import unittest

from rdkit import Chem
from rdkit.Chem import rdAdvancedDepictor


class TestAdvancedDepictor(unittest.TestCase):
    def test_compute_advanced_coords(self):
        mol = Chem.MolFromSmiles("CC(C)(C)OC(=O)N1CCC(CC1)c1ccc(Cl)cc1")
        params = rdAdvancedDepictor.AdvancedDepictionParams()
        params.maxCandidates = 5
        params.randomSamples = 20
        result = rdAdvancedDepictor.Compute2DCoordsAdvanced(mol, params)

        self.assertEqual(result.candidatesEvaluated, 5)
        self.assertEqual(mol.GetNumConformers(), 1)
        self.assertTrue(math.isfinite(result.score.total))

    def test_score_current_depiction(self):
        mol = Chem.MolFromSmiles("c1ccccc1C(=O)NCC")
        params = rdAdvancedDepictor.AdvancedDepictionParams()
        rdAdvancedDepictor.Compute2DCoordsAdvanced(mol, params)
        score = rdAdvancedDepictor.ScoreDepiction(mol)
        self.assertTrue(math.isfinite(score.total))


if __name__ == "__main__":
    unittest.main()
