import unittest

from rdkit import Chem
from rdkit.Chem import rdAdvancedDepiction


class TestAdvancedDepiction(unittest.TestCase):
    def test_commercial_abbreviations(self):
        labels = list(rdAdvancedDepiction.GetCommercialAbbreviationLabels())
        self.assertIn("Boc", labels)
        self.assertIn("Fmoc", labels)
        self.assertIn("Bpin", labels)

    def test_smart_boc_abbreviation(self):
        mol = Chem.MolFromSmiles("CC(C)(C)OC(=O)N1CCC(CC1)c1ccc(Cl)cc1")
        original_atoms = mol.GetNumAtoms()

        params = rdAdvancedDepiction.SmartAbbreviationParams()
        params.maxAbbreviations = 1
        params.minAtomsForAutoAbbreviation = 0
        params.atomReductionReward = 100.0
        params.abbreviationPenalty = 0.0
        params.minimumObjectiveImprovement = 0.0
        params.maxCandidates = 4
        params.randomSamples = 10
        params.useCoordGenCandidate = False

        result = rdAdvancedDepiction.Compute2DCoordsSmart(mol, params)

        self.assertEqual(list(result.abbreviations), ["Boc"])
        self.assertGreater(result.atomsRemoved, 0)
        self.assertLess(mol.GetNumAtoms(), original_atoms)
        self.assertEqual(mol.GetNumConformers(), 1)


if __name__ == "__main__":
    unittest.main()
