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
        params.beamWidth = 4
        params.maxTrials = 32
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
        self.assertLessEqual(result.trialsEvaluated, params.maxTrials)

    def test_beam_search_controls_are_exposed(self):
        params = rdAdvancedDepiction.SmartAbbreviationParams()
        params.beamWidth = 7
        params.maxTrials = 55
        self.assertEqual(params.beamWidth, 7)
        self.assertEqual(params.maxTrials, 55)


if __name__ == "__main__":
    unittest.main()
