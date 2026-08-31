#include <catch2/catch_all.hpp>

#include <GraphMol/AdvancedDepiction/CommercialAbbreviations.h>
#include <GraphMol/AdvancedDepiction/SmartAbbreviationPlanner.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <algorithm>
#include <memory>

TEST_CASE("commercial abbreviation library includes medicinal chemistry groups") {
  const auto defs = RDKit::AdvancedDepiction::getCommercialAbbreviations();
  auto hasLabel = [&defs](const std::string &label) {
    return std::any_of(defs.begin(), defs.end(), [&label](const auto &def) {
      return def.label == label;
    });
  };

  CHECK(hasLabel("Boc"));
  CHECK(hasLabel("Cbz"));
  CHECK(hasLabel("Fmoc"));
  CHECK(hasLabel("Bpin"));
  CHECK(hasLabel("TBS"));
}

TEST_CASE("smart planner can choose Boc when simplification has clear utility") {
  std::unique_ptr<RDKit::ROMol> parsed(
      RDKit::SmilesToMol("CC(C)(C)OC(=O)N1CCC(CC1)c1ccc(Cl)cc1"));
  REQUIRE(parsed);
  RDKit::RWMol mol(*parsed);
  const auto originalAtoms = mol.getNumAtoms();

  RDKit::AdvancedDepiction::SmartAbbreviationParams params;
  params.maxAbbreviations = 1;
  params.minAtomsForAutoAbbreviation = 0;
  params.minAtomsRemoved = 3;
  params.atomReductionReward = 100.0;
  params.abbreviationPenalty = 0.0;
  params.minimumObjectiveImprovement = 0.0;
  params.depictionParams.maxCandidates = 4;
  params.depictionParams.randomSamples = 10;
  params.depictionParams.useCoordGenCandidate = false;

  const auto result =
      RDKit::AdvancedDepiction::compute2DCoordsSmart(mol, params);

  REQUIRE(result.abbreviations.size() == 1);
  CHECK(result.abbreviations.front() == "Boc");
  CHECK(result.atomsRemoved >= 3);
  CHECK(mol.getNumAtoms() < originalAtoms);
  CHECK(mol.getNumConformers() == 1);
}

TEST_CASE("smart planner can leave abbreviation disabled") {
  std::unique_ptr<RDKit::ROMol> parsed(
      RDKit::SmilesToMol("CC(C)(C)OC(=O)N1CCC(CC1)c1ccc(Cl)cc1"));
  REQUIRE(parsed);
  RDKit::RWMol mol(*parsed);

  RDKit::AdvancedDepiction::SmartAbbreviationParams params;
  params.maxAbbreviations = 0;
  params.depictionParams.maxCandidates = 3;
  params.depictionParams.useCoordGenCandidate = false;

  const auto result =
      RDKit::AdvancedDepiction::compute2DCoordsSmart(mol, params);
  CHECK(result.abbreviations.empty());
  CHECK(result.atomsRemoved == 0);
  CHECK(mol.getNumConformers() == 1);
}
