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
  CHECK(hasLabel("Ph"));
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

TEST_CASE("generic phenyl abbreviation is opt-in") {
  // Biphenyl isolates the generic Ph rule: unlike the previous alkylbenzene
  // fixture, it cannot also be interpreted as Bn or another named protecting
  // group from the commercial abbreviation set.
  const std::string smiles = "c1ccccc1-c1ccccc1";

  std::unique_ptr<RDKit::ROMol> parsed1(RDKit::SmilesToMol(smiles));
  REQUIRE(parsed1);
  RDKit::RWMol conservative(*parsed1);

  RDKit::AdvancedDepiction::SmartAbbreviationParams params;
  params.maxAbbreviations = 1;
  params.maxCoverage = 0.80;
  params.minAtomsForAutoAbbreviation = 0;
  params.minAtomsRemoved = 3;
  params.atomReductionReward = 100.0;
  params.abbreviationPenalty = 0.0;
  params.minimumObjectiveImprovement = 0.0;
  params.depictionParams.maxCandidates = 3;
  params.depictionParams.useCoordGenCandidate = false;

  auto conservativeResult =
      RDKit::AdvancedDepiction::compute2DCoordsSmart(conservative, params);
  CHECK(conservativeResult.abbreviations.empty());

  std::unique_ptr<RDKit::ROMol> parsed2(RDKit::SmilesToMol(smiles));
  REQUIRE(parsed2);
  RDKit::RWMol aggressive(*parsed2);
  params.allowPhenylAbbreviation = true;

  auto aggressiveResult =
      RDKit::AdvancedDepiction::compute2DCoordsSmart(aggressive, params);
  REQUIRE(aggressiveResult.abbreviations.size() == 1);
  CHECK(aggressiveResult.abbreviations.front() == "Ph");
}
