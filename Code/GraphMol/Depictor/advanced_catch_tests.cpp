#include <catch2/catch_all.hpp>

#include <GraphMol/Depictor/AdvancedDepictor.h>
#include <GraphMol/Depictor/DepictionScore.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Conformer.h>

TEST_CASE("depiction scorer detects bond crossings") {
  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("C1CCC1"));
  REQUIRE(mol);

  auto *conf = new RDKit::Conformer(mol->getNumAtoms());
  conf->set3D(false);
  conf->setAtomPos(0, RDGeom::Point3D(-1.0, -1.0, 0.0));
  conf->setAtomPos(1, RDGeom::Point3D(1.0, 1.0, 0.0));
  conf->setAtomPos(2, RDGeom::Point3D(-1.0, 1.0, 0.0));
  conf->setAtomPos(3, RDGeom::Point3D(1.0, -1.0, 0.0));
  mol->addConformer(conf, true);

  const auto score = RDDepict::scoreDepiction(*mol);
  CHECK(score.bondCrossings >= 1);
  CHECK(score.total >= 1000.0);
}

TEST_CASE("advanced depictor evaluates multiple layouts") {
  std::unique_ptr<RDKit::ROMol> mol(
      RDKit::SmilesToMol("CC(C)(C)OC(=O)N1CCC(CC1)c1ccc(Cl)cc1"));
  REQUIRE(mol);

  RDDepict::AdvancedDepictionParams params;
  params.maxCandidates = 6;
  params.randomSamples = 25;
  params.useRingTemplates = true;

  const auto result = RDDepict::compute2DCoordsAdvanced(*mol, params);
  CHECK(result.candidatesEvaluated == 6);
  CHECK(mol->getNumConformers() == 1);
  CHECK(std::isfinite(result.score.total));
}
