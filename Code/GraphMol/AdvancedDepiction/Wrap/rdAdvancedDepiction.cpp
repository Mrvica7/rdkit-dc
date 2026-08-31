#include <RDBoost/python.h>

#include <GraphMol/AdvancedDepiction/CommercialAbbreviations.h>
#include <GraphMol/AdvancedDepiction/SmartAbbreviationPlanner.h>
#include <GraphMol/RWMol.h>

namespace python = boost::python;

namespace {

python::list getAbbreviationLabels(
    const RDKit::AdvancedDepiction::SmartAbbreviationResult &result) {
  python::list labels;
  for (const auto &label : result.abbreviations) {
    labels.append(label);
  }
  return labels;
}

python::list getCommercialAbbreviationLabels() {
  python::list labels;
  for (const auto &definition :
       RDKit::AdvancedDepiction::getCommercialAbbreviations()) {
    labels.append(definition.label);
  }
  return labels;
}

unsigned int getMaxCandidates(
    const RDKit::AdvancedDepiction::SmartAbbreviationParams &params) {
  return params.depictionParams.maxCandidates;
}
void setMaxCandidates(
    RDKit::AdvancedDepiction::SmartAbbreviationParams &params,
    unsigned int value) {
  params.depictionParams.maxCandidates = value;
}

unsigned int getRandomSamples(
    const RDKit::AdvancedDepiction::SmartAbbreviationParams &params) {
  return params.depictionParams.randomSamples;
}
void setRandomSamples(
    RDKit::AdvancedDepiction::SmartAbbreviationParams &params,
    unsigned int value) {
  params.depictionParams.randomSamples = value;
}

bool getUseCoordGenCandidate(
    const RDKit::AdvancedDepiction::SmartAbbreviationParams &params) {
  return params.depictionParams.useCoordGenCandidate;
}
void setUseCoordGenCandidate(
    RDKit::AdvancedDepiction::SmartAbbreviationParams &params, bool value) {
  params.depictionParams.useCoordGenCandidate = value;
}

bool getUseRingTemplates(
    const RDKit::AdvancedDepiction::SmartAbbreviationParams &params) {
  return params.depictionParams.useRingTemplates;
}
void setUseRingTemplates(
    RDKit::AdvancedDepiction::SmartAbbreviationParams &params, bool value) {
  params.depictionParams.useRingTemplates = value;
}

double getBaselineScore(
    const RDKit::AdvancedDepiction::SmartAbbreviationResult &result) {
  return result.baselineScore.total;
}
double getFinalScore(
    const RDKit::AdvancedDepiction::SmartAbbreviationResult &result) {
  return result.finalScore.total;
}
unsigned int getFinalBondCrossings(
    const RDKit::AdvancedDepiction::SmartAbbreviationResult &result) {
  return result.finalScore.bondCrossings;
}
unsigned int getFinalAtomAtomClashes(
    const RDKit::AdvancedDepiction::SmartAbbreviationResult &result) {
  return result.finalScore.atomAtomClashes;
}
unsigned int getFinalAtomBondClashes(
    const RDKit::AdvancedDepiction::SmartAbbreviationResult &result) {
  return result.finalScore.atomBondClashes;
}

}  // namespace

BOOST_PYTHON_MODULE(rdAdvancedDepiction) {
  python::scope().attr("__doc__") =
      "Layout-aware high-quality depiction and smart abbreviation support.";

  python::class_<RDKit::AdvancedDepiction::SmartAbbreviationParams>(
      "SmartAbbreviationParams")
      .def_readwrite("maxCoverage",
                     &RDKit::AdvancedDepiction::SmartAbbreviationParams::maxCoverage)
      .def_readwrite(
          "maxAbbreviations",
          &RDKit::AdvancedDepiction::SmartAbbreviationParams::maxAbbreviations)
      .def_readwrite("beamWidth",
                     &RDKit::AdvancedDepiction::SmartAbbreviationParams::beamWidth)
      .def_readwrite("maxTrials",
                     &RDKit::AdvancedDepiction::SmartAbbreviationParams::maxTrials)
      .def_readwrite(
          "minAtomsForAutoAbbreviation",
          &RDKit::AdvancedDepiction::SmartAbbreviationParams::
              minAtomsForAutoAbbreviation)
      .def_readwrite(
          "minAtomsRemoved",
          &RDKit::AdvancedDepiction::SmartAbbreviationParams::minAtomsRemoved)
      .def_readwrite(
          "atomReductionReward",
          &RDKit::AdvancedDepiction::SmartAbbreviationParams::atomReductionReward)
      .def_readwrite(
          "abbreviationPenalty",
          &RDKit::AdvancedDepiction::SmartAbbreviationParams::
              abbreviationPenalty)
      .def_readwrite(
          "minimumObjectiveImprovement",
          &RDKit::AdvancedDepiction::SmartAbbreviationParams::
              minimumObjectiveImprovement)
      .def_readwrite(
          "allowPhenylAbbreviation",
          &RDKit::AdvancedDepiction::SmartAbbreviationParams::
              allowPhenylAbbreviation)
      .def_readwrite(
          "useDefaultAbbreviations",
          &RDKit::AdvancedDepiction::SmartAbbreviationParams::
              useDefaultAbbreviations)
      .add_property("maxCandidates", &getMaxCandidates, &setMaxCandidates)
      .add_property("randomSamples", &getRandomSamples, &setRandomSamples)
      .add_property("useCoordGenCandidate", &getUseCoordGenCandidate,
                    &setUseCoordGenCandidate)
      .add_property("useRingTemplates", &getUseRingTemplates,
                    &setUseRingTemplates);

  python::class_<RDKit::AdvancedDepiction::SmartAbbreviationResult>(
      "SmartAbbreviationResult", python::no_init)
      .def_readonly("baselineObjective",
                    &RDKit::AdvancedDepiction::SmartAbbreviationResult::
                        baselineObjective)
      .def_readonly("finalObjective",
                    &RDKit::AdvancedDepiction::SmartAbbreviationResult::
                        finalObjective)
      .def_readonly("atomsRemoved",
                    &RDKit::AdvancedDepiction::SmartAbbreviationResult::
                        atomsRemoved)
      .def_readonly("trialsEvaluated",
                    &RDKit::AdvancedDepiction::SmartAbbreviationResult::
                        trialsEvaluated)
      .add_property("baselineScore", &getBaselineScore)
      .add_property("finalScore", &getFinalScore)
      .add_property("finalBondCrossings", &getFinalBondCrossings)
      .add_property("finalAtomAtomClashes", &getFinalAtomAtomClashes)
      .add_property("finalAtomBondClashes", &getFinalAtomBondClashes)
      .add_property("abbreviations", &getAbbreviationLabels);

  python::def("GetCommercialAbbreviationLabels", &getCommercialAbbreviationLabels,
              "Return the curated medicinal-chemistry abbreviation labels.");

  python::def(
      "Compute2DCoordsSmart",
      &RDKit::AdvancedDepiction::compute2DCoordsSmart,
      (python::arg("mol"),
       python::arg("params") =
           RDKit::AdvancedDepiction::SmartAbbreviationParams()),
      "Generate a high-quality 2D depiction and automatically select useful "
      "abbreviations using a bounded layout-aware beam search. The molecule is "
      "modified in place.");
}
