#include <RDBoost/python.h>

#include <GraphMol/Depictor/AdvancedDepictor.h>
#include <GraphMol/Depictor/DepictionScore.h>
#include <GraphMol/ROMol.h>

namespace python = boost::python;

namespace {

RDDepict::AdvancedDepictionResult Compute2DCoordsAdvanced(
    RDKit::ROMol &mol, const RDDepict::AdvancedDepictionParams &params) {
  return RDDepict::compute2DCoordsAdvanced(mol, params);
}

RDDepict::DepictionScore ScoreDepiction(
    const RDKit::ROMol &mol, int confId,
    const RDDepict::DepictionScoreWeights &weights) {
  return RDDepict::scoreDepiction(mol, confId, weights);
}

}  // namespace

BOOST_PYTHON_MODULE(rdAdvancedDepictor) {
  python::scope().attr("__doc__") =
      "Experimental high-quality 2D depiction candidate selection and scoring.";

  python::class_<RDDepict::DepictionScoreWeights>("DepictionScoreWeights")
      .def_readwrite("bondCrossing",
                     &RDDepict::DepictionScoreWeights::bondCrossing)
      .def_readwrite("atomAtomClash",
                     &RDDepict::DepictionScoreWeights::atomAtomClash)
      .def_readwrite("atomBondClash",
                     &RDDepict::DepictionScoreWeights::atomBondClash)
      .def_readwrite("congestion", &RDDepict::DepictionScoreWeights::congestion)
      .def_readwrite("angle", &RDDepict::DepictionScoreWeights::angle)
      .def_readwrite("bondLength",
                     &RDDepict::DepictionScoreWeights::bondLength)
      .def_readwrite("aspectRatio",
                     &RDDepict::DepictionScoreWeights::aspectRatio);

  python::class_<RDDepict::DepictionScore>("DepictionScore")
      .def_readonly("total", &RDDepict::DepictionScore::total)
      .def_readonly("bondCrossings", &RDDepict::DepictionScore::bondCrossings)
      .def_readonly("atomAtomClashes",
                    &RDDepict::DepictionScore::atomAtomClashes)
      .def_readonly("atomBondClashes",
                    &RDDepict::DepictionScore::atomBondClashes)
      .def_readonly("congestionPenalty",
                    &RDDepict::DepictionScore::congestionPenalty)
      .def_readonly("anglePenalty", &RDDepict::DepictionScore::anglePenalty)
      .def_readonly("bondLengthPenalty",
                    &RDDepict::DepictionScore::bondLengthPenalty)
      .def_readonly("aspectRatioPenalty",
                    &RDDepict::DepictionScore::aspectRatioPenalty);

  python::class_<RDDepict::AdvancedDepictionParams>("AdvancedDepictionParams")
      .def_readwrite("maxCandidates",
                     &RDDepict::AdvancedDepictionParams::maxCandidates)
      .def_readwrite("randomSamples",
                     &RDDepict::AdvancedDepictionParams::randomSamples)
      .def_readwrite("flipsPerSample",
                     &RDDepict::AdvancedDepictionParams::flipsPerSample)
      .def_readwrite("seed", &RDDepict::AdvancedDepictionParams::seed)
      .def_readwrite("useRingTemplates",
                     &RDDepict::AdvancedDepictionParams::useRingTemplates)
      .def_readwrite("useCoordGenCandidate",
                     &RDDepict::AdvancedDepictionParams::useCoordGenCandidate)
      .def_readwrite("permuteDeg4Nodes",
                     &RDDepict::AdvancedDepictionParams::permuteDeg4Nodes)
      .def_readwrite("canonOrient",
                     &RDDepict::AdvancedDepictionParams::canonOrient)
      .def_readwrite("scoreWeights",
                     &RDDepict::AdvancedDepictionParams::scoreWeights);

  python::class_<RDDepict::AdvancedDepictionResult>("AdvancedDepictionResult")
      .def_readonly("confId", &RDDepict::AdvancedDepictionResult::confId)
      .def_readonly("candidatesEvaluated",
                    &RDDepict::AdvancedDepictionResult::candidatesEvaluated)
      .def_readonly("score", &RDDepict::AdvancedDepictionResult::score);

  python::def("Compute2DCoordsAdvanced", Compute2DCoordsAdvanced,
              (python::arg("mol"),
               python::arg("params") = RDDepict::AdvancedDepictionParams()),
              "Generate several 2D depictions, score them, and keep the best one.");

  python::def("ScoreDepiction", ScoreDepiction,
              (python::arg("mol"), python::arg("confId") = -1,
               python::arg("weights") = RDDepict::DepictionScoreWeights()),
              "Score the current 2D depiction. Lower scores are better.");
}
