#include "AdvancedDepictor.h"

#include <GraphMol/Conformer.h>
#include <GraphMol/ROMol.h>

#include <algorithm>
#include <limits>
#include <vector>

namespace RDDepict {
namespace {

struct Candidate {
  std::vector<RDGeom::Point3D> positions;
  DepictionScore score;
};

Candidate captureCandidate(const RDKit::ROMol &mol, int confId,
                           const DepictionScoreWeights &weights) {
  Candidate candidate;
  const auto &conf = mol.getConformer(confId);
  candidate.positions = conf.getPositions();
  candidate.score = scoreDepiction(mol, confId, weights);
  return candidate;
}

void applyCandidate(RDKit::ROMol &mol, const Candidate &candidate) {
  mol.clearConformers();
  auto *conf = new RDKit::Conformer(mol.getNumAtoms());
  conf->set3D(false);
  for (unsigned int i = 0; i < candidate.positions.size(); ++i) {
    conf->setAtomPos(i, candidate.positions[i]);
  }
  mol.addConformer(conf, true);
}

}  // namespace

AdvancedDepictionResult compute2DCoordsAdvanced(
    RDKit::ROMol &mol, const AdvancedDepictionParams &params) {
  AdvancedDepictionResult result;
  if (!mol.getNumAtoms()) {
    return result;
  }

  const unsigned int maxCandidates = std::max(1u, params.maxCandidates);
  std::vector<Candidate> candidates;
  candidates.reserve(maxCandidates);

  auto evaluate = [&](const Compute2DCoordParameters &input) {
    if (candidates.size() >= maxCandidates) {
      return;
    }
    RDKit::ROMol copy(mol);
    auto cp = input;
    cp.clearConfs = true;
    const auto confId = compute2DCoords(copy, cp);
    candidates.push_back(captureCandidate(copy, confId, params.scoreWeights));
  };

  Compute2DCoordParameters base;
  base.canonOrient = params.canonOrient;
  base.forceRDKit = true;
  evaluate(base);

  if (params.useRingTemplates && candidates.size() < maxCandidates) {
    auto templated = base;
    templated.useRingTemplates = true;
    evaluate(templated);
  }

  // Reserve one candidate slot for CoordGen when requested.
  const unsigned int sampledTarget =
      params.useCoordGenCandidate && maxCandidates > 1 ? maxCandidates - 1
                                                       : maxCandidates;
  unsigned int candidateIndex = 0;
  while (candidates.size() < sampledTarget) {
    auto sampled = base;
    sampled.useRingTemplates = params.useRingTemplates;
    sampled.nFlipsPerSample = params.flipsPerSample;
    sampled.nSamples = params.randomSamples;
    sampled.sampleSeed = params.seed + static_cast<int>(candidateIndex * 7919u);
    sampled.permuteDeg4Nodes = params.permuteDeg4Nodes;
    evaluate(sampled);
    ++candidateIndex;
  }

  if (params.useCoordGenCandidate && candidates.size() < maxCandidates) {
    const bool previousPreference = preferCoordGen;
    preferCoordGen = true;
    Compute2DCoordParameters coordgen;
    coordgen.canonOrient = params.canonOrient;
    coordgen.forceRDKit = false;
    evaluate(coordgen);
    preferCoordGen = previousPreference;
  }

  const auto best = std::min_element(
      candidates.begin(), candidates.end(),
      [](const Candidate &lhs, const Candidate &rhs) {
        return lhs.score.total < rhs.score.total;
      });
  if (best == candidates.end()) {
    result.score.total = std::numeric_limits<double>::infinity();
    return result;
  }

  applyCandidate(mol, *best);
  result.confId = mol.getConformer().getId();
  result.candidatesEvaluated = candidates.size();
  result.score = best->score;
  return result;
}

}  // namespace RDDepict
