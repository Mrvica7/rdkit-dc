//
// Experimental high-quality depiction support.
//
#include "SmartAbbreviationPlanner.h"

#include <GraphMol/Abbreviations/Abbreviations.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>

namespace RDKit {
namespace AdvancedDepiction {
namespace {

bool hasHardDepictionDefect(const RDDepict::DepictionScore &score) {
  return score.bondCrossings || score.atomAtomClashes || score.atomBondClashes;
}

double objectiveFor(const RDDepict::DepictionScore &score,
                    unsigned int abbreviationCount,
                    unsigned int atomsRemoved,
                    const SmartAbbreviationParams &params) {
  return score.total +
         params.abbreviationPenalty * static_cast<double>(abbreviationCount) -
         params.atomReductionReward * static_cast<double>(atomsRemoved);
}

bool labelAllowed(const std::string &label,
                  const SmartAbbreviationParams &params) {
  if (!params.allowPhenylAbbreviation && label == "Ph") {
    return false;
  }
  return true;
}

void refreshRingInfo(RWMol &mol) {
  mol.updatePropertyCache(false);
  mol.getRingInfo()->reset();
  MolOps::symmetrizeSSSR(mol);
}

void moveBackToInput(ROMol &destination, RWMol &source) {
  // ROMol deliberately disables copy assignment, but supports move assignment.
  // RWMol adds no independent topology storage, so moving its ROMol base back
  // into a normal Chem.Mol safely replaces the optimized topology/conformers.
  destination = std::move(static_cast<ROMol &>(source));
}

std::vector<Abbreviations::AbbreviationDefinition> buildDefinitions(
    const SmartAbbreviationParams &params) {
  auto definitions = getCommercialAbbreviations();
  if (params.useDefaultAbbreviations) {
    auto defaults = Abbreviations::Utils::getDefaultAbbreviations();
    definitions.insert(definitions.end(), defaults.begin(), defaults.end());
  }
  return definitions;
}

struct SearchState {
  std::unique_ptr<RWMol> mol;
  RDDepict::DepictionScore score;
  std::vector<std::string> labels;
  unsigned int atomsRemoved = 0;
  double objective = std::numeric_limits<double>::infinity();
};

bool stateLess(const SearchState &lhs, const SearchState &rhs) {
  constexpr double epsilon = 1e-9;
  if (std::abs(lhs.objective - rhs.objective) > epsilon) {
    return lhs.objective < rhs.objective;
  }
  if (lhs.labels.size() != rhs.labels.size()) {
    return lhs.labels.size() < rhs.labels.size();
  }
  if (lhs.atomsRemoved != rhs.atomsRemoved) {
    return lhs.atomsRemoved > rhs.atomsRemoved;
  }
  return lhs.labels < rhs.labels;
}

}  // namespace

SmartAbbreviationResult compute2DCoordsSmart(
    ROMol &mol, const SmartAbbreviationParams &params) {
  SmartAbbreviationResult result;
  if (!mol.getNumAtoms()) {
    return result;
  }

  RWMol current(mol);
  auto baselineDepiction =
      RDDepict::compute2DCoordsAdvanced(current, params.depictionParams);
  result.baselineScore = baselineDepiction.score;
  result.finalScore = baselineDepiction.score;
  result.baselineObjective = baselineDepiction.score.total;
  result.finalObjective = result.baselineObjective;

  if (!params.maxAbbreviations || !params.beamWidth || !params.maxTrials) {
    moveBackToInput(mol, current);
    return result;
  }

  if (current.getNumAtoms() < params.minAtomsForAutoAbbreviation &&
      !hasHardDepictionDefect(baselineDepiction.score)) {
    moveBackToInput(mol, current);
    return result;
  }

  const auto definitions = buildDefinitions(params);

  SearchState initial;
  initial.mol = std::make_unique<RWMol>(current);
  initial.score = baselineDepiction.score;
  initial.objective = result.baselineObjective;

  std::vector<SearchState> frontier;
  frontier.push_back(std::move(initial));

  std::unique_ptr<RWMol> bestMol;
  RDDepict::DepictionScore bestScore = baselineDepiction.score;
  std::vector<std::string> bestLabels;
  unsigned int bestAtomsRemoved = 0;
  double bestObjective = std::numeric_limits<double>::infinity();

  for (unsigned int depth = 0;
       depth < params.maxAbbreviations && !frontier.empty() &&
       result.trialsEvaluated < params.maxTrials;
       ++depth) {
    std::vector<SearchState> nextFrontier;

    for (const auto &parent : frontier) {
      if (result.trialsEvaluated >= params.maxTrials) {
        break;
      }
      if (parent.mol->getNumAtoms() < params.minAtomsForAutoAbbreviation &&
          !hasHardDepictionDefect(parent.score)) {
        continue;
      }

      // Every applicable match is trialed against the same parent state. This
      // handles overlapping definitions without relying on definition order.
      // Keeping several parent states for the next depth also allows useful
      // combinations that a purely greedy first choice would miss.
      for (const auto &definition : definitions) {
        if (result.trialsEvaluated >= params.maxTrials) {
          break;
        }
        if (!labelAllowed(definition.label, params)) {
          continue;
        }

        std::vector<Abbreviations::AbbreviationDefinition> singleDefinition{
            definition};
        auto matches = Abbreviations::findApplicableAbbreviationMatches(
            *parent.mol, singleDefinition, params.maxCoverage);

        for (const auto &match : matches) {
          if (result.trialsEvaluated >= params.maxTrials) {
            break;
          }

          RWMol trialMol(*parent.mol);
          const auto beforeAtoms = trialMol.getNumAtoms();
          Abbreviations::applyMatches(trialMol, {match});
          refreshRingInfo(trialMol);
          const auto afterAtoms = trialMol.getNumAtoms();
          if (afterAtoms >= beforeAtoms) {
            continue;
          }

          const auto removedThisTrial = beforeAtoms - afterAtoms;
          if (removedThisTrial < params.minAtomsRemoved) {
            continue;
          }

          auto depiction = RDDepict::compute2DCoordsAdvanced(
              trialMol, params.depictionParams);
          ++result.trialsEvaluated;

          SearchState candidate;
          candidate.mol = std::make_unique<RWMol>(trialMol);
          candidate.score = depiction.score;
          candidate.labels = parent.labels;
          candidate.labels.push_back(definition.label);
          candidate.atomsRemoved = parent.atomsRemoved + removedThisTrial;
          candidate.objective =
              objectiveFor(candidate.score,
                           static_cast<unsigned int>(candidate.labels.size()),
                           candidate.atomsRemoved, params);

          if (!bestMol || candidate.objective < bestObjective) {
            bestMol = std::make_unique<RWMol>(*candidate.mol);
            bestScore = candidate.score;
            bestLabels = candidate.labels;
            bestAtomsRemoved = candidate.atomsRemoved;
            bestObjective = candidate.objective;
          }

          nextFrontier.push_back(std::move(candidate));
        }
      }
    }

    if (nextFrontier.empty()) {
      break;
    }

    std::sort(nextFrontier.begin(), nextFrontier.end(), stateLess);
    if (nextFrontier.size() > params.beamWidth) {
      nextFrontier.resize(params.beamWidth);
    }
    frontier = std::move(nextFrontier);
  }

  if (bestMol &&
      result.baselineObjective - bestObjective >=
          params.minimumObjectiveImprovement) {
    current = *bestMol;
    result.finalScore = bestScore;
    result.finalObjective = bestObjective;
    result.atomsRemoved = bestAtomsRemoved;
    result.abbreviations = std::move(bestLabels);
  }

  moveBackToInput(mol, current);
  return result;
}

}  // namespace AdvancedDepiction
}  // namespace RDKit
