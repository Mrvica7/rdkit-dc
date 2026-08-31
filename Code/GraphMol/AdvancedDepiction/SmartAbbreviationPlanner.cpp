//
// Experimental high-quality depiction support.
//
#include "SmartAbbreviationPlanner.h"

#include <GraphMol/Abbreviations/Abbreviations.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>

#include <algorithm>
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

std::vector<Abbreviations::AbbreviationDefinition> buildDefinitions(
    const SmartAbbreviationParams &params) {
  auto definitions = getCommercialAbbreviations();
  if (params.useDefaultAbbreviations) {
    auto defaults = Abbreviations::Utils::getDefaultAbbreviations();
    definitions.insert(definitions.end(), defaults.begin(), defaults.end());
  }
  return definitions;
}

struct Trial {
  std::unique_ptr<RWMol> mol;
  RDDepict::DepictionScore score;
  std::string label;
  unsigned int atomsRemoved = 0;
  double objective = std::numeric_limits<double>::infinity();
};

}  // namespace

SmartAbbreviationResult compute2DCoordsSmart(
    RWMol &mol, const SmartAbbreviationParams &params) {
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

  if (!params.maxAbbreviations) {
    mol = current;
    return result;
  }

  const auto definitions = buildDefinitions(params);
  auto currentScore = baselineDepiction.score;
  auto currentObjective = result.baselineObjective;

  for (unsigned int step = 0; step < params.maxAbbreviations; ++step) {
    if (current.getNumAtoms() < params.minAtomsForAutoAbbreviation &&
        !hasHardDepictionDefect(currentScore)) {
      break;
    }

    Trial best;

    // Evaluate each definition independently. RDKit's normal abbreviation
    // function intentionally uses definition order as a greedy priority. Here
    // we want layout quality to determine priority, so each applicable match is
    // trialed against the same current molecule before choosing a winner.
    for (const auto &definition : definitions) {
      if (!labelAllowed(definition.label, params)) {
        continue;
      }
      std::vector<Abbreviations::AbbreviationDefinition> singleDefinition{
          definition};
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          current, singleDefinition, params.maxCoverage);

      for (const auto &match : matches) {
        RWMol trialMol(current);
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

        auto depiction =
            RDDepict::compute2DCoordsAdvanced(trialMol, params.depictionParams);
        ++result.trialsEvaluated;

        const unsigned int totalRemoved =
            result.atomsRemoved + removedThisTrial;
        const unsigned int totalAbbreviations =
            static_cast<unsigned int>(result.abbreviations.size()) + 1;
        const double objective =
            objectiveFor(depiction.score, totalAbbreviations, totalRemoved,
                         params);

        if (objective < best.objective) {
          best.mol = std::make_unique<RWMol>(trialMol);
          best.score = depiction.score;
          best.label = definition.label;
          best.atomsRemoved = removedThisTrial;
          best.objective = objective;
        }
      }
    }

    if (!best.mol ||
        currentObjective - best.objective < params.minimumObjectiveImprovement) {
      break;
    }

    current = *best.mol;
    currentScore = best.score;
    currentObjective = best.objective;
    result.atomsRemoved += best.atomsRemoved;
    result.abbreviations.push_back(best.label);
    result.finalScore = best.score;
    result.finalObjective = best.objective;
  }

  mol = current;
  return result;
}

}  // namespace AdvancedDepiction
}  // namespace RDKit
