//
// Experimental high-quality depiction support.
//
#pragma once

#include "CommercialAbbreviations.h"
#include <GraphMol/Depictor/AdvancedDepictor.h>

#include <string>
#include <vector>

namespace RDKit {
class RWMol;

namespace AdvancedDepiction {

struct RDKIT_ADVANCEDDEPICTION_EXPORT SmartAbbreviationParams {
  //! Maximum fraction of the current molecule a single abbreviation may cover.
  double maxCoverage = 0.45;

  //! Maximum number of abbreviations introduced automatically.
  unsigned int maxAbbreviations = 3;

  //! Do not automatically abbreviate very small molecules unless they have a
  //! hard depiction defect (crossing/clash). This avoids turning ordinary
  //! catalog structures into strings of labels.
  unsigned int minAtomsForAutoAbbreviation = 14;

  //! A candidate must remove at least this many atoms from explicit depiction.
  unsigned int minAtomsRemoved = 3;

  //! Lower objective is better. Atom reduction earns a reward, while every
  //! abbreviation pays a fixed readability/information penalty.
  double atomReductionReward = 7.0;
  double abbreviationPenalty = 14.0;
  double minimumObjectiveImprovement = 8.0;

  //! Generic phenyl abbreviation is intentionally disabled by default. Named
  //! protecting groups containing aryl systems (Cbz, Fmoc, Bn, PMB, Tr, etc.)
  //! remain eligible.
  bool allowPhenylAbbreviation = false;

  //! RDKit defaults contain useful small groups (Et, tBu, OMe, CF3...). They are
  //! opt-in because aggressive automatic use can reduce structural readability.
  bool useDefaultAbbreviations = false;

  RDDepict::AdvancedDepictionParams depictionParams;
};

struct RDKIT_ADVANCEDDEPICTION_EXPORT SmartAbbreviationResult {
  RDDepict::DepictionScore baselineScore;
  RDDepict::DepictionScore finalScore;
  double baselineObjective = 0.0;
  double finalObjective = 0.0;
  unsigned int atomsRemoved = 0;
  unsigned int trialsEvaluated = 0;
  std::vector<std::string> abbreviations;
};

//! Generate a high-quality 2D depiction while choosing abbreviation matches by
//! layout-aware objective optimization. The molecule is modified in place and
//! contains the selected condensed abbreviation atoms and final 2D coordinates.
RDKIT_ADVANCEDDEPICTION_EXPORT SmartAbbreviationResult
compute2DCoordsSmart(RWMol &mol,
                     const SmartAbbreviationParams &params =
                         SmartAbbreviationParams());

}  // namespace AdvancedDepiction
}  // namespace RDKit
