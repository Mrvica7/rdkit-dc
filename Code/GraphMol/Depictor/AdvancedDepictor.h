//
//  Experimental high-quality 2D depiction support.
//
#ifndef RDKIT_ADVANCED_DEPICTOR_H
#define RDKIT_ADVANCED_DEPICTOR_H

#include <GraphMol/Depictor/DepictionScore.h>

namespace RDKit {
class ROMol;
}

namespace RDDepict {

struct RDKIT_DEPICTOR_EXPORT AdvancedDepictionParams {
  unsigned int maxCandidates = 16;
  unsigned int randomSamples = 100;
  unsigned int flipsPerSample = 3;
  int seed = 0xC0FFEE;
  bool useRingTemplates = true;
  bool useCoordGenCandidate = true;
  bool permuteDeg4Nodes = true;
  bool canonOrient = true;
  DepictionScoreWeights scoreWeights;
};

struct RDKIT_DEPICTOR_EXPORT AdvancedDepictionResult {
  unsigned int confId = 0;
  unsigned int candidatesEvaluated = 0;
  DepictionScore score;
};

RDKIT_DEPICTOR_EXPORT AdvancedDepictionResult compute2DCoordsAdvanced(
    RDKit::ROMol &mol,
    const AdvancedDepictionParams &params = AdvancedDepictionParams());

}  // namespace RDDepict

#endif
