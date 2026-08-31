//
//  Experimental high-quality 2D depiction support.
//
#ifndef RDKIT_DEPICTION_SCORE_H
#define RDKIT_DEPICTION_SCORE_H

#include <GraphMol/Depictor/RDDepictor.h>

namespace RDKit {
class ROMol;
}

namespace RDDepict {

struct RDKIT_DEPICTOR_EXPORT DepictionScore {
  double total = 0.0;
  unsigned int bondCrossings = 0;
  unsigned int atomAtomClashes = 0;
  unsigned int atomBondClashes = 0;
  double congestionPenalty = 0.0;
  double anglePenalty = 0.0;
  double bondLengthPenalty = 0.0;
  double aspectRatioPenalty = 0.0;
};

struct RDKIT_DEPICTOR_EXPORT DepictionScoreWeights {
  double bondCrossing = 1000.0;
  double atomAtomClash = 500.0;
  double atomBondClash = 700.0;
  double congestion = 30.0;
  double angle = 25.0;
  double bondLength = 15.0;
  double aspectRatio = 10.0;
};

RDKIT_DEPICTOR_EXPORT DepictionScore scoreDepiction(
    const RDKit::ROMol &mol, int confId = -1,
    const DepictionScoreWeights &weights = DepictionScoreWeights());

}  // namespace RDDepict

#endif
