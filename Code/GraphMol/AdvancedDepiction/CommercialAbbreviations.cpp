//
// Experimental high-quality depiction support.
//
#include "CommercialAbbreviations.h"

namespace RDKit {
namespace AdvancedDepiction {
namespace {

// The first atom in each SMARTS is the atom connected to the rest of the
// molecule. parseAbbreviations() supplies the external dummy attachment atom.
// Keep this list conservative: these are groups for which abbreviated forms are
// widely recognizable in medicinal/synthetic chemistry. Generic Ph is present
// but the smart planner keeps it disabled by default.
const std::string commercialAbbreviations = R"ABBREVS(
Boc C(=O)OC(C)(C)C Boc
Cbz C(=O)OCc1ccccc1 Cbz
Fmoc C(=O)OCC1c2ccccc2-c2ccccc21 Fmoc
Alloc C(=O)OCC=C Alloc
Bn Cc1ccccc1 Bn
Ph c1ccccc1 Ph
PMB Cc1ccc(OC)cc1 PMB
Bz C(=O)c1ccccc1 Bz
Piv C(=O)C(C)(C)C Piv
Ts S(=O)(=O)c1ccc(C)cc1 Ts
Ms S(=O)(=O)C Ms
TBS [Si](C)(C)C(C)(C)C TBS
TIPS [Si](C(C)C)(C(C)C)C(C)C TIPS
Tr C(c1ccccc1)(c1ccccc1)c1ccccc1 Tr
Bpin B1OC(C)(C)C(C)(C)O1 Bpin
B(OH)2 [B]([OH])[OH] B(OH)<sub>2</sub>
THP C1OCCCC1 THP
MOM COC MOM
SEM COCC[Si](C)(C)C SEM
Teoc C(=O)OCC[Si](C)(C)C Teoc
)ABBREVS";

}  // namespace

std::vector<Abbreviations::AbbreviationDefinition>
getCommercialAbbreviations() {
  static const auto defs =
      Abbreviations::Utils::parseAbbreviations(commercialAbbreviations);
  return defs;
}

}  // namespace AdvancedDepiction
}  // namespace RDKit
