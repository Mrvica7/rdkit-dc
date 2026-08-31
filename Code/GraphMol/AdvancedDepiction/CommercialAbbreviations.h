//
// Experimental high-quality depiction support.
//
#pragma once

#include <GraphMol/Abbreviations/Abbreviations.h>
#include <RDGeneral/RDExportMacros.h>

#ifdef RDKIT_ADVANCEDDEPICTION_BUILD
#define RDKIT_ADVANCEDDEPICTION_EXPORT RDKIT_EXPORT_API
#else
#define RDKIT_ADVANCEDDEPICTION_EXPORT RDKIT_IMPORT_API
#endif

namespace RDKit {
namespace AdvancedDepiction {

//! A curated set of abbreviations commonly useful in medicinal chemistry,
//! supplier catalogs, and synthetic schemes. These supplement (rather than
//! replace) RDKit's default abbreviation definitions.
RDKIT_ADVANCEDDEPICTION_EXPORT std::vector<Abbreviations::AbbreviationDefinition>
getCommercialAbbreviations();

}  // namespace AdvancedDepiction
}  // namespace RDKit
