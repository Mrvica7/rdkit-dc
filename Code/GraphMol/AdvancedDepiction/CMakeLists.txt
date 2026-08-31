rdkit_library(AdvancedDepiction
              CommercialAbbreviations.cpp
              SmartAbbreviationPlanner.cpp
              LINK_LIBRARIES Depictor Abbreviations SmilesParse SubstructMatch GraphMol RDGeneral)
target_compile_definitions(AdvancedDepiction PRIVATE RDKIT_ADVANCEDDEPICTION_BUILD)

rdkit_headers(CommercialAbbreviations.h
              SmartAbbreviationPlanner.h
              DEST GraphMol/AdvancedDepiction)

rdkit_catch_test(advancedDepictionCatch catch_tests.cpp
                 LINK_LIBRARIES AdvancedDepiction Depictor Abbreviations SmilesParse GraphMol)

if(RDK_BUILD_BOOST_PYTHON_WRAPPERS)
  add_subdirectory(Wrap)
endif()
