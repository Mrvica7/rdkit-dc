rdkit_library(Abbreviations
              Abbreviations.cpp
              AbbreviationsUtils.cpp
              LINK_LIBRARIES SmilesParse
              SubstructMatch GraphMol RDGeneral)
target_compile_definitions(Abbreviations PRIVATE RDKIT_ABBREVIATIONS_BUILD)

rdkit_headers(Abbreviations.h DEST GraphMol/Abbreviations)

rdkit_catch_test(testAbbreviations catch_tests.cpp
LINK_LIBRARIES Abbreviations SmilesParse FileParsers )

if(RDK_BUILD_BOOST_PYTHON_WRAPPERS)
  add_subdirectory(Wrap)
endif()

if(RDK_BUILD_NANOBIND_WRAPPERS)
  add_subdirectory(nbWrap)
endif()

# AdvancedDepiction depends on both the base Depictor and Abbreviations. Build
# it here so the dependency direction remains one-way and neither lower-level
# library needs to know about the high-level planner.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../AdvancedDepiction
                 ${CMAKE_CURRENT_BINARY_DIR}/../AdvancedDepiction)
