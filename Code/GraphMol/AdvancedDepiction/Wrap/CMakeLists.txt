rdkit_python_extension(rdAdvancedDepiction
                       rdAdvancedDepiction.cpp
                       DEST Chem
                       LINK_LIBRARIES AdvancedDepiction)

add_pytest(pyAdvancedDepiction ${CMAKE_CURRENT_SOURCE_DIR}/testAdvancedDepiction.py)
