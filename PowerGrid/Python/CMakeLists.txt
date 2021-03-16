include_directories(../../Support)
include_directories(../)

pybind11_add_module(PowerGridPy MODULE PowerGridPy.cpp)
target_link_libraries(PowerGridPy PRIVATE ${LIBS} PowerGrid)
