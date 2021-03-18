include_directories(../../Support)
include_directories(../)

pybind11_add_module(PowerGridPy MODULE PowerGridPy.cpp)
target_link_libraries(PowerGridPy PRIVATE ${LIBS} PowerGrid)

# I didnt manage to get the MPI version of PowerGrid running on my system
# So Python bindings are not in the build included yet - MV