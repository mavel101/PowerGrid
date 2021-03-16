include_directories(../../Support)
include_directories(../)

pybind11_add_module(PowerGridPy MODULE PowerGridPy.cpp)
target_link_libraries(PowerGridPy PRIVATE ${LIBS} PowerGrid)

# The MPI version of PowerGrid did not work on my system
# If the MPI version works, the Bindings should work as well, but I didnt test them
# so MPI Support is by default turned of for Python Bindings (in setup.py 's)
if (MPISupport)
    include_directories(../MPI)
    pybind11_add_module(PowerGridPyMPI MODULE PowerGridPyMPI.cpp)
    target_link_libraries(PowerGridPyMPI PRIVATE ${mpiLIBS} PowerGrid)
endif()
