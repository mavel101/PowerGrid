include_directories(../../Support)
include_directories(../)

pybind11_add_module(PowerGridPy MODULE ./PowerGrid/Python/PowerGridPy.cpp)
target_link_libraries(PowerGridPy PRIVATE ${LIBS} PowerGrid)

if (MPISupport)
    include_directories(../MPI)
    pybind11_add_module(PowerGridPyMPI MODULE ../Python/PowerGridPyMPI.cpp)
    target_link_libraries(PowerGridPyMPI PRIVATE ${mpiLIBS} PowerGrid)
endif()
