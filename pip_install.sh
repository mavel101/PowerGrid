#!/bin/bash
# Installs PowerGrid Python bindings depending on configuration
if [ "$#" -lt 1 ]; then
    echo "Provide option: pgi, nonpgi or pgigpu for different build modes"
    exit 1
fi
if [ "$#" -gt 2 ]; then
    echo "Too many arguments passed"
    exit 1
fi

if [ "$2" == "set_hpcsdk" ]; then
    export NVARCH=Linux_x86_64
    export NVCOMPILERS=/opt/nvidia/hpc_sdk
    export NVSDKVER=20.11
    export CUDA_PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKVER}/cuda"
    export MANPATH="${MANPATH}:${NVCOMPILERS}/${NVARCH}/20.11/compilers/man"
    export LD_LIBRARY_PATH="${CUDA_PATH}/lib64:${LD_LIBRARY_PATH}"
    export LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/20.11/compilers/lib:${LD_LIBRARY_PATH}"
    export LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/2020/compilers/lib:${LD_LIBRARY_PATH}"
    export LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/20.11/math_libs/lib64:${LD_LIBRARY_PATH}"
    export LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/2020/math_libs/lib64:${LD_LIBRARY_PATH}"
    export LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/20.11/comm_libs/mpi/lib:${LD_LIBRARY_PATH}"
    export LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/2020/comm_libs/mpi/lib:${LD_LIBRARY_PATH}"
    export LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/20.11/comm_libs/nccl/lib:${LD_LIBRARY_PATH}"
    export LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/2020/comm_libs/nccl/lib:${LD_LIBRARY_PATH}"
    export LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/20.11/comm_libs/nvshmem/lib:${LD_LIBRARY_PATH}"
    export LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/2020/comm_libs/nvshmem/lib:${LD_LIBRARY_PATH}"
    export OPAL_PREFIX="${NVCOMPILERS}/${NVARCH}/20.11/comm_libs/mpi/"
    export PATH="${NVCOMPILERS}/${NVARCH}/20.11/compilers/bin:${PATH}"
    export PATH="${NVCOMPILERS}/${NVARCH}/20.11/comm_libs/mpi/bin:${PATH}"
    export MANPATH="${MANPATH}:${NVCOMPILERS}/${NVARCH}/20.11/comm_libs/mpi/man"
fi

if [ "$1" == "pgi" ]; then
    cp PowerGrid/Python/setup_PGI.py setup.py
    pip install .
    rm setup.py
elif [ "$1" == "nonpgi" ]; then
    cp PowerGrid/Python/setup_nonPGI.py setup.py
    pip install .
    rm setup.py
elif [ "$1" == "pgigpu" ]; then
    cp PowerGrid/Python/setup_PGIGPU.py setup.py
    pip install .
    rm setup.py
else
    echo "Provide option: pgi, nonpgi or pgigpu for different build modes"
    exit 1
fi