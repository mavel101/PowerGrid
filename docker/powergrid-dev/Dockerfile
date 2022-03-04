#FROM nvcr.io/hpc/pgi-compilers:ce as builder
FROM nvcr.io/nvidia/nvhpc:22.1-devel-cuda_multi-ubuntu20.04

# Run everything in /root
WORKDIR /root

# Add needed packages
RUN apt-get -y update && apt-get install -y g++ gcc curl wget libopenblas-dev \
						libarpack2-dev nano wget git \
						libhdf5-dev libfftw3-dev && apt-get clean all;

#ENV LD_LIBRARY_PATH="/opt/pgi/linux86-64/2019/cuda/9.1/lib64:${LD_LIBRARY_PATH}"
ENV BASH_ENV=/usr/share/lmod/lmod/init/bash
RUN ["/bin/bash", "-c", "module load nvhpc; module save default;"]

# Install Latest CMake
RUN curl -O -J -L https://github.com/Kitware/CMake/releases/download/v3.22.2/cmake-3.22.2-linux-x86_64.tar.gz; \
	apt-get -y purge cmake; \
	tar -xvf ./cmake-3.22.2-linux-x86_64.tar.gz; \
	cd ./cmake-3.22.2-linux-x86_64; \
	cp -r bin /usr/; \
	cp -r doc /usr/share/; \
	cp -r man /usr/share/; \
	cp -r share /usr/

	# Install custom compinvcr.io/nvidia/nvhpc:22.1-devel-cuda_multi-ubuntu20.04le of boost for mpi.
RUN wget -O ~/boost_1_70_0.tar.bz2 -c 'https://boostorg.jfrog.io/artifactory/main/release/1.70.0/source/boost_1_70_0.tar.bz2' --max-redirect=100\
 && tar xjf ~/boost_1_70_0.tar.bz2 \
 && cd ~/boost_1_70_0/ \
 && ./bootstrap.sh --prefix=/usr/local --with-toolset=gcc --with-libraries=mpi,program_options,serialization \
 && cp ~/boost_1_70_0/tools/build/example/user-config.jam ~/boost_1_70_0/tools/build/src/user-config.jam \
 && echo "using mpi : mpic++ ; " >> ~/boost_1_70_0/tools/build/src/user-config.jam

RUN cd ~/boost_1_70_0 && ./b2 -j 4 toolset=gcc install || true

# Install SuperLU5
RUN curl -O -J -L https://github.com/xiaoyeli/superlu/archive/v5.2.1.tar.gz; \
	tar xvf superlu-5.2.1.tar.gz; \
	mkdir /root/superlu-5.2.1/build; \
	cd /root/superlu-5.2.1/build; \
	cmake ../ -DBUILD_SHARED_LIBS=ON; \
	make; \
	make install;

# Install Armadillo
RUN git clone https://gitlab.com/conradsnicta/armadillo-code; \
	cd armadillo-code; \
	git checkout 9.200.x; \
	mkdir build; \
	cd build; \
	echo $PATH; \
	cmake ../ -DCMAKE_INSTALL_PREFIX=/opt; \
	make; \
	make install; \
	ldconfig

ENV LD_LIBRARY_PATH="/opt/lib:/opt/lib64:${LD_LIBRARY_PATH}"

RUN git clone https://github.com/acerjanic/ismrmrd.git; \
    cd ismrmrd; \
    git fetch; \
    git checkout matlab; \
    mkdir ./build; \
    cd build; \
    cmake ../ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DBoost_FOUND=OFF; \
    make;  \
    make install

ENV LD_LIBRARY_PATH="/opt/PowerGrid/lib:${LD_LIBRARY_PATH}"
ENV PATH="${PATH}:/opt/PowerGrid/bin"

# For singularity compatibility
RUN ldconfig
