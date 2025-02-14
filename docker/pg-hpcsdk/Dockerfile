FROM nvcr.io/nvidia/nvhpc:25.1-devel-cuda_multi-ubuntu24.04
# Run everything in /root
WORKDIR /root

# Add needed packages
RUN apt-get -y update && apt-get install -y pkg-config g++ gcc curl wget libopenblas-dev \
						libarpack2-dev nano wget git make libssl-dev \
						libhdf5-dev libfftw3-dev libnl-route-3-200 && apt-get clean all;

# Build and Install Latest CMake
RUN curl -O -J -L https://github.com/Kitware/CMake/releases/download/v3.31.1/cmake-3.31.1.tar.gz; \
	apt-get -y purge cmake; \
	tar -xvf ./cmake-3.31.1.tar.gz; \
	cd ./cmake-3.31.1; \
	./bootstrap --prefix=/usr/local/; \
	make -j4; \
	make install; \
	cd ~; \
	rm -rf ./cmake-3.31.1; \
	rm ./cmake-3.31.1.tar.gz;

# Fix errors during install of HPC SDK
#RUN chown -Rv _apt:root /var/cache/apt/archives/partial/; \
#chmod -Rv 700 /var/cache/apt/archives/partial/

#RUN curl https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg

#RUN echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list

#RUN apt-get update -y

#RUN apt-get install -y nvhpc-25-1

# Download and install NVIDIA HPC SDK 25.1
#RUN  wget https://developer.download.nvidia.com/hpc-sdk/21.9/nvhpc-21-9_21.9_amd64.deb \
#	https://developer.download.nvidia.com/hpc-sdk/21.9/nvhpc-2021_21.9_amd64.deb; \
#	apt-get install -y ./nvhpc-2021_21.9_amd64.deb ./nvhpc-21-9_21.9_amd64.deb; \
#	rm ./nvhpc-21-9_21.9_amd64.deb ./nvhpc-2021_21.9_amd64.deb;

# Handle NVIDIA HPC SDK specific stuff
#ENV NVARCH=Linux_x86_64
#ENV NVCOMPILERS=/opt/nvidia/hpc_sdk
#ENV NVSDKVER=25.1
#ENV NVSDKYEAR=2025
#ENV CUDA_PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKVER}/cuda"
#ENV MANPATH="${MANPATH}:${NVCOMPILERS}/${NVARCH}/${NVSDKVER}/compilers/man"
#ENV LD_LIBRARY_PATH="${CUDA_PATH}/lib64:${LD_LIBRARY_PATH}"
#ENV LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKVER}/compilers/lib:${LD_LIBRARY_PATH}"
#ENV LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKYEAR}/compilers/lib:${LD_LIBRARY_PATH}"
#ENV LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKVER}/math_libs/lib64:${LD_LIBRARY_PATH}"
#ENV LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKYEAR}/math_libs/lib64:${LD_LIBRARY_PATH}"
#ENV LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKVER}/comm_libs/mpi/lib:${LD_LIBRARY_PATH}"
#ENV LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKYEAR}/comm_libs/mpi/lib:${LD_LIBRARY_PATH}"
#ENV LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKVER}/comm_libs/nccl/lib:${LD_LIBRARY_PATH}"
#ENV LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKYEAR}/comm_libs/nccl/lib:${LD_LIBRARY_PATH}"
#ENV LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKVER}/comm_libs/nvshmem/lib:${LD_LIBRARY_PATH}"
#ENV LD_LIBRARY_PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKYEAR}/comm_libs/nvshmem/lib:${LD_LIBRARY_PATH}"
#ENV OPAL_PREFIX="${NVCOMPILERS}/${NVARCH}/NVSDKVER/comm_libs/mpi/"
#ENV PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKVER}/compilers/bin:${PATH}"
#ENV PATH="${NVCOMPILERS}/${NVARCH}/${NVSDKVER}/comm_libs/mpi/bin:${PATH}"
#ENV MANPATH="${MANPATH}:${NVCOMPILERS}/${NVARCH}/${NVSDKVER}/comm_libs/mpi/man"

# Install custom compile of boost for mpi.
RUN wget -O ~/boost_1_75_0.tar.bz2 -c 'https://archives.boost.io/release/1.75.0/source/boost_1_75_0.tar.bz2' --max-redirect=100;\
	tar xjf ~/boost_1_75_0.tar.bz2; \
	cd ~/boost_1_75_0/; \
	./bootstrap.sh --prefix=/usr/local --with-toolset=pgi;  \
	cp ~/boost_1_75_0/tools/build/example/user-config.jam ~/user-config.jam; \
	echo "using mpi : mpic++ ; " >> ~/user-config.jam; \
	cd ~/boost_1_75_0/tools/build; \
	./bootstrap.sh && ./b2 install; \
	cd ~/boost_1_75_0 && b2 toolset=pgi --with-mpi --with-program_options --with-serialization stage install; \
	cd ~/; \
	rm -rf ./boost_1_75_0; \
	rm ./boost_1_75_0.tar.bz2;

# Install SuperLU5
RUN curl -O -J -L https://github.com/xiaoyeli/superlu/archive/v5.2.1.tar.gz; \
	tar xvf superlu-5.2.1.tar.gz; \
	mkdir /root/superlu-5.2.1/build; \
	cd /root/superlu-5.2.1/build; \
	cmake ../ -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_C_COMPILER=nvc -DBUILD_SHARED_LIBS=ON; \
	make; \
	make install; \
	cd ~; \
	rm -rf ./superlu-5.2.1/; \
	rm ./superlu-5.2.1.tar.gz;

# Install Armadillo
RUN git clone https://gitlab.com/conradsnicta/armadillo-code; \
	cd armadillo-code; \
	git checkout --track origin/9.900.x; \
	mkdir build; \
	cd build; \
	echo $PATH; \
	cmake ../ -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_C_COMPILER=nvc -DCMAKE_INSTALL_PREFIX=/usr/local; \
	make; \
	make install; \
	ldconfig; \
	cd ~; \
	rm -rf ./armadillo-code;

ENV LD_LIBRARY_PATH="/opt/lib:/opt/lib64:${LD_LIBRARY_PATH}"

RUN git clone https://github.com/acerjanic/ismrmrd.git; \
    cd ismrmrd; \
    git fetch; \
    git checkout matlab; \
    mkdir ./build; \
    cd build; \
    cmake ../ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DBoost_FOUND=OFF; \
    make;  \
    make install; \
	cd ~; \
	rm -rf ./ismrmrd;


ENV LD_LIBRARY_PATH="/opt/PowerGrid/lib:${LD_LIBRARY_PATH}"
ENV PATH="${PATH}:/opt/PowerGrid/bin"

