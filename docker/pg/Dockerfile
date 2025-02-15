FROM ghcr.io/acerjanic/powergrid-hpcsdk:latest

# Run everything in /root
WORKDIR /root

RUN git clone https://github.com/acerjanic/PowerGrid.git; \
    cd PowerGrid; \
    git fetch; \
    git checkout HPCSDK_25.1; \
    mkdir build; \
    cd build; \
    cmake ../ -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_INSTALL_PREFIX=/opt/PowerGrid; \
    make -j4; make install

# Add paths to LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH="/opt/PowerGrid/lib:${LD_LIBRARY_PATH}"
ENV PATH="${PATH}:/opt/PowerGrid/bin"

# For singularity compatibility
RUN ldconfig





