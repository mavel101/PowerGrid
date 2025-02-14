FROM pg-hpcsdk:25.1

# Run everything in /root
WORKDIR /root

# Looks like we lost the PGI environmental variable
ENV PGI=/opt/pgi
RUN echo $PGI; 

RUN git clone https://github.com/mrfil/PowerGrid.git; \
    cd PowerGrid; \
    git fetch; \
    git checkout master; \
    mkdir build; \
    cd build; \
    cmake ../ -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_INSTALL_PREFIX=/opt/PowerGrid; \
    make -j4; make install

# For singularity compatibility
RUN ldconfig





