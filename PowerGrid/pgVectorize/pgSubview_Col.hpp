// pgSubview_Col.hpp

#include "PGIncludes.h"

template<typename T>
class subview_col<T> {

private:

T *mem; //Pointer to raw data
bool isInitialized; 
bool isOnGPU; 

// Number of elements in array
unsigned int n_elem;

public:

unsigned int n_elem() const{
    return n_elem;
}

void zeros() {
    #pragma acc parallel loop present(mem[0:n_elem])
    for(unsigned int ii = 0; ii < n_elem; ii++) {
        mem[ii] = T();
    }
}

void ones() {
    #pragma acc parallel loop present(mem[0:n_elem])
    for(unsigned int ii = 0; ii < n_elem; ii++) {
        mem[ii] = T(1.0);
    }
}

operator=(const pgCol<T>& d) {
    if()
}

operator()(const unsigned int d) {
    return
}

const operator()(const unsigned int d) const {
    return
}

};


