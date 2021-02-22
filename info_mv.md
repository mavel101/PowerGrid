# PowerGrid Compiling & Debugging Infos

## Compiling

For Compiling install packages on https://github.com/mrfil/PowerGrid.  
Armadillo should be version 9.900.x 
Compiling is working with newer versions of ISMRMRD (e.g. 1.4.2).  
Cmake should be compiled with OPENACC_GPU=OFF und OPENACC_MP=ON/OFF if no CUDA is available (change in CMakeLists.txt).  

For Compiling:  
mkdir build  
cd build  
cmake ..  
make  
sudo make install  

For Compiling with PGI compilers run:
cmake .. -DCMAKE_CXX_COMPILER=pgc++ instead of cmake ..

The PGI Compiler has to be in the PATH. See docker/pg-hpcsdk/Dockerfile for information on how to set the PATH variables.

## Compile Python wrapper

Python wrapper depends on Pybind11. Pybind11 is added as a submodule and can be installed with
cd pybind11
mkdir build
cd build
cmake ..
make
sudo make install

For installing the wrapper just run `pip install .` from the PowerGrid root directory. This installs the module PowerGridPy..  
For information run:

```python
from PowerGridPy import PowerGridIsmrmrd
PowerGridIsmrmrd.__doc__
```

If the Wrapper should be compiled with PGI compilers, use setupPGI.py instead of setup.py. The PGI Compiler has to be in the PATH.
IMPORTANT: For compiling with Pybind11 & PGI, the CXX Standard in CMakeLists.txt has to be CMAKE_CXX_STANDARD 14.

## Debugging in VS Code

in tasks.json:

```cpp
{
    "version": "2.0.0",
    "tasks": [
        {
            "label": "build",
            "type": "shell",
            "command": "cd build && cmake .. && make",
            "problemMatcher": [
                "$gcc"
            ]
        }
    ]
}
```

in launch.json (Dateien jeweils ersetzen):

```cpp
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "(gdb) Starten",
            "type": "cppdbg",
            "request": "launch",
            "program": "/home/veldmannm/Software/PowerGrid/build/PowerGridIsmrmrd",
            "args": [
                "-i",
                "/home/veldmannm/Projects/SpiralGirfProject/RawData/20201210/in_vivo/power_grid_test_MID251.h5",
                "-o",
                "/home/veldmannm/Projects/SpiralGirfProject/RawData/20201210/in_vivo/powergridimages",
                "-F",
                "NUFFT",
                "-t",
                "1",
                "-B",
                "0.1"],
            "stopAtEntry": false,
            "cwd": "${workspaceFolder}",
            "environment": [],
            "externalConsole": false,
            "MIMode": "gdb",
            "setupCommands": [
                {
                    "description": "Automatische Strukturierung und Einrückung für \"gdb\" aktivieren",
                    "text": "-enable-pretty-printing",
                    "ignoreFailures": true
                }
            ]
        }
    ]
}
```
