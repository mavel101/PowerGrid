# PowerGrid Compiling & Debugging Infos

## Compiling

For Compiling install packages listed in docker/powergrid-dev.  
Armadillo should be version 9.900.x .
Compiling is working with newer versions of ISMRMRD (e.g. 1.4.2).  

For Compiling:  
mkdir build  
cd build  
cmake ..  
make  
sudo make install  

## Compile Python wrapper

Python wrapper depends on Pybind11. Pybind11 is added as a submodule, which has to be pulled with git submodule update --init. 
There is no need to install Pybind11.

For installing the wrapper use the setup.py's from PowerGrid/Python. There are 3 different options for compiling:

- "setup_CPU.py" : compile for multi-core CPU support
- "setup_GPU.py": compile for GPU support (only on Nvidia GPU)

For compiling with PGI, Nvidia HPC-SDK 20.11 has to be installed and the compiler has to be in the PATH.
Pass "set_hpcsdk" as 2nd argument to pip_install.sh to set environment variables for the PGI compiler (might have to check the compiler paths in the script)
Boost, SuperLU5 and Armadillo can optionally also be compiled with PGI compilers, see Dockerfile for instructions.
For compiling with Pybind11 & PGI, the CXX Standard in CMakeLists.txt has to be CMAKE_CXX_STANDARD 14.

For an easier build, install Docker & Nvidia Docker and use Dockerfiles from https://github.com/pehses/python-ismrmrd-server/tree/pulseq/docker.

Import the module inside Python with:

```python
from PowerGridPy import PowerGridIsmrmrd
PowerGridIsmrmrd.__doc__
```


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

in launch.json:

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
