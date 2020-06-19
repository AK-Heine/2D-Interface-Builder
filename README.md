# 2D-interface-builder

Builds interface of 2D structures via coincidene lattice theory.

## Installation

Requires a python3 installation with python header files, for example with the [Anaconda package manager](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html).

Requires a compiler for C++11 and [cmake](https://cmake.org/) at least with version 3.17.2. On Unix and WSL, you can install a compiler with:

```bash
sudo apt-get update
sudo apt-get install build-essential gdb
whereis g++
whereis gdb
```

Note that installing cmake with most package managers gives an out-dated version. I recommend downloading and installing the newest version of [cmake](https://cmake.org/) directly, even though the installation might take a while.

Then, download or clone the .zip from here and run:
```bash
pip install 2D-interface-builder.zip
```

## Documentation

Documentation is available at [ReadTheDocs](https://2d-interface-builder.readthedocs.io).

## Testing

Tests can be run in the project directory with

```bash
pytest -v tests
```

More tests will follow.

## Requirements

- [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/) ase 3.19 or higher
- [Space Group Libary](https://atztogo.github.io/spglib/python-spglib.html) spglib
- [Scientific Python](https://www.scipy.org/) scipy (including numpy, matplotlib and pandas)
- [networkx](https://networkx.github.io/documentation/stable/install.html)
