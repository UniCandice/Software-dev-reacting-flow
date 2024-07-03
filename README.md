
# Code_Saturne Software

## Introduction

**Code_Saturne** is an open-source computational fluid dynamics (CFD) software developed by EDF (Électricité de France). It is designed to handle a wide range of fluid flow simulations, including laminar and turbulent flows, heat transfer, and multiphase flows. The software is robust, highly flexible, and can be used for industrial applications, research, and educational purposes.

## Features

- **Multi-physics Capabilities**: Supports simulations involving fluid dynamics, thermal transfers, and coupling with other physics such as electromagnetism and structure mechanics.
- **Turbulence Modeling**: Includes various turbulence models such as RANS (Reynolds-Averaged Navier-Stokes), LES (Large Eddy Simulation), and DNS (Direct Numerical Simulation).
- **Multiphase Flows**: Ability to simulate flows involving multiple fluid phases (liquid, gas, particles).
- **Mesh Flexibility**: Can handle structured and unstructured meshes, with support for various mesh formats.
- **Parallel Computing**: Efficient parallel computing capabilities to handle large-scale simulations on high-performance computing (HPC) systems.
- **User Customization**: Allows for user-defined functions and advanced customization of simulations.

## Installation

### Prerequisites

- **Operating System**: Linux, macOS, Windows (via WSL)
- **Compilers**: GCC, Intel, or Clang
- **Dependencies**: CMake, MPI library (such as OpenMPI or MPICH), HDF5, MEDCoupling, and optionally ParMETIS, Scotch, and others for specific features.

### Installation Steps

1. **Clone the Repository**:
    ```sh
    git clone https://github.com/code-saturne/code_saturne.git
    cd code_saturne
    ```

2. **Create a Build Directory**:
    ```sh
    mkdir build
    cd build
    ```

3. **Configure the Build**:
    ```sh
    cmake ..
    ```

4. **Build the Software**:
    ```sh
    make
    ```

5. **Install the Software**:
    ```sh
    sudo make install
    ```

6. **Set Up Environment Variables**:
    Add the installation path to your `PATH` and `LD_LIBRARY_PATH` in your shell configuration file (e.g., `.bashrc` or `.zshrc`):
    ```sh
    export PATH=/path/to/code_saturne/bin:$PATH
    export LD_LIBRARY_PATH=/path/to/code_saturne/lib:$LD_LIBRARY_PATH
    ```

## Getting Started

1. **Create a New Case**:
    ```sh
    code_saturne create --study my_study --case my_case
    cd my_study/my_case/DATA
    ```

2. **Configure the Simulation**:
    - Edit the configuration files (`.xml` and `.cfg` files cs_user_module.f90 and cs_user_parameters.f90) to define the simulation parameters, such as mesh, physical properties, boundary conditions, and numerical settings.

3. **compile the code**:
    ```sh
    code_saturne compile
    ```

4. **Run the Simulation**:
    ```sh
    code_saturne run --initialization
    ```

5. **Post-process the Results**:
    - Use visualization tools like ParaView or post-processing scripts to analyze the simulation results.

## Documentation and Support

- **User Guide**: Comprehensive user guide and tutorials available at the [Code_Saturne Documentation](http://code-saturne.org/documentation)
- **Community Support**: Join the [Code_Saturne Forum](http://code-saturne.org/forum) for community support and discussions.
- **Contributing**: Contributions are welcome! Please read the contributing guidelines and follow the code of conduct.


## Contributions

Liyuan LIU has developed the SRC for reacting flow, focusing on the combustion models module development. The methodology can refer to:

- **Assessment of Bray Moss Libby formulation for premixed flame-wall interaction within turbulent boundary layers: Influence of flow configuration**
- **Statistical Behaviour and Modelling of Variances of Reaction Progress Variable and Temperature During Flame‑Wall Interaction of Premixed Flames Within Turbulent Boundary Layers**
