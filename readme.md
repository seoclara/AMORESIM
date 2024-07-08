# AMORE simulation

## Set your simulation environment
The requirement of this project,
    GEANT4 version = v11.1.3
    ROOT version > v6.22

## Clone and build
    >> git clone https://github.com/seoclara/AMORESIM.git
    >> cd AMORESIM
    >> mkdir build; cd build;
    >> cmake ../.
    >> make -j[NCPU]

## Execution
You have to set LD_LIBRARY_PATH before running.

    >> export LD_LIBRARY_PATH=[build directory/MCObjs]:$LD_LIBRARY_PATH

And, just in the build directory

    >>./amoresim gui.mac