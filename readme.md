# AMORE simulation

## Set your simulation environment
The requirement of this project,

    GEANT4 version = v11.1.3
    ROOT version > v6.22

## Clone and build
    >> git clone https://github.com/seoclara/AMORESIM.git
    >> cd AMORESIM
    >> mkdir build; cd build;
    >> ccmake ../.
    >> cd AmoresSim
    >> source environment.sh
    >> make -j[NCPU]

## Execution
There are some macros in the AmoreSim directory. In order to check the geometry, you can turn on the gui window and draw the detectors by 

    >>./amoresim gui.mac

If you want run the simulation, you need to move run directory and execute one of the macro. 

run_II_[simulation mode].sh is for AMoRE-II, run_I_[simulation mode].sh is for AMoRE-I, run_Pilot_[simulation mode].sh is for AMoRE-Pilot simulation .

    >>./run_I_decay.sh 0 10