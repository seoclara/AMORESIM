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
    >> cd AmoreSim
    >> source environment.sh
    >> make -j[NCPU]

## Execution
There are some macros in the AmoreSim directory. In order to check the geometry, you can turn on the gui window and draw the detectors by 

    >>./amoresim gui.mac

If you want run the simulation, you need to move run directory and excute one of the macro. There are three AMoRE phase (AMoRE-II, AMoRE-I and AMoRE-Pilot) and three simulation mode (decay, muon and neutron).

run_II_[simulation_mode].sh is for AMoRE-II,   
run_I_[simulation_mode].sh is for AMoRE-I,   
run_Pilot_[simulation_mode].sh is for AMoRE-Pilot simulation.   

For example, you can do radio active background event simulation for AMoRE-I as:

    >>./run_I_decay.sh 0 10
