#ifndef AmoreVetoSD_h
#define AmoreVetoSD_h 1

#include "CupSim/CupVetoSD.hh"
#include "AmoreSim/AmoreDetectorConstruction.hh"

class G4TouchableHistory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class AmoreVetoSD: public CupVetoSD 
{
    public:
        AmoreVetoSD(G4String name, int arg_max_tgs = 1000);
        ~AmoreVetoSD();
    public:
        virtual G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);

    protected:
        using eDetGeometry = AmoreDetectorConstruction::eDetGeometry;
};

#endif