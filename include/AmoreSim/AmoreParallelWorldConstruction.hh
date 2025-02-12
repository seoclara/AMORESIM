#ifndef AmoreParallelWorldConstruction_HH
#define AmoreParallelWorldConstruction_HH 1

#include "AmoreSim/AmoreDetectorConstruction.hh"
#include "G4VUserParallelWorld.hh"
#include "globals.hh"

#include "CupSim/CupDetectorConstruction.hh"

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class CupPMTSD;

class AmoreParallelWorldConstruction : public G4VUserParallelWorld {
    public:
        AmoreParallelWorldConstruction(G4String& parallelWorldName);
        virtual ~AmoreParallelWorldConstruction();

        virtual void Construct();
        virtual void ConstructSD();

    private:
        //G4LogicalVolume* fWCPMT_Logi;
        //G4VPhysicalVolume* fWCPMT_phys;
        //G4LogicalVolume* ghostLogical;
        G4Material *fWater;
        G4Material *fBlackAcryl ;
        G4Material *fGlass;
        G4Material *fStainless;
        G4Material *fPMT_Vac;
        G4Material *dummyMat;

    protected:
        using eDetGeometry  = AmoreDetectorConstruction::eDetGeometry;
        using eSimulationType = AmoreDetectorConstruction::eSimulationType;
        void ConstructAMoRE200_ParallelWorldForWCMD();

        CupPMTSD* pmtSDWC;
        G4int maxWCPMTNo;

};

#endif
