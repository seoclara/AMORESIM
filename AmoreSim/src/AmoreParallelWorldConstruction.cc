#include "AmoreSim/AmoreParallelWorldConstruction.hh"
#include "AmoreSim/AmoreDetectorConstruction.hh"

#include "CupSim/CupParam.hh"
#include "CupSim/CupPMTSD.hh" // for "sensitive detector"
#include "CupSim/CupPMTOpticalModel.hh"
#include "CupSim/CupDetectorConstruction.hh"

#include "G4SDManager.hh"
#include "G4RegionStore.hh"
#include "G4FastSimulationManagerProcess.hh"

#include "G4Box.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4OpticalSurface.hh"
#include "G4Region.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"


AmoreParallelWorldConstruction::AmoreParallelWorldConstruction(G4String& parallelWorldName)
: G4VUserParallelWorld(parallelWorldName) {;}

AmoreParallelWorldConstruction::~AmoreParallelWorldConstruction() {;}


////////////////////////////////////////////////////////////////////
// CONSTRUCT PARALLAL DETECTOR WORLD
//
void AmoreParallelWorldConstruction::Construct(){
    G4cout << "========================================" << G4endl;
    G4cout << "Construct Parallel World....." << G4endl;
    G4cout << "========================================" << G4endl;
    G4cout << G4endl;

    CupParam &db(CupParam::GetDB());
    G4double nShield_hatGapFromLead = db["nShield_hatGapFromLead"];
    eSimulationType SimType = AmoreDetectorConstruction::GetSimType();
 
    // Call the realworld to get the world volume
    G4VPhysicalVolume* ghostWorld = GetWorld();
    G4LogicalVolume* ghostLogical = ghostWorld->GetLogicalVolume();
    //ghostLogical->SetMaterial(0); // The material of real world will be used for ghost volume.

    // Materials
    dummyMat  = 0;
    fWater = G4Material::GetMaterial("Water");
    fBlackAcryl = G4Material::GetMaterial("Acrylic");
    fGlass = G4Material::GetMaterial("Glass");
    fStainless = G4Material::GetMaterial("StainlessSteel");
    fPMT_Vac = G4Material::GetMaterial("PMT_Vac");

    // Create Para volume
    G4Box* paraBox = new G4Box("paraBox", 6600./2.*mm, 6600./2.*mm, 3300./2.*mm);
    G4LogicalVolume* paraLogical = new G4LogicalVolume(paraBox, dummyMat, "paraLogical");
    if (SimType == eSimulationType::kRockGammaMode){
        new G4PVPlacement(0, G4ThreeVector(0,0,paraBox->GetZHalfLength()+nShield_hatGapFromLead), paraLogical, "paraPhysical", ghostLogical, false, 0);
    } else if (SimType == eSimulationType::kNeutronMode){
        new G4PVPlacement(0, G4ThreeVector(0,0,paraBox->GetZHalfLength()+nShield_hatGapFromLead), paraLogical, "paraPhysical", ghostLogical, false, 0);
    } else {
        new G4PVPlacement(0, G4ThreeVector(0,0,5880), paraLogical, "paraPhysical", ghostLogical, false, 0);
    }

    
 
    ////////////////////////////////////////////////////
    // --- put in the Inner Detector PMTs, and details
    switch (AmoreDetectorConstruction::GetDetGeometryType()){
        case eDetGeometry::kDetector_AMoRE200:
            G4cout << "Constructing AMoRE200 parallel world" << G4endl;
            ConstructAMoRE200_ParallelWorldForWCMD();
            break;
        default:
            break;
    }
}

////////////////////////////////////////////////////////////////////
// CONSTRUCT PARALLAL WORLD SD
//
void AmoreParallelWorldConstruction::ConstructSD(){
        eDetGeometry whichDetGeometry = AmoreDetectorConstruction::GetDetGeometryType();
    	switch (whichDetGeometry){
          case eDetGeometry::kDetector_AMoRE200: {
            G4cout << "Constructing AMoRE200 parallel world SD" << G4endl; 

            // Set PMT sensitive detector
            G4SDManager* SDman = G4SDManager::GetSDMpointer();
            G4String pmtSDWCname = "/cupdet/pmt/inner";
            pmtSDWC = new CupPMTSD(pmtSDWCname, maxWCPMTNo, 0, 10); 
            SDman->AddNewDetector(pmtSDWC);

	        // Set PMT Optical modle for AMoRE-II WCMD
            G4cout << "Set PMT Optical Model for AMoRE-II WCMD" << G4endl;
            G4VPhysicalVolume* PMT_body_phys = CupDetectorConstruction::GetPhysicalVolumeByName("WCPMT_body_phys");
            if (PMT_body_phys == NULL) {
                G4Exception(" ", " ", JustWarning, "Could not find PMT_body_phys!");
            } else G4cout << "found PMT_body_phys: " << PMT_body_phys->GetName() << G4endl;

            G4LogicalVolume* PMT_body_log = PMT_body_phys->GetLogicalVolume();
            PMT_body_log->SetSensitiveDetector(pmtSDWC);

            G4VPhysicalVolume* PMT_inner1_phys = CupDetectorConstruction::GetPhysicalVolumeByName("WCPMT_inner1_phys");
            G4LogicalVolume* PMT_inner1_log = PMT_inner1_phys->GetLogicalVolume();
            PMT_inner1_log->SetSensitiveDetector(pmtSDWC);

            new CupPMTOpticalModel("WCPMT_optical_model", PMT_body_phys);
            G4cout << "PMT Optical Model done" << G4endl;
            break;
        }
        default: { break; }
    }
    
}