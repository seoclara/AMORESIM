#include "globals.hh"

#include "AmoreSim/AmoreDetectorConstruction.hh" // the DetectorConstruction class header
#include "AmoreSim/AmoreVetoSD.hh"
#include "CupSim/CupPMTSD.hh"
// #include "CupSim/CupParam.hh"
#include "CupSim/CupScintSD.hh"            // for making sensitive photocathodes
// #include "CupSim/CupVetoSD.hh"             // for making sensitive photocathodes
// #include "CupSim/CupTorusStack.hh"         // for making the cavern and LSC envelope
// #include "CupSim/Cup_PMT_LogicalVolume.hh" // for making PMT assemblies

/*
#include "G4Box.hh"
#include "G4Element.hh"
#include "G4IntersectionSolid.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4Sphere.hh" // for making spheres
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
*/

//#include "CLHEP/Matrix/Matrix.h"
// #include "G4Colour.hh"
// #include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
// #include "G4PVPlacement.hh"
// #include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
// #include "G4SmartVoxelHeader.hh"
#include "G4ThreeVector.hh"
// #include "G4VPhysicalVolume.hh"
// #include "G4VisAttributes.hh"

// #include "G4ios.hh"

// #include "G4Types.hh"

#include <sstream>

//#include "Randomize.hh"                     // for G4UniformRand()

//#include "fstream"
//#include "strstream"
//#include <stdio.h>

// == Construct Geometry for AMoRE-200 ======================================
// uses parameters from database or file
void AmoreDetectorConstruction::ConstructAMoRE200() {
	// --- put in the Inner Detector tanks, PMTs, and details
	G4LogicalVolume* odLV = ConstructAMoRE200_OD();

	// --- put in the Inner Detector tanks, PMTs, and details
	ConstructAMoRE200_ID(odLV);

	//if(whichSimType!=kRockGammaMode){
		// --- put the Plastics Scintillator for Veto
	ConstructAMoRE200_PSMD();
	//}

}

void AmoreDetectorConstruction::ConstructAMoRE200_SDandField() {
	//////////////////////////////
	// --- TG sensitive detector
	//////////////////////////////
	// get pointer to logical volume

	G4SDManager *SDman = G4SDManager::GetSDMpointer();
	G4String SDname;

	CupScintSD *TGSD = new CupScintSD(SDname = "/CupDet/TGSD", f200_TotCrystalNum);
	// TGSD = new CupScintSD(SDname = "/cupdet/TGSD", f200_TotCrystalNum);
	SDman->AddNewDetector(TGSD);
	//f200_logiCrystalCell->SetSensitiveDetector(TGSD);
	for(int i = 0; i < f200_TotCrystalNum; i++){
		f200_logiCrystalCell[i]->SetSensitiveDetector(TGSD);
	}

	// CupVetoSD *PSOMD = new CupVetoSD(SDname = "/CupDet/MuonVetoSD", f200_TotPSNum);
	AmoreVetoSD *PSMD = new AmoreVetoSD(SDname = "/CupDet/MuonVetoSD", f200_TotPSNum*2);
	SDman->AddNewDetector(PSMD);
	auto PSO_PV = CupDetectorConstruction::GetPhysicalVolumeByName("PlasticScintO_PV");
	auto PSO_LV = PSO_PV->GetLogicalVolume();
	PSO_LV->SetSensitiveDetector(PSMD);

	// CupVetoSD *PSIMD = new CupVetoSD(SDname = "/CupDet/MuonVetoSD", f200_TotPSNum);
	// SDman->AddNewDetector(PSIMD);
	auto PSI_PV = CupDetectorConstruction::GetPhysicalVolumeByName("PlasticScintI_PV");
	auto PSI_LV = PSI_PV->GetLogicalVolume();
	PSI_LV->SetSensitiveDetector(PSMD);

/*
	CupPMTSD *WCMD = new CupPMTSD(SDname = "/cupdet/pmt/inner", f200_TotPMTNum, 0, 10);
	SDman->AddNewDetector(WCMD);
	f200_logiPMTbody->SetSensitiveDetector(WCMD);
	f200_logiPMTinner->SetSensitiveDetector(WCMD);
*/
}
