#include "globals.hh"

#include "AmoreSim/AmoreDetectorConstruction.hh" // the DetectorConstruction class header
#include "AmoreSim/AmoreVetoSD.hh"
#include "CupSim/CupPMTSD.hh"
#include "CupSim/CupScintSD.hh"            // for making sensitive photocathodes

#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include <sstream>


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
	SDman->SetVerboseLevel(2);
	G4String SDname;

	CupScintSD *TGSD = new CupScintSD(SDname = "/CupDet/TGSD", f200_TotCrystalNum);
	// TGSD = new CupScintSD(SDname = "/cupdet/TGSD", f200_TotCrystalNum);
	SDman->AddNewDetector(TGSD);
	//f200_logiCrystalCell->SetSensitiveDetector(TGSD);
	for(int i = 0; i < f200_TotCrystalNum; i++){
		f200_logiCrystalCell[i]->SetSensitiveDetector(TGSD);
	}

	// CupVetoSD *PSOMD = new CupVetoSD(SDname = "/CupDet/MuonVetoSD", f200_TotPSNum);
	// AmoreVetoSD *PSMD = new AmoreVetoSD(SDname = "/CupDet/MuonVetoSD", f200_TotPSNum*2);
	PSMD = new AmoreVetoSD(SDname = "/CupDet/MuonVetoSD", f200_TotPSNum*2);
	SDman->AddNewDetector(PSMD);
	auto PSO_PV_long = CupDetectorConstruction::GetPhysicalVolumeByName("PlasticScintO_PVlong");
	auto PSI_PV_long = CupDetectorConstruction::GetPhysicalVolumeByName("PlasticScintI_PVlong");
	auto PSO_PV_short = CupDetectorConstruction::GetPhysicalVolumeByName("PlasticScintO_PVshort");
	auto PSI_PV_short = CupDetectorConstruction::GetPhysicalVolumeByName("PlasticScintI_PVshort");

	auto PSO_LV_long = PSO_PV_long->GetLogicalVolume();
	G4cout << "JW: SD check: " << PSO_LV_long->GetName() << G4endl;
	PSO_LV_long->SetSensitiveDetector(PSMD);
	auto PSI_LV_long = PSI_PV_long->GetLogicalVolume();
	PSI_LV_long->SetSensitiveDetector(PSMD);

	auto PSO_LV_short = PSO_PV_short->GetLogicalVolume();
	G4cout << "JW: SD check: " << PSO_LV_short->GetName() << G4endl;
	PSO_LV_short->SetSensitiveDetector(PSMD);
	auto PSI_LV_short = PSI_PV_short->GetLogicalVolume();
	PSI_LV_short->SetSensitiveDetector(PSMD);


	G4Region *PSRegion = new G4Region("PSVeto");
	PSRegion->AddRootLogicalVolume(PSO_LV_long);
	PSRegion->AddRootLogicalVolume(PSI_LV_long);
	PSRegion->AddRootLogicalVolume(PSO_LV_short);
	PSRegion->AddRootLogicalVolume(PSI_LV_short);
	
	

	// CupVetoSD *PSIMD = new CupVetoSD(SDname = "/CupDet/MuonVetoSD", f200_TotPSNum);
	// SDman->AddNewDetector(PSIMD);
	// auto PSI_PV = CupDetectorConstruction::GetPhysicalVolumeByName("PlasticScintI_PV");
	// auto PSI_LV = PSI_PV->GetLogicalVolume();
	// PSI_LV->SetSensitiveDetector(PSMD);

/*
	CupPMTSD *WCMD = new CupPMTSD(SDname = "/cupdet/pmt/inner", f200_TotPMTNum, 0, 10);
	SDman->AddNewDetector(WCMD);
	f200_logiPMTbody->SetSensitiveDetector(WCMD);
	f200_logiPMTinner->SetSensitiveDetector(WCMD);
*/
}
