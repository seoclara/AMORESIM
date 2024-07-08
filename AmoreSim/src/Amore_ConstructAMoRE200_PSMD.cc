// Written by Jeewon Seo (2024.02.13.)

#include "AmoreSim/AmoreDetectorConstruction.hh" // the DetectorConstruction class header
#include "AmoreSim/AmoreDetectorStaticInfo.hh"

#include "CupSim/CupParam.hh"
#include "CupSim/CupVetoSD.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

// #include "TString.h"

#include <fstream>
#include <sstream>

using namespace CLHEP;
using namespace std;
using namespace AmoreDetectorStaticInfo;
using namespace AmoreDetectorStaticInfo::AMoRE_200;

static G4LogicalVolume *MakePS(G4String type);

void AmoreDetectorConstruction::ConstructAMoRE200_PSMD()
{
    CupParam &db ( CupParam::GetDB() );


	/////////////////////////////////////////////
	// --- Getting the physical volume of the veto housing
	/////////////////////////////////////////////
    G4VPhysicalVolume* VetoHousingSide_PV = GetPhysicalVolumeByName("PlasticVetoHousing_PV"); 
    G4VPhysicalVolume* VetoHousingBottom_PV = GetPhysicalVolumeByName("PlasticVetoHousing2_PV");


	/////////////////////////////////////////////
    // --- Retrieve various geometry parameters
	/////////////////////////////////////////////
    //G4double rockshell_thickness       = db["rockshell_thickness"];
	//G4double airbuffer_radius          = db["airbuffer_radius"];
	G4double airbuffer_height          = db["airbuffer_height"];

	G4double lead_shield_thickness     = db["lead_shield_thickness"];
	G4double lead_housing_thickness    = db["lead_housing_thickness"];

	G4double boricacid_thickness       = db["boricacid_thickness"];

	// G4double plastic_scintillator_thickness = db["plastic_scintillator_thickness"];
	G4double plastic_veto_thickness    = db["plastic_veto_thickness"];
	// G4double plastic_veto_width        = db["plastic_veto_width"];
	// G4double al_plate_thickness        = db["al_plate_thickness"];

	G4double PE_shield_thickness       = db["PE_shield_thickness"];
	G4double thin_lead_shield_thickness= db["thin_lead_shield_thickness"];

	//G4double nShield_hatGapFromLead    = db["nShield_hatGapFromLead"];
	//G4double nShield_GapFromCeiling    = db["nShield_GapFromCeiling"];

	//G4double cavern_sphere_radius      = db["cavern_sphere_radius"];
	//G4double rock_shell_radius         = db["rock_shell_radius"];

	// G4double lead_shield_halfsize    =  // Lead shield  half x,y
		// airbuffer_radius + thin_lead_shield_thickness + lead_shield_thickness; 

	// G4double PS_housing_halfsize =  // Barrel size
		// longPSlength + profile_thickness*2.5 + plastic_veto_thickness;

	G4double PS_housing_height =  airbuffer_height + thin_lead_shield_thickness + lead_shield_thickness 
		+ lead_housing_thickness*2 + boricacid_thickness + PE_shield_thickness;


	/////////////////////////////////////////////
	// --- Visualization attributes
	/////////////////////////////////////////////
	G4VisAttributes *aluminiumVis = new G4VisAttributes(G4Colour(0.6, 0.6, 0.7, 0.1));
	aluminiumVis->SetForceSolid(true);


	/////////////////////////////////////////////
	// --- Read the coordinates of the PSs from the file
	/////////////////////////////////////////////
	ifstream wherePS;
	const char *basic_fn = "pspos_amore200.dat";
	if (getenv("AmoreDATA") != NULL){
		wherePS.open((G4String(getenv("AmoreDATA")) + "/" + G4String(basic_fn)).c_str());
	}
	// print error message on failure of file open
	if (wherePS.fail()) {
		G4cerr << "Error, " << basic_fn << " could not be opened.\n";
		if (getenv("AmoreDATA") == NULL) {
			G4cerr << 
					"AmoreDATA environment variable is not set, so I was looking for"
					<< basic_fn << " in the current directory." 
			<< G4endl;
		} else {
			G4cerr << 
					"I was looking for it in the AmoreDATA directory, "
					<< getenv("AmoreDATA") 
			<< G4endl;
		}
		G4Exception(" ", " ", JustWarning, "Error, ps coordinates file could not be opened.\n");
	}
	// read max number of ps
	int maxPSNo;
	wherePS >> maxPSNo;

	f200_TotPSNum = maxPSNo;

    ////////////////////////////////////////////
	// --- position the plastic scintillators
	////////////////////////////////////////////
	G4int region, ps_id;
	G4double cood_x, cood_y, cood_z, ps_size;

	G4RotationMatrix *psRotMtx = new G4RotationMatrix();
	psRotMtx->rotateZ(-90*deg);

	G4double grad = 3;
	G4RotationMatrix *vpsRotMtx = new G4RotationMatrix();
	vpsRotMtx->rotateY(90*deg);
	vpsRotMtx->rotateX((90+grad)*deg);

	G4RotationMatrix *bpsRotMtx = new G4RotationMatrix();
	bpsRotMtx->rotateY(90*deg);

	while (wherePS.good()) {
		// get a line from the file
		char linebuffer[128];
		wherePS.getline(linebuffer, sizeof(linebuffer) - 1);
		if (wherePS.fail()) break;

		// skip blank lines and lines beginning with '#'
		if (linebuffer[0] == '#' || linebuffer[0] == '\0') continue;

		// put the line in an istrstream for convenient parsing
		istringstream lineStream(linebuffer);

		// parse out region, coordinates,
		region = ps_id = -1;
		cood_x = cood_y = cood_z = ps_size = 0.0;
		lineStream >> region >> ps_id >> cood_x >> cood_y >> cood_z >> ps_size;

		// check for bad data
		if (lineStream.fail() || region < 0 || (cood_x == 0. && cood_y == 0. && cood_z == 0.) || ps_size == 0.) {
			G4cerr << "BAD DATA in PS position file:  line=\"" << linebuffer << "\"\n";
			G4cerr.flush();
			continue;
		}

		// a possible DAQ numbering
		int ith = ps_id;

		// name this PS
		char PSname[64];
		sprintf(PSname, "physMuonVetoPS_Envelope%d", ith);

		// position the PS
		G4ThreeVector pos(cood_x, cood_y, cood_z);

		// place the PS
		if (region < 9){ // side PS
			new G4PVPlacement(G4Transform3D(*psRotMtx, pos), PSname, MakePS("long"), VetoHousingSide_PV, false, ith, OverlapCheck);
		} else if (region < 11){ // bottom PS
			new G4PVPlacement(G4Transform3D(*bpsRotMtx, pos), PSname, MakePS("long"), VetoHousingBottom_PV, false, ith, OverlapCheck);
		} else if (ps_size == longPSlength) { // vertical side PS
			new G4PVPlacement(G4Transform3D(*vpsRotMtx, pos), PSname, MakePS("long"), VetoHousingSide_PV, false, ith, OverlapCheck);
		} else{ // vertical side PS
			new G4PVPlacement(G4Transform3D(*vpsRotMtx, pos), PSname, MakePS("short"), VetoHousingSide_PV, false, ith, OverlapCheck);
		}

		if (ith == 12 || ith == 3*12 || ith == 5*12+17 || ith == 7*12+17) {
			psRotMtx->rotateZ(-90*deg);
		}
		if (ith == 132 || ith == 135 || ith == 138 || ith == 141){
			vpsRotMtx->rotateZ(90*deg);
		}
	}
	wherePS.close();


	/////////////////////////////////////////////
	// --- PS's supporter (aluminium profile)
	/////////////////////////////////////////////
	// veto supporter (aluminium profile)
	G4Box *plasticVetoSupporterV1Box = new G4Box("PlasticVetoSupporterV1_Box",
			profile_thickness/2., profile_thickness/2., PS_housing_height/2.);
	G4Box *plasticVetoSupporterV2Box = new G4Box("PlasticVetoSupporterV2_Box",
			(profile_thickness - plastic_veto_thickness)/2., 
			(profile_thickness - plastic_veto_thickness)/2., 
			(PS_housing_height-(profile_thickness - plastic_veto_thickness))/2.);
	G4Box *plasticVetoSupporterH1Box = new G4Box("PlasticVetoSupporterH1_Box",
			profile_thickness/2., shortPSlength/2., profile_thickness/2.);
	G4Box *plasticVetoSupporterH2Box = new G4Box("PlasticVetoSupporterH2_Box",
			(profile_thickness - plastic_veto_thickness)/2.,
			shortPSlength/4.,
			(profile_thickness - plastic_veto_thickness)/2.);

	G4LogicalVolume *plasticVetoSupporterV1LV = new G4LogicalVolume(plasticVetoSupporterV1Box, 
			_alprofile, "PlasticVetoSupporterV1_LV");
	G4LogicalVolume *plasticVetoSupporterV2LV = new G4LogicalVolume(plasticVetoSupporterV2Box, 
			_alprofile, "PlasticVetoSupporterV2_LV");
	G4LogicalVolume *plasticVetoSupporterH1LV = new G4LogicalVolume(plasticVetoSupporterH1Box, 
			_alprofile, "PlasticVetoSupporterH1_LV");
	G4LogicalVolume *plasticVetoSupporterH2LV = new G4LogicalVolume(plasticVetoSupporterH2Box, 
			_alprofile, "PlasticVetoSupporterH2_LV");
	plasticVetoSupporterV1LV->SetVisAttributes(aluminiumVis);
	plasticVetoSupporterV2LV->SetVisAttributes(aluminiumVis);
	plasticVetoSupporterH1LV->SetVisAttributes(aluminiumVis);
	plasticVetoSupporterH2LV->SetVisAttributes(aluminiumVis);

	/////////////////////////////////////////////
	// --- Read the coordinates of the PS's supporter
	/////////////////////////////////////////////
	ifstream whereProf;
	const char *prof_fn = "profpos_amore200.dat";
	if (getenv("AmoreDATA") != NULL){
		whereProf.open((G4String(getenv("AmoreDATA")) + "/" + G4String(prof_fn)).c_str());
	}
	// print error message on failure of file open
	if (whereProf.fail()) {
		G4cerr << "Error, " << prof_fn << " could not be opened.\n";
		if (getenv("AmoreDATA") == NULL) {
			G4cerr << 
					"AmoreDATA environment variable is not set, so I was looking for"
					<< prof_fn << " in the current directory." 
			<< G4endl;
		} else {
			G4cerr << 
					"I was looking for it in the AmoreDATA directory, "
					<< getenv("AmoreDATA") 
			<< G4endl;
		}
		G4Exception(" ", " ", JustWarning, "Error, profile coordinates file could not be opened.\n");
	}

	G4int prof_id;
	G4double prof_size;
	G4RotationMatrix *profRotMtx = new G4RotationMatrix();
	profRotMtx->rotateZ(-90*deg);

	while (whereProf.good()) {
		// get a line from the file
		char linebuffer[128];
		whereProf.getline(linebuffer, sizeof(linebuffer) - 1);
		if (whereProf.fail()) break;

		// skip blank lines and lines beginning with '#'
		if (linebuffer[0] == '#' || linebuffer[0] == '\0') continue;

		// put the line in an istrstream for convenient parsing
		istringstream lineStream(linebuffer);

		// parse out region, coordinates,
		region = prof_id = -1;
		cood_x = cood_y = cood_z = prof_size = 0.0;
		lineStream >> region >> prof_id >> cood_x >> cood_y >> cood_z >> prof_size;

		// check for bad data
		if (lineStream.fail() || region < 0 || (cood_x == 0. && cood_y == 0. && cood_z == 0.) || prof_size == 0.) {
			G4cerr << "BAD DATA in profile position file:  line=\"" << linebuffer << "\"\n";
			G4cerr.flush();
			continue;
		}

		// a possible DAQ numbering
		int ith = prof_id;

		// name this PS
		char profname[64];
		sprintf(profname, "physPSSupporter%d", ith);

		// position the PS
		G4ThreeVector pos(cood_x, cood_y, cood_z);

		if (prof_size==1){
			new G4PVPlacement(nullptr, pos, profname, plasticVetoSupporterV1LV, VetoHousingSide_PV, false, 0, OverlapCheck);
		} else if (prof_size==2){
			new G4PVPlacement(nullptr, pos, profname, plasticVetoSupporterV2LV, VetoHousingSide_PV, false, 0, OverlapCheck);
		} else if (prof_size==21){
			new G4PVPlacement(G4Transform3D(*profRotMtx, pos), profname, plasticVetoSupporterH1LV, VetoHousingSide_PV, false, 0, OverlapCheck);
		} else if (prof_size==22){
			new G4PVPlacement(G4Transform3D(*profRotMtx, pos), profname, plasticVetoSupporterH2LV, VetoHousingSide_PV, false, 0, OverlapCheck);
		}

		if (ith == 2 || ith == 10 || ith == 16 || ith == 24) {
			profRotMtx->rotateZ(-90*deg);
		}
	}
	whereProf.close();

	/////////////////////////////////////////////
	// --- Print out the geometry information
	/////////////////////////////////////////////
/*
	G4cout << " ======      Muon Veto      ======" << G4endl;
	G4cout << " PS supporter (Aluminium profile)" << G4endl;
	G4cout << "     mass: " << plasticVetoSupporterV1LV->GetMass(true,false)/kg << "(" 
		<< plasticVetoSupporterV1Box->GetXHalfLength()*2. << "x" 
		<< plasticVetoSupporterV1Box->GetYHalfLength()*2. << "x" 
		<< plasticVetoSupporterV1Box->GetZHalfLength()*2. << ") \n" 
		<< "     " << plasticVetoSupporterV2LV->GetMass(true,false)/kg << "(" 
		<< plasticVetoSupporterV2Box->GetXHalfLength()*2. << "x" 
		<< plasticVetoSupporterV2Box->GetYHalfLength()*2. << "x" 
		<< plasticVetoSupporterV2Box->GetZHalfLength()*2. << ") \n" 
		<< "     " << plasticVetoSupporterH1LV->GetMass(true,false)/kg << "(" 
		<< plasticVetoSupporterH1Box->GetXHalfLength()*2. << "x" 
		<< plasticVetoSupporterH1Box->GetYHalfLength()*2. << "x" 
		<< plasticVetoSupporterH1Box->GetZHalfLength()*2. << ") \n" 
		<< "     " << plasticVetoSupporterH2LV->GetMass(true,false)/kg << "(" 
		<< plasticVetoSupporterH2Box->GetXHalfLength()*2. << "x" 
		<< plasticVetoSupporterH2Box->GetYHalfLength()*2. << "x" 
		<< plasticVetoSupporterH2Box->GetZHalfLength()*2. << ")" << G4endl;

	G4cout << " PS veto stainless flame" << G4endl;
	G4cout << "     mass       : " << plasticVetoLV[0]->GetMass(true,false)/kg << "("
		<< plasticVetoBox[0]->GetXHalfLength()*2 << " x "
		<< plasticVetoBox[0]->GetYHalfLength()*2 << " x "
		<< plasticVetoBox[0]->GetZHalfLength()*2 << ") \n"
		<< "     coordinate : " <<  VetoHousingSide_PV->GetTranslation() 
			+ GetPhysicalVolumeByName("physMuonVetoPS0")->GetTranslation() << G4endl;

	G4cout << " PS veto Aluminium plate " << G4endl;
	G4cout << "     mass       : " << aluminiumHolderLV[0]->GetMass(true,false)/kg << "("
		<< plasticScintHolderBox[0]->GetXHalfLength()*2 << "x"
		<< plasticScintHolderBox[0]->GetYHalfLength()*2 << "x"
		<< plasticScintHolderBox[0]->GetZHalfLength()*2 << ")" << G4endl;
	
	G4cout << " Plastic scintillator " << G4endl;
	G4cout << "     mass       : " << plasticScintOLV[0]->GetMass(true,false)/kg << "("
		<< plasticScintBox[0]->GetXHalfLength()*2 << "x"
		<< plasticScintBox[0]->GetYHalfLength()*2 << "x"
		<< plasticScintBox[0]->GetZHalfLength()*2 << ")" << G4endl;
		*/
	
	G4cout << " The total number of PS veto detector = " << maxPSNo << G4endl; 

}

static G4LogicalVolume *MakePS(G4String type)
{
	////////////////////////////////////////////
	// --- Muon Veto (Plastic scintillator)
	////////////////////////////////////////////
	G4bool OverlapCheck = false;
    CupParam &db ( CupParam::GetDB() );
	G4double plastic_scintillator_thickness = db["plastic_scintillator_thickness"];
	G4double plastic_veto_width        = db["plastic_veto_width"];
	G4double al_plate_thickness        = db["al_plate_thickness"];
	G4double plastic_veto_thickness    = db["plastic_veto_thickness"];

	G4double PSlength = 0;
	if (type == "long"){
		PSlength = longPSlength/2.;
	}
	else{
		PSlength = shortPSlength/2.;
	}
	// plastic scintillator
	auto plasticScintBox = new G4Box("PlasticScint_Box", 
				plastic_scintillator_thickness/2., PSlength, plastic_veto_width/2. - veto_frame_thickness);
	auto plasticScintOLV = new G4LogicalVolume(plasticScintBox, G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE"), "PlasticScintO_LV");
	auto plasticScintILV = new G4LogicalVolume(plasticScintBox, G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE"), "PlasticScintI_LV");

	auto plasticScintVis = new G4VisAttributes(G4Colour(0.8, 1.0, 0.725, 0.4));
	plasticScintVis->SetForceSolid(true);
	plasticScintOLV -> SetVisAttributes(plasticScintVis);
	plasticScintILV -> SetVisAttributes(plasticScintVis);

	// aluminium plate
	auto plasticScintHolderBox = new G4Box("AluminiumHolder_Box",
			al_plate_thickness/2., PSlength, plastic_veto_width/2. - veto_frame_thickness);

	auto aluminiumHolderLV = new G4LogicalVolume(plasticScintHolderBox, G4Material::GetMaterial("Aluminium"), "AluminiumHolder_LV");
	aluminiumHolderLV->SetVisAttributes(G4Colour(0.6, 0.6, 0.7, 0.1));

	// air volume
	auto plasticVetoAirBox = new G4Box("PlasticVetoAir_Box",
			plastic_veto_thickness/2. - veto_frame_thickness, PSlength, plastic_veto_width/2. - veto_frame_thickness);
	auto plasticVetoAirLV = new G4LogicalVolume(plasticVetoAirBox, G4Material::GetMaterial("Air"), "PlasticVetoAir_LV");

	// muon veto box
	auto plasticVetoBox = new G4Box("PlasticVeto_Box",
			plastic_veto_thickness/2., PSlength + veto_frame_thickness, plastic_veto_width/2. );
	auto plasticVetoCut = new G4Box("PlasticVetoCut",
			plastic_veto_thickness/2., PSlength + veto_frame_thickness - veto_frame_width,
			plastic_veto_width/2. - veto_frame_width);

	auto plasticVetoSolid1 = new G4SubtractionSolid("PlasticVetoSolid1",
			plasticVetoBox, plasticVetoCut,0,
			{plastic_veto_thickness - veto_frame_thickness/2. + solidBooleanTol,0,0});
	auto plasticVetoSolid = new G4SubtractionSolid("PlasticVetoSolid",
			plasticVetoSolid1, plasticVetoCut,0,
			{-plastic_veto_thickness + veto_frame_thickness/2. - solidBooleanTol,0,0});

	auto PSLV = new G4LogicalVolume(plasticVetoSolid, G4Material::GetMaterial("StainlessSteel"), "PlasticVetoEnvelope_LV");
	auto stainlessVis = new G4VisAttributes(G4Colour(0.6, 0.6, 0.7, 0.9));
	stainlessVis->SetForceSolid(true);
	PSLV->SetVisAttributes(stainlessVis);

	// place the inner components of the Veto Box
	new G4PVPlacement(nullptr, G4ThreeVector(0.,0.,0.), plasticVetoAirLV, 
		"PlasticVetoAir_PV", PSLV, false, 0, OverlapCheck);

	new G4PVPlacement(nullptr, G4ThreeVector(
		plastic_veto_thickness/2. - veto_frame_thickness - al_plate_thickness/2.,0,0), 
		aluminiumHolderLV, "AluminiumHolder_PV", plasticVetoAirLV, false, 0, OverlapCheck);

	new G4PVPlacement(nullptr, G4ThreeVector(
		-plastic_veto_thickness/2. + veto_frame_thickness + al_plate_thickness/2.,0,0), 
		aluminiumHolderLV, "AluminiumHolder_PV", plasticVetoAirLV, false, 0, OverlapCheck);

	new G4PVPlacement(nullptr, G4ThreeVector(
		plastic_veto_thickness/2. - veto_frame_thickness - al_plate_thickness - plastic_scintillator_thickness/2.,0, 0), 
		plasticScintOLV, "PlasticScintO_PV", plasticVetoAirLV, false, 0, OverlapCheck);

	new G4PVPlacement(nullptr, G4ThreeVector(
		-plastic_veto_thickness/2. + veto_frame_thickness + al_plate_thickness + plastic_scintillator_thickness/2.,0, 0),
		plasticScintILV, "PlasticScintI_PV", plasticVetoAirLV, false, 0, OverlapCheck);

	//auto PSRegion = new G4Region("PS_muonveto");
	// PSRegion->AddRootLogicalVolume(plasticScintOLV);
	// PSRegion->AddRootLogicalVolume(plasticScintILV);

	return PSLV;
}