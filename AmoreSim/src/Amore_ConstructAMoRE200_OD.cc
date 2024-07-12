//
//  Original by G. Horton-Smith 2004/12/02
//
//  Modified by E.J.Jeon 2007/06/14
//  Modified by Y.S.Yoon 2015/06/15
//  Updated by J.Seo 2024/04/01

#include "globals.hh"

#include "AmoreSim/AmoreDetectorConstruction.hh" // the DetectorConstruction class header
#include "AmoreSim/AmoreDetectorStaticInfo.hh"
#include "AmoreSim/AmoreEventAction.hh"
#include "CupSim/CupPMTSD.hh" // for "sensitive detector"
#include "CupSim/CupParam.hh"
#include "CupSim/CupScintSD.hh"
#include "CupSim/Cup_PMT_LogicalVolume.hh" // for making PMT assemblies

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4RotationMatrix.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4ExtrudedSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4Colour.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"

#include <fstream>
#include <sstream>

using namespace std;

////////////////////////////////////////////////////////////////
// declaration of "private" static utility functions that we
// don't need in class definition
static G4LogicalVolume *MakeHBeam1(G4Material *mat);
static G4LogicalVolume *MakePit(G4double pitBox_x, G4double pitBox_z, G4Material *mat, G4double HatBarrel_gap);

////////////////////////////////////////////////////////////////
/** ConstructXeDetector_OD() constructs the GenericLAND inner detector.
	It is a member of the AmoreDetectorConstruction class, so it has
	direct access to all the member variables in that class, including
	saved material pointers and such.
	*/
G4LogicalVolume *AmoreDetectorConstruction::ConstructAMoRE200_OD()
{
	using namespace CLHEP;
	using namespace AmoreDetectorStaticInfo;
	using namespace AmoreDetectorStaticInfo::ColorTable;
	using namespace AmoreDetectorStaticInfo::AMoRE_200;

	// -- database
	CupParam &db(CupParam::GetDB());

	// -- Visualization Attributes
	G4VisAttributes *visRock = new G4VisAttributes(G4Colour(0.0, 0.8, 0.8, 0.3));
	G4VisAttributes *visCavern = new G4VisAttributes(orangel);
	G4VisAttributes *visAir = new G4VisAttributes(G4Colour(1, 1, 1, 0.4));
	G4VisAttributes *leadVis = new G4VisAttributes(G4Colour(0.0, 0.0, 0.7, 0.3));
	G4VisAttributes *stainlessVis = new G4VisAttributes(G4Colour(0.6, 0.6, 0.7, 0.9));
	G4VisAttributes *aluminiumVis = new G4VisAttributes(G4Colour(0.6, 0.6, 0.7, 0.1));
	G4VisAttributes *ironVis = new G4VisAttributes(G4Colour(0, 0, 0, 0.1));
	G4VisAttributes *CuShieldVisAttr = new G4VisAttributes(G4Colour(23 / 255., 189 / 255., 79 / 255., 0.8));
	G4VisAttributes *boricAcidVisAttr = new G4VisAttributes(G4Colour(0.4, 0.4, 0.4, 0.2));
	G4VisAttributes *boricAcidCaseVisAttr = new G4VisAttributes(G4Colour(0.7, 0.2, 0.4, 0.2)); // su-yeon
	G4VisAttributes *plasticScintVis = new G4VisAttributes(G4Colour(0.8, 1.0, 0.725, 0.4));
	G4VisAttributes *plasticVetoVis = new G4VisAttributes(G4Colour(1.0, 0.855, 0.725, 0.4));
	G4VisAttributes *shieldPE_VisAttr = new G4VisAttributes(G4Colour(129 / 255., 193 / 255., 71 / 255., 0.9));
	G4VisAttributes *shieldWaterTank_VisAttr = new G4VisAttributes(G4Colour(0, 0, 1, 0.1));
	G4VisAttributes *boricAcidCoverBoxVisAttr = new G4VisAttributes(G4Colour(0., 0., 1., 0.3)); // su-yeon
	G4VisAttributes *boricAcidCoverVisAttr = new G4VisAttributes(G4Colour(1, 0., 0., 0.7));		// su-yeon

	visAir->SetForceSolid(true);
	visRock->SetForceSolid(true);
	visCavern->SetForceSolid(true);
	leadVis->SetForceSolid(true);
	stainlessVis->SetForceSolid(true);
	aluminiumVis->SetForceSolid(true);
	CuShieldVisAttr->SetForceSolid(true);
	boricAcidVisAttr->SetForceSolid(true);
	boricAcidCaseVisAttr->SetForceSolid(true);	   // su-yeon
	boricAcidCoverBoxVisAttr->SetForceSolid(true); // su-yeon
	boricAcidCoverVisAttr->SetForceSolid(true);	   // su-yeon
	plasticScintVis->SetForceSolid(true);
	plasticVetoVis->SetForceSolid(true);
	shieldPE_VisAttr->SetForceSolid(true);
	shieldWaterTank_VisAttr->SetForceSolid(true);

	////////////////////////////////////////////////////
	// Primitive values retrived from database
	////////////////////////////////////////////////////
	G4double rockshell_thickness = db["rockshell_thickness"];
	G4double airbuffer_radius = db["airbuffer_radius"];
	G4double airbuffer_height = db["airbuffer_height"];

	G4double lead_shield_thickness = db["lead_shield_thickness"];
	G4double lead_housing_thickness = db["lead_housing_thickness"];

	G4double boricacid_thickness = db["boricacid_thickness"];
	G4double boricacidcase_thickness = db["boricacidcase_thickness"]; // su-yeon

	G4double plastic_veto_thickness = db["plastic_veto_thickness"];

	G4double PE_shield_thickness = db["PE_shield_thickness"];
	G4double thin_lead_shield_thickness = db["thin_lead_shield_thickness"];

	G4double nShield_hatGapFromLead = db["nShield_hatGapFromLead"];
	G4double nShield_GapFromCeiling = db["nShield_GapFromCeiling"];

	G4double cavern_sphere_radius = db["cavern_sphere_radius"];
	G4double rock_shell_radius = db["rock_shell_radius"];

	G4LogicalVolume *retvalLV = nullptr;

	////////////////////////////////////////////////////
	// Composite values calculated from primitive values (Temporary)
	////////////////////////////////////////////////////
	G4ThreeVector roomTlate = G4ThreeVector(-room_dist_x, room_dist_y, 0);

	// Barrel shield
	// const int nPStype = 2;
	G4double lead_shield_height = // Lead shield height
		airbuffer_height + thin_lead_shield_thickness + lead_shield_thickness;
	G4double lead_shield_halfsize = // Lead shield  half x,y
		airbuffer_radius + thin_lead_shield_thickness + lead_shield_thickness;

	G4double PS_housing_height = // Barrel size
		lead_shield_height + lead_housing_thickness * 2 +
		boricacid_thickness + PE_shield_thickness;

	G4double PS_housing_halfsize = // Barrel size
		longPSlength + profile_thickness * 2.5 + plastic_veto_thickness;

	// Rock
	G4double rock_floor_thickness = pitBox_z / 2.;

	// World
	G4double floor_spacing = rock_shell_radius;
	G4double world_size = rock_shell_radius * sqrt(2);
	G4double world_size_Z = rock_shell_radius * sqrt(2) + rock_floor_thickness + floor_spacing;

	////////////////////////////////////////////////////
	// Make the world
	////////////////////////////////////////////////////
	G4Box *worldSolid = new G4Box("World_solid", world_size, world_size, world_size_Z / 2.);
	G4LogicalVolume *logiWorld = new G4LogicalVolume(worldSolid, _air, "logiWorld");
	world_phys = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logiWorld,
								   "physWorld", nullptr, false, OverlapCheck);

	////////////////////////////////////////////////////
	// Build the rock geometry (Hemisphere and Floor)
	////////////////////////////////////////////////////
	G4Tubs *floorSolid = new G4Tubs("Rock_Floor", 0, rock_shell_radius, rock_floor_thickness, 0, 360 * deg);
	G4LogicalVolume *logiFloor = new G4LogicalVolume(floorSolid, _rock, "logiFloor");
	logiFloor->SetVisAttributes(visRock);

	G4Sphere *solidRock = new G4Sphere("Rock_Solid", 0, rock_shell_radius, 0, 360 * deg, 0, 90 * deg);
	G4LogicalVolume *logiRock = new G4LogicalVolume(solidRock, _rock, "logiRock", 0, 0, 0);
	logiRock->SetVisAttributes(visRock);
	// G4LogicalVolume *logiVirtualRock   = new G4LogicalVolume(solidRock, _rock, "logiVirtualRock", 0, 0, 0);
	// G4Sphere *solidVirtualRock = new G4Sphere("VirtualRock_Solid", 0, 15 * m, 0, 360 * deg, 0, 90 * deg);
	// G4LogicalVolume *logiRock = new G4LogicalVolume(solidVirtualRock, _rock, "logiRock", 0, 0, 0);

	G4Box *RockBox = new G4Box("RockBox",  // for rock gamma simulation
								HatInnerX + waterhousing_thickness * 2 + watertank_thickness + rockshell_thickness,
								HatInnerY + waterhousing_thickness * 2 + watertank_thickness + rockshell_thickness,
								(HatInnerZ + waterhousing_thickness * 2 + watertank_top_thickness + PMTroom_thickness) / 2.
									+ PS_housing_height / 2. + rockshell_thickness + ovc_gap / 2. + sst_zsize_half
									+ nShield_hatGapFromLead + nShield_GapFromCeiling);
	G4LogicalVolume *logiRockShell = new G4LogicalVolume(RockBox, _rock, "logiRockShell");
	logiRockShell->SetVisAttributes(visRock);

	////////////////////////////////////////////
	// Make Shield Solid and Logical volumes
	////////////////////////////////////////////
	G4Box *IDspaceBox = new G4Box("IDspace_Box",
								  airbuffer_radius - boricacid_thickness,
								  airbuffer_radius - boricacid_thickness,
								  airbuffer_height / 2. - boricacid_thickness / 2.);

	// PE shield ----------
	G4Box *shieldPEBox = new G4Box("IPEShield_Box",
								   lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness,
								   lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness,
								   PS_housing_height / 2.);
	G4VSolid *shieldPESolid = new G4SubtractionSolid("IPEShield_Solid", shieldPEBox, IDspaceBox, 0,
													 G4ThreeVector(0, 0, shieldPEBox->GetZHalfLength() - IDspaceBox->GetZHalfLength()));

	G4LogicalVolume *shieldPELV = new G4LogicalVolume(shieldPESolid, _polyethylene, "IPEShield_LV");
	shieldPELV->SetVisAttributes(shieldPE_VisAttr);

	// Boric Acid shield -------
	G4Box *shieldBoricAcidBox = new G4Box("shieldBoricAcid_Box",
										  lead_shield_halfsize + lead_housing_thickness + boricacid_thickness,
										  lead_shield_halfsize + lead_housing_thickness + boricacid_thickness,
										  (lead_shield_height + boricacid_thickness) / 2. + lead_housing_thickness);
	G4VSolid *shieldBoricAcidSolid = new G4SubtractionSolid("shieldBoricAcid_Solid", shieldBoricAcidBox, IDspaceBox, 0,
															G4ThreeVector(0, 0, shieldBoricAcidBox->GetZHalfLength() - IDspaceBox->GetZHalfLength()));

	G4LogicalVolume *shieldBoricAcidLV = new G4LogicalVolume(shieldBoricAcidSolid, _BoricAcidRubber, "shieldBoricAcid_LV");
	shieldBoricAcidLV->SetVisAttributes(boricAcidVisAttr);

	// Lead shield ----------
	G4Box *leadShieldBox = new G4Box("OuterVetoShield_Box",
									 lead_shield_halfsize, lead_shield_halfsize, lead_shield_height / 2.);
	G4VSolid *leadShieldSolid = new G4SubtractionSolid("LeadShield_Solid", leadShieldBox, IDspaceBox, 0,
													   G4ThreeVector(0, 0, leadShieldBox->GetZHalfLength() + lead_housing_thickness - IDspaceBox->GetZHalfLength()));

	G4LogicalVolume *leadShieldLV = new G4LogicalVolume(leadShieldSolid, _lead, "OuterVetoShield_LV");
	leadShieldLV->SetVisAttributes(leadVis);

	G4Box *leadHousingBox = new G4Box("OuterVetoHousing_Box",
									  lead_shield_halfsize + lead_housing_thickness,
									  lead_shield_halfsize + lead_housing_thickness,
									  lead_shield_height / 2. + lead_housing_thickness);
	G4VSolid *leadHousingSolid = new G4SubtractionSolid("LeadHousing_Solid", leadHousingBox, IDspaceBox, 0,
														G4ThreeVector(0, 0, leadHousingBox->GetZHalfLength() - IDspaceBox->GetZHalfLength()));

	G4LogicalVolume *leadHousingLV = new G4LogicalVolume(leadHousingSolid, _stainless, "OuterVetoHousing_LV");
	leadHousingLV->SetVisAttributes(stainlessVis);

	// Thin lead shield -------
	// Has been dicided to remove copper shield. Aug.2022.
	// Lead shield divided two part.
	// (For low 214Bi contaminated lead(20 cm) and low 210Pb contaimnated lead(5 cm).)
	G4Box *ThinLeadShieldBox = new G4Box("ThinLead_Box",
										 airbuffer_radius + thin_lead_shield_thickness,
										 airbuffer_radius + thin_lead_shield_thickness,
										 (airbuffer_height + thin_lead_shield_thickness) / 2.);
	G4VSolid *ThinLeadShieldSolid = new G4SubtractionSolid("ThinLead_Solid",
														   ThinLeadShieldBox, IDspaceBox, 0,
														   G4ThreeVector(0, 0,
																		 ThinLeadShieldBox->GetZHalfLength() + lead_housing_thickness - IDspaceBox->GetZHalfLength()));
	G4LogicalVolume *ThinLeadShieldLV = new G4LogicalVolume(ThinLeadShieldSolid, _lead, "ThinLeadShield_LV");
	ThinLeadShieldLV->SetVisAttributes(CuShieldVisAttr);

	// Boric Acid ----------
	G4Box *boricAcidBox = new G4Box("BoricAcid_Box",
									airbuffer_radius, airbuffer_radius, airbuffer_height / 2.);
	G4VSolid *boricAcidSolid = new G4SubtractionSolid("BoricAcid_Solid", boricAcidBox, IDspaceBox, 0,
													  G4ThreeVector(0, 0,
																	boricAcidBox->GetZHalfLength() + lead_housing_thickness - IDspaceBox->GetZHalfLength()));

	// G4LogicalVolume *boricAcidLV = new G4LogicalVolume(boricAcidSolid, _BoricAcidRubber, "BoricAcid_LV"); // su-yeon
	// G4LogicalVolume *boricAcidLV = new G4LogicalVolume(boricAcidSolid, _BoricAcidPowder, "BoricAcid_LV"); // su-yeon
	G4LogicalVolume *boricAcidLV = new G4LogicalVolume(boricAcidSolid, _air, "BoricAcid_LV"); // su-yeon
	boricAcidLV->SetVisAttributes(boricAcidVisAttr);


	////////////////////////////////////////////
	// Muon Veto Housing(Plastic scintillator)
	////////////////////////////////////////////
	// veto housing
	G4Box *plasticVetoHousing1Box = new G4Box("PlasticVetoHousing1_Box",
											  PS_housing_halfsize, PS_housing_halfsize, PS_housing_height / 2.);
	G4VSolid *plasticVetoHousingSolid = new G4SubtractionSolid("PlasticVetoHousing_Solid", plasticVetoHousing1Box, IDspaceBox, 0,
															   G4ThreeVector(0, 0, plasticVetoHousing1Box->GetZHalfLength() - IDspaceBox->GetZHalfLength()));
	G4LogicalVolume *plasticVetoHousing1LV = new G4LogicalVolume(plasticVetoHousingSolid, _air, "PlasticVetoHousing1_LV");
	plasticVetoHousing1LV->SetVisAttributes(G4VisAttributes::Invisible);

	G4Box *plasticVetoHousing2Box = new G4Box("PlasticVetoHousing2_Box",
											  bottom_veto_housingX / 2., bottom_veto_housingY / 2., plastic_veto_thickness / 2.);
	G4LogicalVolume *plasticVetoHousing2LV = new G4LogicalVolume(plasticVetoHousing2Box, _air, "PlasticVetoHousing2_LV");
	plasticVetoHousing2LV->SetVisAttributes(G4VisAttributes::Invisible);

	////////////////////////////////////////////
	// Muon Veto (Water Cerenkov)
	////////////////////////////////////////////
	// HAT Air ----------
	G4Box *shieldHatAirBox = new G4Box("HatAir_Box", HatInnerX, HatInnerY, HatInnerZ / 2.);
	G4Box *shieldHatAirBox1 = new G4Box("HatAir_Box1",
										HatInnerX + waterhousing_thickness, HatInnerY + waterhousing_thickness, HatInnerZ + waterhousing_thickness);
	G4Box *shieldHatAirBox2 = new G4Box("HatAir_Box2",
										HatInnerX + waterhousing_thickness + tyvek_thickness, HatInnerY + waterhousing_thickness + tyvek_thickness, HatInnerZ + waterhousing_thickness + tyvek_thickness);

	// Hat Water Cerenkov Housing  -------------------
	G4Box *shieldWaterHousingBox = new G4Box("WaterHousing_Box",
											 HatInnerX + waterhousing_thickness * 2 + watertank_thickness,
											 HatInnerY + waterhousing_thickness * 2 + watertank_thickness,
											 (HatInnerZ + waterhousing_thickness * 2 + watertank_top_thickness + PMTroom_thickness) / 2.);
	G4VSolid *HatWaterHousingSolid = new G4SubtractionSolid("HatWaterHousing_Solid", shieldWaterHousingBox, shieldHatAirBox, 0,
															G4ThreeVector(0, 0, -shieldWaterHousingBox->GetZHalfLength() + shieldHatAirBox->GetZHalfLength()));
	G4LogicalVolume *shieldWaterHousingLV = new G4LogicalVolume(HatWaterHousingSolid, _stainless, "HatWaterHousing_LV");
	shieldWaterHousingLV->SetVisAttributes(stainlessVis);

	// Hat reflector --------------------------------
	G4Box *shieldWaterReflectorBox = new G4Box("WCReflector_Box",
											   HatInnerX + waterhousing_thickness + watertank_thickness + tyvek_thickness,
											   HatInnerY + waterhousing_thickness + watertank_thickness + tyvek_thickness,
											   (HatInnerZ + watertank_top_thickness + PMTroom_thickness + tyvek_thickness*2) / 2.);
	G4VSolid *HatWaterReflectorSolid = new G4SubtractionSolid("WCReflector_Solid", shieldWaterReflectorBox, shieldHatAirBox1, 0,
														G4ThreeVector(0, 0, -shieldWaterReflectorBox->GetZHalfLength()));
	G4LogicalVolume *shieldWaterReflectorLV = new G4LogicalVolume(HatWaterReflectorSolid, _tyvek, "HatWaterReflector_LV");
	shieldWaterReflectorLV->SetVisAttributes("G4_WHITE");

	// Air in Water Cerenkov Housing ----------------
	G4Box *shieldWaterAirBox = new G4Box("WCAir_Box",
										 HatInnerX + waterhousing_thickness + watertank_thickness,
										 HatInnerY + waterhousing_thickness + watertank_thickness,
										 (HatInnerZ + watertank_top_thickness + PMTroom_thickness) / 2.);
	G4VSolid *HatWaterAirSolid = new G4SubtractionSolid("WCAir_Solid", shieldWaterAirBox, shieldHatAirBox2, 0,
														G4ThreeVector(0, 0, -shieldWaterAirBox->GetZHalfLength()));
	G4LogicalVolume *shieldWaterTankAirLV = new G4LogicalVolume(HatWaterAirSolid, _air, "HatWaterTankAir_LV");
	shieldWaterTankAirLV->SetVisAttributes(G4VisAttributes::Invisible);

	// Hat Water Tank --------------------------------
	G4Box *shieldWaterTankBox = new G4Box("WaterTank_Box",
										  HatInnerX + waterhousing_thickness + watertank_thickness,
										  HatInnerY + waterhousing_thickness + watertank_thickness,
										  (HatInnerZ + watertank_top_thickness) / 2.);
	G4VSolid *HatWaterTankSolid = new G4SubtractionSolid("HatWaterTank_Solid", shieldWaterTankBox, shieldHatAirBox2, 0,
														 G4ThreeVector(0, 0, -shieldWaterTankBox->GetZHalfLength()));
	G4LogicalVolume *shieldWaterTankLV = new G4LogicalVolume(HatWaterTankSolid, _water, "HatWaterTank_LV");
	shieldWaterTankLV->SetVisAttributes(shieldWaterTank_VisAttr);

	// PMT ROOM -------------------------------------
	/*
	G4Box *shieldPMTroomBox = new G4Box("PMTroom_Box",
			HatInnerX + waterhousing_thickness + watertank_thickness,
			HatInnerY + waterhousing_thickness + watertank_thickness,
			PMTroom_thickness/2);
	G4LogicalVolume *shieldPMTroomLV = new G4LogicalVolume(shieldPMTroomBox,_air,"HatPMTroom_LV");
	*/

	// HAT H-beam ----------------
	G4Box *HatBeamHousingInBox = new G4Box("HatBeamHousing_Box",
										   HatInnerX - HatHBeam_size, HatInnerY - HatHBeam_size, (HatInnerZ - HatHBeam_size) / 2.);
	G4VSolid *HatBeamHousingSolid = new G4SubtractionSolid("HatBeamHousingSolid",
														   shieldHatAirBox, HatBeamHousingInBox, 0, G4ThreeVector(0, 0, -HatHBeam_size / 2.));
	G4LogicalVolume *HatBeamHousingLV = new G4LogicalVolume(HatBeamHousingSolid, _air, "HatBeamHousing_LV");
	HatBeamHousingLV->SetVisAttributes(visAir);

	G4Box *HatBeamLong1Box = new G4Box("HatBeamLong1_Box", HatHBeam_size / 2., HatHBeam_size / 2., HatHBeam_heightL / 2.);
	G4Box *HatBeamLong2Box = new G4Box("HatBeamLong2_Box", HatHBeam_size / 2., HatHBeam_size / 2.,
									   HatHBeam_heightL / 2. - HatHBeam_size);
	G4Box *HatBeamShortBox = new G4Box("HatBeamShort_Box", HatHBeam_size / 2., HatHBeam_size / 2., HatHBeam_heightS / 2.);
	G4Box *HatBeamSpaceBox = new G4Box("HatBeamSpace_Box",
									   HatHBeam_size / 2. - HatHBeam_thickness, HatHBeam_size / 2., HatHBeam_heightL);
	G4VSolid *HatBeamLong1Solid = new G4SubtractionSolid("HatBeamLong1_Solid0",
														 HatBeamLong1Box, HatBeamSpaceBox, 0, {0, -HatHBeam_size / 2. - HatHBeam_in_thickness / 2., 0});
	HatBeamLong1Solid = new G4SubtractionSolid("HatBeamLong1_Solid",
											   HatBeamLong1Solid, HatBeamSpaceBox, 0, {0, HatHBeam_size / 2. + HatHBeam_in_thickness / 2., 0});
	G4VSolid *HatBeamLong2Solid = new G4SubtractionSolid("HatBeamLong2_Solid0",
														 HatBeamLong2Box, HatBeamSpaceBox, 0, {0, -HatHBeam_size / 2. - HatHBeam_in_thickness / 2., 0});
	HatBeamLong2Solid = new G4SubtractionSolid("HatBeamLong2_Solid",
											   HatBeamLong2Solid, HatBeamSpaceBox, 0, {0, HatHBeam_size / 2. + HatHBeam_in_thickness / 2., 0});
	G4VSolid *HatBeamShortSolid = new G4SubtractionSolid("HatBeamShort_Solid0",
														 HatBeamShortBox, HatBeamSpaceBox, 0, {0, -HatHBeam_size / 2. - HatHBeam_in_thickness / 2., 0});
	HatBeamShortSolid = new G4SubtractionSolid("HatBeamShort_Solid",
											   HatBeamShortSolid, HatBeamSpaceBox, 0, {0, HatHBeam_size / 2. + HatHBeam_in_thickness / 2., 0});
	G4LogicalVolume *HatBeamLong1LV = new G4LogicalVolume(HatBeamLong1Solid, _iron2, "HatBeamLong1_LV");
	G4LogicalVolume *HatBeamLong2LV = new G4LogicalVolume(HatBeamLong2Solid, _iron2, "HatBeamLong2_LV");
	G4LogicalVolume *HatBeamShortLV = new G4LogicalVolume(HatBeamShortSolid, _iron2, "HatBeamShort_LV");
	HatBeamLong1LV->SetVisAttributes(ironVis);
	HatBeamLong2LV->SetVisAttributes(ironVis);
	HatBeamShortLV->SetVisAttributes(ironVis);

	// HAT Aluminium plate -------
	G4Box *HatAlPlateInBox = new G4Box("HatAlPlateIn_Box",
									   HatInnerX - HatHBeam_size - HatAlPlate_thickness,
									   HatInnerY - HatHBeam_size - HatAlPlate_thickness,
									   (HatInnerZ - HatHBeam_size - HatAlPlate_thickness) / 2.);
	G4VSolid *HatAlPlateSolid = new G4SubtractionSolid("HatAlPlateSolid",
													   HatBeamHousingInBox, HatAlPlateInBox, 0, G4ThreeVector(0, 0, -HatAlPlate_thickness / 2.));

	G4LogicalVolume *HatAlPlateLV = new G4LogicalVolume(HatAlPlateSolid, _aluminium, "HatAlPlate_LV");
	HatAlPlateLV->SetVisAttributes(aluminiumVis);

	// HAT Boric Acid ------------
	G4Box *shieldHatSpaceBox = new G4Box("HatBoric_Box",
										 HatInnerX - HatHBeam_size - HatAlPlate_thickness - boricacid_thickness,
										 HatInnerY - HatHBeam_size - HatAlPlate_thickness - boricacid_thickness,
										 (HatInnerZ - HatHBeam_size - HatAlPlate_thickness - boricacid_thickness) / 2.);
	G4VSolid *HatBoricAcidSolid = new G4SubtractionSolid("HatBoricAcid_Solid", HatAlPlateInBox, shieldHatSpaceBox, 0,
														 G4ThreeVector(0, 0, -boricacid_thickness / 2.));

	G4LogicalVolume *shieldHatBoricLV = new G4LogicalVolume(HatBoricAcidSolid, _BoricAcidRubber, "HatBoric_LV");
	shieldHatBoricLV->SetVisAttributes(boricAcidVisAttr);

	// Detector supporting H-beam
	G4Box *DetHbeamHousingBox = new G4Box("DetHbeamHousingBox",
										  shieldHatSpaceBox->GetXHalfLength(), shieldHatSpaceBox->GetYHalfLength(), shieldHatSpaceBox->GetZHalfLength() - DetHbeam_size / 2.);
	G4LogicalVolume *DetHbeamHousingLV = new G4LogicalVolume(DetHbeamHousingBox, _air, "DetHbeamHousing_LV");
	DetHbeamHousingLV->SetVisAttributes(G4VisAttributes::Invisible);

	G4Box *DetHbeamHBox = new G4Box("DetHbeamH_Box", DetHbeam_size / 2., DetHbeam_size / 2., shieldHatSpaceBox->GetYHalfLength());
	G4Box *DetHbeamVBox = new G4Box("DetHbeamV_Box", DetHbeam_size / 2., DetHbeam_size / 2., PS_housing_height / 2.);
	G4Box *DetHbeamBox = new G4Box("DetHbeam_Box", DetHbeam_size / 2., DetHbeam_size / 2., DetHbeam_height / 2.);
	G4Box *DetHbeamSBox = new G4Box("DetHbeamS_Box", DetHbeam_size / 2., DetHbeam_size / 2., DetHbeamS_height / 2.);
	G4Box *DetHbeamMBox = new G4Box("DetHbeamM_Box", DetHbeam_size / 2., DetHbeam_size / 2., DetHbeamS_height / 4.);
	G4Box *DetHbeamSpaceBox = new G4Box("DetHbeamSpace_Box",
										DetHbeam_size / 2. - DetHbeam_thickness / 2., DetHbeam_size / 2. - DetHbeam_thickness, shieldHatSpaceBox->GetYHalfLength() * 2);

	G4VSolid *DetHbeamSolid = new G4SubtractionSolid("DetHbeam1_Solid", DetHbeamBox, DetHbeamSpaceBox, 0,
													 G4ThreeVector(-DetHbeam_size / 2., 0, 0));
	DetHbeamSolid = new G4SubtractionSolid("DetHbeam_Solid", DetHbeamSolid, DetHbeamSpaceBox, 0,
										   G4ThreeVector(DetHbeam_size / 2., 0, 0));
	G4LogicalVolume *DetHbeamLV = new G4LogicalVolume(DetHbeamSolid, _iron2, "DetHbeam_LV");
	DetHbeamLV->SetVisAttributes(ironVis);

	G4VSolid *DetHbeamSSolid = new G4SubtractionSolid("DetHbeam2_Solid", DetHbeamSBox, DetHbeamSpaceBox, 0,
													  G4ThreeVector(-DetHbeam_size / 2., 0, 0));
	DetHbeamSSolid = new G4SubtractionSolid("DetHbeamS_Solid", DetHbeamSSolid, DetHbeamSpaceBox, 0,
											G4ThreeVector(DetHbeam_size / 2., 0, 0));
	G4LogicalVolume *DetHbeamSLV = new G4LogicalVolume(DetHbeamSSolid, _iron2, "DetHbeamS_LV");
	DetHbeamSLV->SetVisAttributes(ironVis);

	G4VSolid *DetHbeamMSolid = new G4SubtractionSolid("DetHbeam3_Solid", DetHbeamMBox, DetHbeamSpaceBox, 0,
													  G4ThreeVector(-DetHbeam_size / 2., 0, 0));
	DetHbeamMSolid = new G4SubtractionSolid("DetHbeamM_Solid", DetHbeamMSolid, DetHbeamSpaceBox, 0,
											G4ThreeVector(DetHbeam_size / 2., 0, 0));
	G4LogicalVolume *DetHbeamMLV = new G4LogicalVolume(DetHbeamMSolid, _iron2, "DetHbeamM_LV");
	DetHbeamMLV->SetVisAttributes(ironVis);

	G4VSolid *DetHbeamHSolid = new G4SubtractionSolid("DetHbeam4_Solid", DetHbeamHBox, DetHbeamSpaceBox, 0,
													  G4ThreeVector(-DetHbeam_size / 2., 0, 0));
	DetHbeamHSolid = new G4SubtractionSolid("DetHbeamH_Solid", DetHbeamHSolid, DetHbeamSpaceBox, 0,
											G4ThreeVector(DetHbeam_size / 2., 0, 0));
	G4LogicalVolume *DetHbeamHLV = new G4LogicalVolume(DetHbeamHSolid, _iron2, "DetHbeamH_LV");
	DetHbeamHLV->SetVisAttributes(ironVis);

	G4VSolid *DetHbeamVSolid = new G4SubtractionSolid("DetHbeam5_Solid", DetHbeamVBox, DetHbeamSpaceBox, 0,
													  G4ThreeVector(-DetHbeam_size / 2., 0, 0));
	DetHbeamVSolid = new G4SubtractionSolid("DetHbeamV_Solid", DetHbeamVSolid, DetHbeamSpaceBox, 0,
											G4ThreeVector(DetHbeam_size / 2., 0, 0));
	G4LogicalVolume *DetHbeamVLV = new G4LogicalVolume(DetHbeamVSolid, _iron2, "DetHbeamV_LV");
	DetHbeamVLV->SetVisAttributes(ironVis);

	//////////////////////////////////////////
	// Cavern and Rock geometry options
	//////////////////////////////////////////
	// G4VSolid *cavern_solid = nullptr;
	G4LogicalVolume *logiCavern = nullptr;
	G4VPhysicalVolume *physCavern = nullptr;

	G4Sphere *neutronmodeCavern = nullptr;
	G4Sphere *hemiSphereCavern = nullptr;
	G4VSolid *realisticCavern = nullptr;

	G4Box *cavernLoafBox = nullptr;
	G4Tubs *cavernMainPizza = nullptr;
	G4Tubs *cavernSubPizza = nullptr;

	switch (whichSimType)
	{
		case kRockGammaMode: {// for rock gamma simulation
			fRockPhysical = new G4PVPlacement(NULL, {0,0,0}, logiRockShell, 
											"physRock", logiWorld, false, 0, OverlapCheck);
			G4Box *rockgammaCavern = new G4Box("cavern_solid",
								HatInnerX + waterhousing_thickness * 2 + watertank_thickness,
								HatInnerY + waterhousing_thickness * 2 + watertank_thickness,
								(HatInnerZ + waterhousing_thickness * 2 + watertank_top_thickness + PMTroom_thickness) / 2.
									+ PS_housing_height / 2. + ovc_gap / 2. + sst_zsize_half + nShield_hatGapFromLead + nShield_GapFromCeiling);
			logiCavern = new G4LogicalVolume(rockgammaCavern, _air, "logiCavern");
			physCavern = new G4PVPlacement(nullptr, {0, 0, 0}, logiCavern, "physCavern", logiRockShell, false, 0, OverlapCheck);
			logiCavern->SetVisAttributes(visCavern);

			break;}

		case kNeutronMode: // for neutron simulation  (There isn't rock in neutron mode)
			neutronmodeCavern = new G4Sphere("cavern_solid", 0, cavern_nMode_radius, 0., 360. * deg, 0, 180. * deg);
			// cavern_solid = neutronmodeCavern;
			logiCavern = new G4LogicalVolume(neutronmodeCavern, _air, "logiCavern");
			logiCavern->SetVisAttributes(G4VisAttributes::Invisible);
			physCavern = new G4PVPlacement(nullptr, {}, logiCavern, "physCavern", logiWorld, false, 0, OverlapCheck);
			break;

		case kIdealMode: // for muon simulation with ideal geometry (hemi-sphere cavern)
			// Rock ---
			fRockPhysical = new G4PVPlacement(NULL, {0, 0, 0}, logiRock, "physRock", logiWorld, false, 0, OverlapCheck);
			fFloorPhysical = new G4PVPlacement(nullptr, {0, 0, -rock_floor_thickness}, 
								logiFloor, "physFloor", logiWorld, false, 0, OverlapCheck);
			// AmoreEventAction::SetPrimSkew(G4ThreeVector(0,0,0));

			// Cavern ---
			hemiSphereCavern = new G4Sphere("cavern_solid", 0, cavern_sphere_radius, 0, 360 * deg, 0, 90 * deg);
			// cavern_solid = hemiSphereCavern;
			logiCavern = new G4LogicalVolume(hemiSphereCavern, _air, "logiCavern");
			logiCavern->SetVisAttributes(visCavern);
			physCavern = new G4PVPlacement(nullptr, {0, 0, 0},
								logiCavern, "physCavern", logiRock, false, OverlapCheck);
			break;
		case kRealMode: // for general simulation with realistic geometry
		{
			// Rock ---
			fRockPhysical = new G4PVPlacement(NULL, {0, 0, 0}, logiRock,
											"physRock", logiWorld, false, 0, OverlapCheck);
			fFloorPhysical = new G4PVPlacement(nullptr, {0, 0, -rock_floor_thickness}, logiFloor,
										   	"physFloor", logiWorld, false, 0, OverlapCheck);
			// AmoreEventAction::SetPrimSkew(G4ThreeVector(0,0,0));

			// Cavern ---
			G4double cavern_loaf_width = (cavern_subpizza_radius 
									+ (cavern_pizza_radius - cavern_subpizza_radius) 
									* std::sin(cavern_pizza_angle_real / 2.)) * 2.;
			G4double cavern_loaf_totalheight = cavern_loaf_height;

			cavernLoafBox = new G4Box("CavernLoafBox", 
									cavern_loaf_thickness / 2., 
									cavern_loaf_width / 2., 
									cavern_loaf_totalheight / 2.);
			cavernMainPizza = new G4Tubs("CavernMainPizza", 0, cavern_pizza_radius, 
									cavern_pizza_thickness / 2., 
									0, 
									cavern_pizza_angle_real + cavern_pizza_angle_tol * 2.);
			cavernSubPizza = new G4Tubs("CavernSubPizza", 0, cavern_subpizza_radius,
										cavern_pizza_thickness / 2., 
										0,
										cavern_pizza_angle_real + cavern_pizza_angle_tol * 2.);

			G4RotationMatrix *pizzaAlignRotMtx = new G4RotationMatrix();
			pizzaAlignRotMtx->rotateY(90 * deg);
			pizzaAlignRotMtx->rotateZ(cavern_subpizza_angle / 2. + cavern_pizza_angle_tol);

			G4RotationMatrix *pizzaPart1RotMtx = new G4RotationMatrix();
			pizzaPart1RotMtx->rotateZ(cavern_subpizza_angle + cavern_pizza_angle_tol);

			G4RotationMatrix *pizzaPart2RotMtx = new G4RotationMatrix();
			pizzaPart2RotMtx->rotateZ(-cavern_subpizza_angle - cavern_pizza_angle_tol);

			G4ThreeVector pizzaTlateInLoaf = G4ThreeVector(0, 0,
										-cavern_loaf_totalheight / 2. - cavern_totalpizza_dist);

			G4UnionSolid *pizzaUnionStage1 = new G4UnionSolid("PizzaUnion1_Solid",
													cavernMainPizza, cavernSubPizza, pizzaPart1RotMtx,
													G4ThreeVector(cavern_pizza_radius - cavern_subpizza_radius, 0, 0)
														.rotateZ(cavern_pizza_angle_tol));
			G4UnionSolid *pizzaUnionStage2 = new G4UnionSolid("PizzaUnion2_Solid",
													pizzaUnionStage1, cavernSubPizza, pizzaPart2RotMtx,
													G4ThreeVector(cavern_pizza_radius - cavern_subpizza_radius, 0, 0)
														.rotateZ(cavern_pizza_angle_real + cavern_pizza_angle_tol));
			G4SubtractionSolid *pizzaUnionStage3 = new G4SubtractionSolid("PizzaUnion3_Solid",
													pizzaUnionStage2, cavernLoafBox, 
													pizzaAlignRotMtx, {0, 0, 0});

			G4UnionSolid *cavernFinalSolid = new G4UnionSolid("CavernFinalSolid",
													cavernLoafBox, pizzaUnionStage3, 
													pizzaAlignRotMtx, pizzaTlateInLoaf);
			// cavern_solid = cavernFinalSolid;
			realisticCavern = cavernFinalSolid;
			logiCavern = new G4LogicalVolume(realisticCavern, _air, "logiCavern");
			logiCavern->SetVisAttributes(visCavern);
			physCavern = new G4PVPlacement(nullptr, roomTlate + G4ThreeVector(0, 0, cavernLoafBox->GetZHalfLength()),
								logiCavern, "physCavern", logiRock, false, OverlapCheck);
			break;
		}
		default:
			G4cout << "Error: Invalid simulation type" << G4endl;
		fCavernPhysical = physCavern;
	}

	// Upper part positioning --------------------------
	G4PVPlacement *shieldWaterHousingPV = nullptr;
	// G4PVPlacement *shieldWaterTankAirPV = nullptr;
	G4PVPlacement *shieldWaterTankPV = nullptr;
	G4PVPlacement *HatAlPlatePV = nullptr;
	G4PVPlacement *shieldHatBoricPV = nullptr;
	G4PVPlacement *HatBeamHousingPV = nullptr;
	G4PVPlacement *DetHbeamHousingPV = nullptr;
	// G4PVPlacement *DetHbeamHPV[4] = {nullptr};
	// G4PVPlacement *DetHbeamVPV[8] = {nullptr};

	G4PVPlacement *vetoHousing1PV = nullptr;
	// G4PVPlacement *vetoHousing2PV = nullptr;
	
	G4PVPlacement *shieldPEPV = nullptr;

	G4ThreeVector hatTlate = G4ThreeVector(0);
	G4ThreeVector barrelTlate = G4ThreeVector(0);
	G4ThreeVector bottomTlate = G4ThreeVector(0);

	switch(whichSimType)
	{
		case kRockGammaMode:
			hatTlate = G4ThreeVector(0, 0, nShield_hatGapFromLead);
			barrelTlate = G4ThreeVector(0, 0, -plasticVetoHousing1Box->GetZHalfLength() - nShield_GapFromCeiling);
			bottomTlate = G4ThreeVector(0, 0, -plasticVetoHousing1Box->GetZHalfLength() * 2 - nShield_GapFromCeiling 
												- HBeam_size - plasticVetoHousing2Box->GetZHalfLength());
			break;
		case kNeutronMode:
			hatTlate = G4ThreeVector(0, 0, nShield_hatGapFromLead);
			barrelTlate = G4ThreeVector(0, 0, -plasticVetoHousing1Box->GetZHalfLength() - nShield_GapFromCeiling);
			bottomTlate = G4ThreeVector(0, 0, -plasticVetoHousing1Box->GetZHalfLength() * 2 - nShield_GapFromCeiling 
												- HBeam_size - plasticVetoHousing2Box->GetZHalfLength());
			break;
		case kIdealMode:
			hatTlate = G4ThreeVector(0, 0, plasticVetoHousing1Box->GetZHalfLength() * 2 + nShield_hatGapFromLead + nShield_GapFromCeiling);
			barrelTlate = G4ThreeVector(0, 0, plasticVetoHousing1Box->GetZHalfLength());
			bottomTlate = G4ThreeVector(0, 0, rock_floor_thickness - HBeam_size - plasticVetoHousing2Box->GetZHalfLength());
			// bottomTlate = G4ThreeVector(0, 0, - HBeam_size - plasticVetoHousing2Box->GetZHalfLength());
			break;
		case kRealMode:
			hatTlate = -roomTlate + G4ThreeVector(0, 0, -cavernLoafBox->GetZHalfLength() + plasticVetoHousing1Box->GetZHalfLength() * 2 
															+ nShield_hatGapFromLead + nShield_GapFromCeiling);
			barrelTlate = -roomTlate + G4ThreeVector(0, 0, -cavernLoafBox->GetZHalfLength() + plasticVetoHousing1Box->GetZHalfLength());
			// bottomTlate = G4ThreeVector(room_dist_x, -room_dist_y,
			bottomTlate = G4ThreeVector(0, 0, rock_floor_thickness - HBeam_size - plasticVetoHousing2Box->GetZHalfLength());
			// bottomTlate = -roomTlate + G4ThreeVector(0, 0, -cavernLoafBox->GetZHalfLength() - HBeam_size - plasticVetoHousing2Box->GetZHalfLength());
			break;
		default:
			G4cout << "Error: Invalid simulation type" << G4endl;
	}

	// Hat Water cerenkov detector ----------
	shieldWaterHousingPV = new G4PVPlacement(nullptr, 
									hatTlate + G4ThreeVector(0,0,shieldWaterHousingBox->GetZHalfLength()+solidBooleanTol),
									shieldWaterHousingLV, "HatWaterHousing_PV", logiCavern, false, 0, OverlapCheck);
	new G4PVPlacement(nullptr, {0, 0, 0},
			shieldWaterReflectorLV, "HatWaterReflector_PV", shieldWaterHousingLV, false, 0, OverlapCheck);
	new G4PVPlacement(nullptr, {0, 0, 0},
			shieldWaterTankAirLV, "HatWaterTankAir_PV", shieldWaterReflectorLV, false, 0, OverlapCheck);
	shieldWaterTankPV = new G4PVPlacement(nullptr, {0, 0,
								-shieldWaterHousingBox->GetZHalfLength() + waterhousing_thickness 
									+ shieldWaterTankBox->GetZHalfLength()},
								shieldWaterTankLV, "HatWaterTank_PV", shieldWaterTankAirLV, false, 0, OverlapCheck);
	// Hat Aluminium plate -------------------
	HatAlPlatePV = new G4PVPlacement(nullptr,
							hatTlate + G4ThreeVector(0, 0, HatBeamHousingInBox->GetZHalfLength()),
							HatAlPlateLV, "HatAlPlate_PV", logiCavern, false, 0, OverlapCheck);
	// Hat Boric Acid ------------------------
	shieldHatBoricPV = new G4PVPlacement(nullptr,
							hatTlate + G4ThreeVector(0, 0, HatAlPlateInBox->GetZHalfLength()),
							shieldHatBoricLV, "HatBoricAcid_PV", logiCavern, false, 0, OverlapCheck);
	
	// PS veto housing ------------------------
	vetoHousing1PV = new G4PVPlacement(nullptr, barrelTlate, 
							plasticVetoHousing1LV, "PlasticVetoHousing_PV", logiCavern, false, 0, OverlapCheck);
	// PE shield ------------------------------
	shieldPEPV = new G4PVPlacement(nullptr, {0, 0, 0},
									shieldPELV, "PEShield_PV", plasticVetoHousing1LV, false, 0, OverlapCheck);

	switch(whichSimType)
	{
		case kRockGammaMode:
			fEnable_Gantry = false;
			break;
		case kNeutronMode:
			new G4PVPlacement(nullptr, bottomTlate,
								plasticVetoHousing2LV, "PlasticVetoHousing2_PV", logiCavern, false, 0, OverlapCheck);
			break;
		case kIdealMode:
		case kRealMode:
			new G4PVPlacement(nullptr, bottomTlate,
								plasticVetoHousing2LV, "PlasticVetoHousing2_PV", logiFloor, false, 0, OverlapCheck);
			break;
		default:
			break;
	}


	// H-beams and detector supporting H-beams
	if (fEnable_Gantry)
	{
		// Hat H-beam ----------------------------
		HatBeamHousingPV = new G4PVPlacement(nullptr,
									hatTlate + G4ThreeVector(0, 0, shieldHatAirBox->GetZHalfLength()),
									HatBeamHousingLV, "HatBeamHousing_PV", logiCavern, false, 0, OverlapCheck);

		double beamdist_y = (HatInnerY * 2 - HatHBeam_size) / 5.;
		double beamdist_x = (HatInnerX * 2 - HatHBeam_size) / 5.;
		for (int ih = 0; ih < 6; ih++)
		{
			new G4PVPlacement(nullptr, 
					{-HatInnerX + HatHBeam_size / 2., 
					-HatInnerY + HatHBeam_size / 2. + ih * beamdist_y, 
					HatInnerZ / 2. - HatHBeam_size - HatHBeam_heightS / 2.},
					HatBeamShortLV, "HatBeamShort_PV", HatBeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(nullptr, 
					{HatInnerX - HatHBeam_size / 2., 
					-HatInnerY + HatHBeam_size / 2. + ih * beamdist_y, 
					HatInnerZ / 2. - HatHBeam_size - HatHBeam_heightS / 2.},
					HatBeamShortLV, "HatBeamShort_PV", HatBeamHousingLV, false, 0, OverlapCheck);
			if (0 < ih && ih < 5)
			{
				new G4PVPlacement(nullptr, 
						{-HatInnerX + HatHBeam_size / 2. + ih * beamdist_x, 
						-HatInnerY + HatHBeam_size / 2., 
						HatInnerZ / 2. - HatHBeam_size - HatHBeam_heightS / 2.},
						HatBeamShortLV, "HatBeamShort_PV", HatBeamHousingLV, false, 0, OverlapCheck);
				new G4PVPlacement(nullptr, 
						{-HatInnerX + HatHBeam_size / 2. + ih * beamdist_x, 
						HatInnerY - HatHBeam_size / 2., 
						HatInnerZ / 2. - HatHBeam_size - HatHBeam_heightS / 2.},
						HatBeamShortLV, "HatBeamShort_PV", HatBeamHousingLV, false, 0, OverlapCheck);
			}
		}

		G4RotationMatrix *hatbeamRotMtx = new G4RotationMatrix();
		hatbeamRotMtx->rotateY(90 * deg);
		new G4PVPlacement(G4Transform3D(*hatbeamRotMtx, G4ThreeVector(
									0, 
									-HatInnerY + HatHBeam_size / 2., 
									HatInnerZ / 2. - HatHBeam_size / 2.)),
					HatBeamLong1LV, "HatBeamLong1_PV", HatBeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*hatbeamRotMtx, G4ThreeVector(
									0, 
									HatInnerY - HatHBeam_size / 2., 
									HatInnerZ / 2. - HatHBeam_size / 2.)),
					HatBeamLong1LV, "HatBeamLong1_PV", HatBeamHousingLV, false, 0, OverlapCheck);
			hatbeamRotMtx->rotateZ(90 * deg);
		new G4PVPlacement(G4Transform3D(*hatbeamRotMtx, G4ThreeVector(
									-HatInnerX + HatHBeam_size / 2., 
									0, 
									HatInnerZ / 2. - HatHBeam_size / 2.)),
					HatBeamLong2LV, "HatBeamLong2_PV", HatBeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*hatbeamRotMtx, G4ThreeVector(
									HatInnerX - HatHBeam_size / 2., 
									0, 
									HatInnerZ / 2. - HatHBeam_size / 2.)),
					HatBeamLong2LV, "HatBeamLong2_PV", HatBeamHousingLV, false, 0, OverlapCheck);

		// Detector supporting H-beam
		DetHbeamHousingPV = new G4PVPlacement(nullptr, 
										hatTlate + G4ThreeVector(0, 0, shieldHatSpaceBox->GetZHalfLength() + DetHbeam_size / 2.),
										DetHbeamHousingLV, "DetHbeamHousing_PV", logiCavern, false, 0, OverlapCheck);

		G4RotationMatrix *detbeamRotMtx = new G4RotationMatrix();
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									-DetHbeam_height / 2. - DetHbeam_size / 2., 
									-DetHbeam_height / 2. + DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. + solidBooleanTol)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									-DetHbeam_height / 2. - DetHbeam_size / 2., 
									DetHbeam_height / 2. - DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. + solidBooleanTol)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									DetHbeam_height / 2. + DetHbeam_size / 2., 
									-DetHbeam_height / 2. + DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. + solidBooleanTol)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									DetHbeam_height / 2. + DetHbeam_size / 2., 
									DetHbeam_height / 2. - DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. + solidBooleanTol)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);

		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									-DetHbeamS_height / 4. - DetHbeam_size / 2., 
									-DetHbeamS_height / 4. + DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. + solidBooleanTol + DetHbeam_height / 2. - DetHbeamS_height / 2.)),
					DetHbeamSLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									-DetHbeamS_height / 4. - DetHbeam_size / 2., 
									DetHbeamS_height / 4. - DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. + solidBooleanTol + DetHbeam_height / 2. - DetHbeamS_height / 2.)),
					DetHbeamSLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									DetHbeamS_height / 4. + DetHbeam_size / 2., 
									-DetHbeamS_height / 4. + DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. + solidBooleanTol + DetHbeam_height / 2. - DetHbeamS_height / 2.)),
					DetHbeamSLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									DetHbeamS_height / 4. + DetHbeam_size / 2., 
									DetHbeamS_height / 4. - DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. + solidBooleanTol + DetHbeam_height / 2. - DetHbeamS_height / 2.)),
					DetHbeamSLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);

		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, hatTlate
												+ G4ThreeVector(
													-DetHbeam_height / 2. - DetHbeam_size / 2.,
													-DetHbeamHousingBox->GetYHalfLength() + DetHbeam_size / 2.,
													-DetHbeamVBox->GetZHalfLength())
												+ G4ThreeVector(0, 0,
													-DetHbeamHousingBox->GetZHalfLength() - DetHbeam_size/2. - HBeam_size
													+ shieldHatSpaceBox->GetZHalfLength())),
										DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, hatTlate 
												+ G4ThreeVector(
													-DetHbeam_height / 2. + DetHbeam_size / 2.,
													DetHbeamHousingBox->GetYHalfLength() - DetHbeam_size / 2.,
													-DetHbeamVBox->GetZHalfLength())
												+ G4ThreeVector(0, 0,
													-DetHbeamHousingBox->GetZHalfLength() - DetHbeam_size/2. - HBeam_size
													+ shieldHatSpaceBox->GetZHalfLength())),
										DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, hatTlate 
												+ G4ThreeVector(
													+DetHbeam_height / 2. + DetHbeam_size / 2.,
													-DetHbeamHousingBox->GetYHalfLength() + DetHbeam_size / 2.,
													-DetHbeamVBox->GetZHalfLength())
												+ G4ThreeVector(0, 0,
													-DetHbeamHousingBox->GetZHalfLength() - DetHbeam_size/2. - HBeam_size
													+ shieldHatSpaceBox->GetZHalfLength())),
										DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, hatTlate
												+ G4ThreeVector(
														+DetHbeam_height / 2. - DetHbeam_size / 2.,
													DetHbeamHousingBox->GetYHalfLength() - DetHbeam_size / 2.,
													-DetHbeamVBox->GetZHalfLength())
												+ G4ThreeVector(0, 0,
													-DetHbeamHousingBox->GetZHalfLength() - DetHbeam_size/2. - HBeam_size
													+ shieldHatSpaceBox->GetZHalfLength())),
										DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, hatTlate
												+ G4ThreeVector(
													-DetHbeam_height / 2. - DetHbeam_size / 2.,
													-DetHbeamHousingBox->GetYHalfLength() + DetHbeam_size / 2. 
														- DetHbeamS_height,
													-DetHbeamVBox->GetZHalfLength())
												+ G4ThreeVector(0, 0,
													-DetHbeamHousingBox->GetZHalfLength() - DetHbeam_size/2. - HBeam_size
													+ shieldHatSpaceBox->GetZHalfLength())),
										DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, hatTlate
												+ G4ThreeVector(
													-DetHbeam_height / 2. + DetHbeam_size / 2.,
													DetHbeamHousingBox->GetYHalfLength() - DetHbeam_size / 2. 
														+ DetHbeamS_height,
													-DetHbeamVBox->GetZHalfLength())
												+ G4ThreeVector(0, 0,
													-DetHbeamHousingBox->GetZHalfLength() - DetHbeam_size/2. - HBeam_size
													+ shieldHatSpaceBox->GetZHalfLength())),
										DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, hatTlate
												+ G4ThreeVector(
													+DetHbeam_height / 2. + DetHbeam_size / 2.,
													-DetHbeamHousingBox->GetYHalfLength() + DetHbeam_size / 2. 
														- DetHbeamS_height,
													-DetHbeamVBox->GetZHalfLength())
												+ G4ThreeVector(0, 0,
													-DetHbeamHousingBox->GetZHalfLength() - DetHbeam_size/2. - HBeam_size
													+ shieldHatSpaceBox->GetZHalfLength())),
										DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, hatTlate
												+ G4ThreeVector( 
													+DetHbeam_height / 2. - DetHbeam_size / 2.,
													DetHbeamHousingBox->GetYHalfLength() - DetHbeam_size / 2. 
														+ DetHbeamS_height,
													-DetHbeamVBox->GetZHalfLength())
												+ G4ThreeVector(0, 0,
													-DetHbeamHousingBox->GetZHalfLength() - DetHbeam_size/2. - HBeam_size
													+ shieldHatSpaceBox->GetZHalfLength())),
									DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, OverlapCheck);

		detbeamRotMtx->rotateX(90 * deg);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									-DetHbeam_height / 2. - DetHbeam_size / 2., 
									0,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. 
										+ solidBooleanTol + DetHbeam_height / 2. + DetHbeam_size / 2.)),
						DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									DetHbeam_height / 2. + DetHbeam_size / 2., 
									0,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. 
										+ solidBooleanTol + DetHbeam_height / 2. + DetHbeam_size / 2.)),
						DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);

		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									-DetHbeamS_height / 4. - DetHbeam_size / 2., 
									0,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. 
										+ solidBooleanTol + DetHbeam_height / 2. - DetHbeamS_height - DetHbeam_size / 2.)),
						DetHbeamMLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									DetHbeamS_height / 4. + DetHbeam_size / 2., 
									0,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. 
										+ solidBooleanTol + DetHbeam_height / 2. - DetHbeamS_height - DetHbeam_size / 2.)),
					DetHbeamMLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);

		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, hatTlate
												+ G4ThreeVector(
													-DetHbeam_height / 2. - DetHbeam_size / 2., 
													0, 
													+nShield_hatGapFromLead + DetHbeam_size / 2.)
												+ G4ThreeVector(0, 0,
													-DetHbeamHousingBox->GetZHalfLength() - DetHbeam_size/2. 
													+ shieldHatSpaceBox->GetZHalfLength() - nShield_hatGapFromLead)),
										DetHbeamHLV, "HatBeamH_PV", logiCavern, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, hatTlate
												+ G4ThreeVector(
													-DetHbeam_height / 2. + DetHbeam_size / 2., 
													0, 
													+nShield_hatGapFromLead + DetHbeam_size / 2.)
												+ G4ThreeVector(0, 0,
													-DetHbeamHousingBox->GetZHalfLength() - DetHbeam_size/2. 
													+ shieldHatSpaceBox->GetZHalfLength() - nShield_hatGapFromLead)),
										DetHbeamHLV, "HatBeamH_PV", logiCavern, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, hatTlate
												+ G4ThreeVector(
													+DetHbeam_height / 2. + DetHbeam_size / 2., 
													0, 
													+nShield_hatGapFromLead + DetHbeam_size / 2.)
												+ G4ThreeVector(0, 0,
													-DetHbeamHousingBox->GetZHalfLength() - DetHbeam_size/2. 
													+ shieldHatSpaceBox->GetZHalfLength() - nShield_hatGapFromLead)),
										DetHbeamHLV, "HatBeamH_PV", logiCavern, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, hatTlate
												+ G4ThreeVector(
													+DetHbeam_height / 2. - DetHbeam_size / 2., 
													0, 
													+nShield_hatGapFromLead + DetHbeam_size / 2.)
												+ G4ThreeVector(0, 0,
													-DetHbeamHousingBox->GetZHalfLength() - DetHbeam_size/2. 
													+ shieldHatSpaceBox->GetZHalfLength() - nShield_hatGapFromLead)),
										DetHbeamHLV, "HatBeamH_PV", logiCavern, false, 0, OverlapCheck);

		detbeamRotMtx->rotateZ(90 * deg);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									0, 
									-DetHbeam_height / 2. + DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. 
										+ solidBooleanTol + DetHbeam_height / 2. + DetHbeam_size / 2.)),
						DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									0, 
									DetHbeam_height / 2. - DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. 
										+ solidBooleanTol + DetHbeam_height / 2. + DetHbeam_size / 2.)),
						DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									0, 
									DetHbeamS_height / 4. - DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. 
										+ solidBooleanTol + DetHbeam_height / 2. + DetHbeam_size / 2.)),
						DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									0, 
									-DetHbeamS_height / 4. + DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. 
										+ solidBooleanTol + DetHbeam_height / 2. + DetHbeam_size / 2.)),
						DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);

		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									0, 
									-DetHbeamS_height / 4. + DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. 
										+ solidBooleanTol + DetHbeam_height / 2. - DetHbeamS_height - DetHbeam_size / 2.)),
						DetHbeamMLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(G4Transform3D(*detbeamRotMtx, G4ThreeVector(
									0, 
									DetHbeamS_height / 4. - DetHbeam_size / 2.,
									-DetHbeamHousingBox->GetZHalfLength() + DetHbeam_height / 2. 
										+ solidBooleanTol + DetHbeam_height / 2. - DetHbeamS_height - DetHbeam_size / 2.)),
						DetHbeamMLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
	}


	//////////////////////////////////////////
	// Boric aicd powder/rubber cover surrounding OVC, su-yeon ----------
	//////////////////////////////////////////
	G4Box *boricAcidCoverBox = nullptr;
	G4Box *boricAcidCoverBox1 = nullptr;
	G4Box *boricAcidCoverBox2 = nullptr;
	G4Tubs *boricAcidCoverBoxHole = nullptr;
	G4VSolid *boricAcidCoverBoxSolid = nullptr;
	G4VSolid *boricAcidCoverBoxFin = nullptr;
	G4LogicalVolume *boricAcidCoverBoxLV = nullptr;
	// G4PVPlacement *boricAcidCoverBoxPV = nullptr;

	G4Box *boricAcidCover = nullptr;
	G4Box *boricAcidCover1 = nullptr;
	G4Box *boricAcidCover2 = nullptr;
	G4VSolid *boricAcidCoverHole = nullptr;
	G4VSolid *boricAcidCoverSolid = nullptr;
	G4VSolid *boricAcidCoverFin = nullptr;
	G4LogicalVolume *boricAcidCoverLV = nullptr;
	// G4PVPlacement *boricAcidCoverPV = nullptr;

	G4RotationMatrix *CoverHoleZ = nullptr;
	switch (whichNShieldDesign)
	{
		case kNSDesign1: // su-yeon's idea
			// Cover desing 1, boru cover series, su-yeon's idea
			// boric acid powder container
			boricAcidCoverBox1 = new G4Box("BoricAcidCoverBox1_Box",
										lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness,
										lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness,
										boricacidcase_thickness/2.);
			boricAcidCoverBox2 = new G4Box("BoricAcidCoverBox2_Box",
										airbuffer_radius - boricacid_thickness, 
										airbuffer_radius - boricacid_thickness, 
										11./2.*cm);
			boricAcidCoverBoxHole = new G4Tubs("BoricAcidCoverBox2_Tubs", 0, ss_radius + 1.*cm, 13./2.*cm, 0.0*deg, 360.0*deg);

			boricAcidCoverBoxSolid = new G4UnionSolid("BoricAcidCoverBox_Solid", boricAcidCoverBox1, boricAcidCoverBox2, 
												0, G4ThreeVector(0, 0, -11./2.*cm + 2./2.*cm));
			boricAcidCoverBoxFin = new G4SubtractionSolid("BoricAcidCoverBoxFin", boricAcidCoverBoxSolid, boricAcidCoverBoxHole, 
												0, G4ThreeVector(0, 0, -10./2.*cm));
			boricAcidCoverBoxLV = new G4LogicalVolume(boricAcidCoverBoxFin, _polyethylene, "BoricAcidCoverBox_LV");
	//				G4LogicalVolume *boricAcidCoverBoxLV = new G4LogicalVolume(boricAcidCoverBoxFin, _air, "BoricAcidCoverBox_LV");
			boricAcidCoverBoxLV -> SetVisAttributes(boricAcidCoverBoxVisAttr);

			// boric acid powder
			boricAcidCover1 = new G4Box("BoricAcidCover_Box",
									lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness - 1./2.*cm,
									lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness - 1./2.*cm,
									1./2.*cm);
			boricAcidCover2 = new G4Box("BoricAcidCover2_Box",
									airbuffer_radius - boricacid_thickness - 1./2.*cm, 
									airbuffer_radius - boricacid_thickness - 1./2.*cm, 
									10./2.*cm);
			boricAcidCoverHole = new G4Tubs("BoricAcidCover_Tubs", 0, ss_radius + 1.*cm + 1./2.*cm, 12./2.*cm, 0.0*deg, 360.0*deg);

			boricAcidCoverSolid = new G4UnionSolid("BoricAcidCover_Solid", boricAcidCover1, boricAcidCover2, 
											0, G4ThreeVector(0, 0, -10./2.*cm + 1./2.*cm));
			boricAcidCoverFin = new G4SubtractionSolid("BoricAcidCover_Solid", boricAcidCoverSolid, boricAcidCoverHole, 
											0, G4ThreeVector(0, 0, -10./2.*cm));
			boricAcidCoverLV = new G4LogicalVolume(boricAcidCoverFin, _BoricAcidPowder, "BoricAcidCover_LV");
			// G4LogicalVolume *boricAcidCoverLV = new G4LogicalVolume(boricAcidCoverFin, _BoricAcidRubber, "BoricAcidCover_LV");
			// G4LogicalVolume *boricAcidCoverLV = new G4LogicalVolume(boricAcidCoverFin, _air, "BoricAcidCover_LV");
			boricAcidCoverLV -> SetVisAttributes(boricAcidCoverVisAttr);

			// Positioning
			new G4PVPlacement(nullptr,
					barrelTlate + G4ThreeVector(0, 0, plasticVetoHousing1Box->GetZHalfLength() + boricAcidCoverBox1->GetZHalfLength()),
										// {0, 0, - nShield_GapFromCeiling + 2./2.*cm},
										boricAcidCoverBoxLV, "BoricAcidCoverBox_PV", logiCavern, false, 0, OverlapCheck);
			new G4PVPlacement(nullptr, {0, 0, 0},
									boricAcidCoverLV, "BoricAcidCover_PV", boricAcidCoverBoxLV, false, 0, OverlapCheck);
			break;
		case kNSDesign2: // Jaison's idea
			CoverHoleZ = new G4RotationMatrix();
			CoverHoleZ->rotateZ(45 * deg);

			boricAcidCover = new G4Box("BoricAcidCover_Box",
									airbuffer_radius - boricacid_thickness - 1./2.*cm, 
									airbuffer_radius - boricacid_thickness - 1./2.*cm, 
									10./2.*cm);
			boricAcidCoverHole = new G4Box("BoricAcidCover_BoxHole",
										ss_radius + 1.*cm + 1./2.*cm, ss_radius + 1.*cm + 1./2.*cm, 10./2.*cm);
			boricAcidCoverSolid = new G4SubtractionSolid("BoricAcidCover_Solid", boricAcidCover, boricAcidCoverHole, 
										CoverHoleZ, G4ThreeVector(0, 0, 0));

			boricAcidCoverLV = new G4LogicalVolume(boricAcidCoverSolid, _BoricAcidRubber, "BoricAcidCover_LV");
			boricAcidCoverLV -> SetVisAttributes(boricAcidCoverVisAttr);

			new G4PVPlacement(nullptr, 
					barrelTlate + G4ThreeVector(0, 0, plasticVetoHousing1Box->GetZHalfLength() - boricAcidCover->GetZHalfLength() + ovc_gap/2.),
						// {0, 0, - nShield_GapFromCeiling - 10./2.*cm},
									boricAcidCoverLV, "BoricAcidCover_PV", logiCavern, false, 0, OverlapCheck);
			break;
		case kNSDesign3: // su-yeon's idea
			// Cover design 3, lead side cover, su-yeon's idea

			boricAcidCoverBox = new G4Box("BoricAcidCoverBox_Box",
										airbuffer_radius - boricacid_thickness, 
										airbuffer_radius - boricacid_thickness, 
										11./2.*cm);
			boricAcidCoverBoxHole = new G4Tubs("BoricAcidCoverBox_Tubs", 0, ss_radius + 1.*cm, 11./2.*cm, 0.0*deg, 360.0*deg);
			boricAcidCoverBoxSolid = new G4SubtractionSolid("BoricAcidCoverBox_Solid", boricAcidCoverBox, boricAcidCoverBoxHole, 
												0, G4ThreeVector(0, 0, 0));
			boricAcidCoverBoxLV = new G4LogicalVolume(boricAcidCoverBoxSolid, _polyethylene, "BoricAcidCoverBox_LV");
			boricAcidCoverBoxLV -> SetVisAttributes(boricAcidCoverBoxVisAttr);

			boricAcidCover = new G4Box("BoricAcidCover_Box",
									airbuffer_radius - boricacid_thickness - 1./2.*cm, 
									airbuffer_radius - boricacid_thickness - 1./2.*cm, 
									10./2.*cm);
			boricAcidCoverHole = new G4Tubs("BoricAcidCover_Tubs", 0, ss_radius + 1.*cm + 1./2.*cm, 10./2.*cm, 0.0*deg, 360.0*deg);
			boricAcidCoverSolid = new G4SubtractionSolid("BoricAcidCover_Solid", boricAcidCover, boricAcidCoverHole, 0, G4ThreeVector(0, 0, 0));
			boricAcidCoverLV = new G4LogicalVolume(boricAcidCoverSolid, _BoricAcidRubber, "BoricAcidCover_LV");
			// boricAcidCoverLV = new G4LogicalVolume(boricAcidCoverSolid, _BoricAcidPowder, "BoricAcidCover_LV");
			//	boricAcidCoverLV = new G4LogicalVolume(boricAcidCoverSolid, _air, "BoricAcidCover_LV");
			boricAcidCoverLV -> SetVisAttributes(boricAcidCoverVisAttr);

			new G4PVPlacement(nullptr,
					barrelTlate + G4ThreeVector(0, 0, plasticVetoHousing1Box->GetZHalfLength() - boricAcidCoverBox->GetZHalfLength() + ovc_gap/2. + boricacid_thickness/2.),
										// {0, 0, - nShield_GapFromCeiling - 11./2.*cm},
										boricAcidCoverBoxLV, "BoricAcidCoverBox_PV", logiCavern, false, 0, OverlapCheck);
			new G4PVPlacement(nullptr, {0, 0, 0}, boricAcidCoverLV, "BoricAcidCover_PV", boricAcidCoverBoxLV, false, 0, OverlapCheck);
			break;
		case kNSDesign4: // su-yeon's idea
		{// boric acid powder container
			boricAcidCoverBox1 = new G4Box("BoricAcidCoverBox1_Box",
										lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness,
										lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness,
										2./2.*cm);
			boricAcidCoverBox2 = new G4Box("BoricAcidCoverBox2_Box",
										airbuffer_radius - boricacid_thickness, airbuffer_radius - boricacid_thickness, 11./2.*cm);
			boricAcidCoverBoxHole = new G4Tubs("BoricAcidCoverBox2_Tubs", 0, ss_radius + 1.*cm, 13./2.*cm, 0.0*deg, 360.0*deg);

			boricAcidCoverBoxSolid = new G4UnionSolid("BoricAcidCoverBox_Solid", boricAcidCoverBox1, boricAcidCoverBox2, 
											0, G4ThreeVector(0, 0, -11./2.*cm + 2./2.*cm));
			boricAcidCoverBoxFin = new G4SubtractionSolid("BoricAcidCoverBoxFin", boricAcidCoverBoxSolid, boricAcidCoverBoxHole, 
											0, G4ThreeVector(0, 0, -10./2.*cm));
			boricAcidCoverBoxLV = new G4LogicalVolume(boricAcidCoverBoxFin, _polyethylene, "BoricAcidCoverBox_LV");
			//	G4LogicalVolume *boricAcidCoverBoxLV = new G4LogicalVolume(boricAcidCoverBoxFin, _air, "BoricAcidCoverBox_LV");
			boricAcidCoverBoxLV -> SetVisAttributes(boricAcidCoverBoxVisAttr);

			G4cout << "BoricAcidCoverBox2 width = " << airbuffer_radius - boricacid_thickness << G4endl;
			G4cout << "BoricAcidCoverBoxHole diameter = " << ss_radius + 1.*cm << G4endl;

			// Outer and Inner boric acid powder
			boricAcidCover1 = new G4Box("BoricAcidCover_Box", // outer boric acid powder
									lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness - 1./2.*cm,
									lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness - 1./2.*cm,
									1./2.*cm);
			G4Box *boricAcidCover1Hole = new G4Box("BoricAcidCover1Hole_Tubs",
										airbuffer_radius - boricacid_thickness - 1./2.*cm, 
										airbuffer_radius - boricacid_thickness - 1./2.*cm, 
										2./2.*cm);

			boricAcidCover2 = new G4Box("BoricAcidCover2_Box", // inner boric acid powder
									airbuffer_radius - boricacid_thickness - 1./2.*cm, 
									airbuffer_radius - boricacid_thickness - 1./2.*cm, 
									1./2.*cm);
			G4Tubs *boricAcidCover2Hole = new G4Tubs("BoricAcidCover2Hole_Tubs", 0, ss_radius + 1.*cm + 1./2.*cm, 2./2.*cm, 0.0*deg, 360.0*deg);


			G4VSolid *OuterBoricAcidCoverFin = new G4SubtractionSolid("OuterBoricAcidCover_Solid", 
														boricAcidCover1, boricAcidCover1Hole, 0, G4ThreeVector(0, 0, 0));
			G4LogicalVolume *OuterBoricAcidCoverLV = new G4LogicalVolume(OuterBoricAcidCoverFin, _BoricAcidPowder, "OuterBoricAcidCover_LV");
			// G4LogicalVolume *boricAcidCoverLV = new G4LogicalVolume(boricAcidCoverFin, _air, "BoricAcidCover_LV");
			OuterBoricAcidCoverLV -> SetVisAttributes(boricAcidCoverVisAttr);

			G4VSolid *InnerBoricAcidCoverFin = new G4SubtractionSolid("InnerBoricAcidCover_Solid", 
														boricAcidCover2, boricAcidCover2Hole, 0, G4ThreeVector(0, 0, 0));
			G4LogicalVolume *InnerBoricAcidCoverLV = new G4LogicalVolume(InnerBoricAcidCoverFin, _BoricAcidPowder, "BoricAcidCover_LV");
			//	G4LogicalVolume *boricAcidCoverLV = new G4LogicalVolume(boricAcidCoverFin, _air, "BoricAcidCover_LV");
			InnerBoricAcidCoverLV -> SetVisAttributes(boricAcidCoverVisAttr);

			new G4PVPlacement(nullptr, 
					barrelTlate + G4ThreeVector(0, 0, plasticVetoHousing1Box->GetZHalfLength() + boricAcidCoverBox1->GetZHalfLength()),
										// {0, 0, - nShield_GapFromCeiling + 2./2.*cm},
					boricAcidCoverBoxLV, "BoricAcidCoverBox_PV", logiCavern, false, 0, OverlapCheck);
			new G4PVPlacement(nullptr, {0, 0, 0},
									OuterBoricAcidCoverLV, "OuterBoricAcidCover_PV", boricAcidCoverBoxLV, false, 0,  OverlapCheck);
			new G4PVPlacement(nullptr, {0, 0, -11.*cm + 2*cm},
										InnerBoricAcidCoverLV, "InnerBoricAcidCover_PV", boricAcidCoverBoxLV, false, 0, OverlapCheck);
			break;
		}
		case kNSDesign5: // su-yeon's idea
		{
			// boric acid powder container
			boricAcidCoverBox1 = new G4Box("BoricAcidCoverBox1_Box",
										  lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness,
										  lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness,
										  2. / 2. * cm);
			boricAcidCoverBox2 = new G4Box("BoricAcidCoverBox2_Box",
										  airbuffer_radius - boricacid_thickness, airbuffer_radius - boricacid_thickness, 11. / 2. * cm);
			boricAcidCoverBoxHole = new G4Tubs("BoricAcidCoverBox2_Tubs", 0, ss_radius + 1. * cm, 13. / 2. * cm, 0.0 * deg, 360.0 * deg);

			boricAcidCoverBoxSolid = new G4UnionSolid("BoricAcidCoverBox_Solid", boricAcidCoverBox1, boricAcidCoverBox2, 
											0, G4ThreeVector(0, 0, -11. / 2. * cm + 2. / 2. * cm));
			boricAcidCoverBoxFin = new G4SubtractionSolid("BoricAcidCoverBoxFin", boricAcidCoverBoxSolid, boricAcidCoverBoxHole, 
											0, G4ThreeVector(0, 0, -10. / 2. * cm));
			boricAcidCoverBoxLV = new G4LogicalVolume(boricAcidCoverBoxFin, _polyethylene, "BoricAcidCoverBox_LV");
			// G4LogicalVolume *boricAcidCoverBoxLV = new G4LogicalVolume(boricAcidCoverBoxFin, _air, "BoricAcidCoverBox_LV");
			boricAcidCoverBoxLV->SetVisAttributes(boricAcidCoverBoxVisAttr);

			G4cout << "BoricAcidCoverBox2 width = " << airbuffer_radius - boricacid_thickness << G4endl;
			G4cout << "BoricAcidCoverBoxHole diameter = " << ss_radius + 1. * cm << G4endl;

			// Inner boric acid powder
			boricAcidCover2 = new G4Box("BoricAcidCover2_Box", // inner boric acid powder
									   airbuffer_radius - boricacid_thickness - 1. / 2. * cm, 
									   airbuffer_radius - boricacid_thickness - 1. / 2. * cm, 
									   1. / 2. * cm);
			G4Tubs *boricAcidCover2Hole = new G4Tubs("BoricAcidCover2Hole_Tubs", 0, ss_radius + 1. * cm + 1. / 2. * cm, 2. / 2. * cm, 
												0.0 * deg, 360.0 * deg);
			G4VSolid *InnerBoricAcidCoverFin = new G4SubtractionSolid("InnerBoricAcidCover_Solid", boricAcidCover2, boricAcidCover2Hole, 
														0, G4ThreeVector(0, 0, 0));
			G4LogicalVolume *InnerBoricAcidCoverLV = new G4LogicalVolume(InnerBoricAcidCoverFin, _BoricAcidPowder, "BoricAcidCover_LV");
			// G4LogicalVolume *boricAcidCoverLV = new G4LogicalVolume(boricAcidCoverFin, _air, "BoricAcidCover_LV");
			InnerBoricAcidCoverLV->SetVisAttributes(boricAcidCoverVisAttr);

			new G4PVPlacement(nullptr,
										barrelTlate + G4ThreeVector(0, 0, plasticVetoHousing1Box->GetZHalfLength() + boricAcidCoverBox1->GetZHalfLength()),
										// {0, 0, -nShield_GapFromCeiling + boricacidcase_thickness / 2.},
										boricAcidCoverBoxLV, "BoricAcidCoverBox_PV", logiCavern, false, 0, OverlapCheck);
			new G4PVPlacement(nullptr, {0, 0, -11. * cm + boricacidcase_thickness},
										InnerBoricAcidCoverLV, "InnerBoricAcidCover_PV", boricAcidCoverBoxLV, false, 0, OverlapCheck);
			break;
		}
		default:
			break;
	}

	// Outer BoricAcid Shield -------------
	G4PVPlacement *shieldBoricAcidPV = new G4PVPlacement(nullptr,
														 G4ThreeVector(0, 0, PE_shield_thickness / 2.),
														 shieldBoricAcidLV, "BoricAcidShield_PV", shieldPELV, false, 0, OverlapCheck);

	// Lead shield -----------------
	f200_VetoMaterialPhysical = new G4PVPlacement(nullptr,
												  {0, 0, boricacid_thickness / 2.},
												  leadHousingLV, "LeadShield_PV", shieldBoricAcidLV, false, 0, OverlapCheck);

	G4VPhysicalVolume *BufferMother = new G4PVPlacement(nullptr, {},
														leadShieldLV, "LeadShield_PV", leadHousingLV, false, 0, OverlapCheck);

	// Thin Lead shield -----------
	new G4PVPlacement(nullptr, {0, 0, lead_shield_thickness / 2.},
					  ThinLeadShieldLV, "ThinLeadShield_PV", leadShieldLV, false, 0, OverlapCheck);

/*
	// Boric Acid ----------
	new G4PVPlacement(nullptr, {0, 0, thin_lead_shield_thickness / 2.},
					  boricAcidLV, "InnerBoricAcid_PV", ThinLeadShieldLV, false, 0, OverlapCheck);
					  */


	//////////////////////////////////////////
	// H-beam and Pit option
	//////////////////////////////////////////
	// WC supporting H-beam housing----------
	G4double dist[7] = {2150, 380, 3800, 565, 465, 1680, 1680};
	G4double len[4] = {1192, 2208, 2218, 882};
	G4Box *HbeamHousingOut = nullptr;
	G4Box *HbeamHousingIn = nullptr;
	G4VSolid *HbeamHousing = nullptr;
	G4Box *HbeamHousingCut = new G4Box("hbeamHousingCut", sst_xsize_half, sst_ysize_half, PS_housing_height);

	if (whichSimType == kRockGammaMode){
		HbeamHousingOut = new G4Box("hbeamhousingOut_Box", 
								HatInnerX + waterhousing_thickness * 2 + watertank_thickness,
								HatInnerX + waterhousing_thickness * 2 + watertank_thickness,
								(PS_housing_height + HBeam_size + nShield_GapFromCeiling) /2.);
		HbeamHousingIn = new G4Box("hbeamhousingIn_Box",
								HatInnerX + waterhousing_thickness * 2 + watertank_thickness - HBeam_size,
								HatInnerX + waterhousing_thickness * 2 + watertank_thickness - HBeam_size,
								(PS_housing_height + nShield_GapFromCeiling) /2.);
		HbeamHousing = new G4SubtractionSolid("hbeamhousing_Solid", HbeamHousingOut, HbeamHousingIn, 0,
														// G4ThreeVector(0,0, - HbeamHousingOut->GetZHalfLength() - HBeam_size));
														G4ThreeVector(0, 0, -HbeamHousingOut->GetZHalfLength() + HbeamHousingIn->GetZHalfLength()));
		HbeamHousing = new G4SubtractionSolid("hbeamhousing_Solid", HbeamHousing, HbeamHousingCut, 0,
											G4ThreeVector(0, 0, 0));
	} else{
		HbeamHousingOut = new G4Box("hbeamhousingOut_Box",
									   HBeam_housingX / 2., HBeam_housingY / 2., (PS_housing_height + HBeam_size + nShield_GapFromCeiling) / 2.);
		HbeamHousingIn = new G4Box("hbeamhousingIn_Box",
									  HBeam_housingX / 2. - HBeam_size,
									  HBeam_housingY / 2. - HBeam_size,
									  (PS_housing_height + nShield_GapFromCeiling) / 2.);
		HbeamHousing = new G4SubtractionSolid("hbeamhousing_Solid", HbeamHousingOut, HbeamHousingIn, 0,
														// G4ThreeVector(0,0, - HbeamHousingOut->GetZHalfLength() - HBeam_size));
														G4ThreeVector(0, 0, -HbeamHousingOut->GetZHalfLength() + HbeamHousingIn->GetZHalfLength()));
		HbeamHousing = new G4SubtractionSolid("hbeamhousing_Solid", HbeamHousing, HbeamHousingCut, 0,
											G4ThreeVector(HBeam_housingDist, 0, 0));
	}
	G4LogicalVolume *HbeamHousingLV = new G4LogicalVolume(HbeamHousing, _air, "HbeamHousing_LV");
	HbeamHousingLV->SetVisAttributes(G4VisAttributes::Invisible);


	G4LogicalVolume *HbeamBotLV = nullptr;
	G4LogicalVolume *PitLV = nullptr;
	if (fEnable_Gantry)
	{
		// H-beam Bottom ----------
		HbeamBotLV = MakeHBeam1(_iron1);
		HbeamBotLV->SetVisAttributes(ironVis);

		// Pit ----------
		PitLV = MakePit(pitBox_x, pitBox_z, _rebar, HBeam_size);
		PitLV->SetVisAttributes(ironVis);

		// WC supporting H-beam -------------------
		G4Box *beamCut = new G4Box("beamCut", HBeam_size, HBeam_size / 2. - HBeam_thickness, HBeam_housingX);
		G4RotationMatrix *beamRotMtx = new G4RotationMatrix();

		// H-beam vertical
		G4Box *HbeamLong = new G4Box("HbeamLong", HBeam_size / 2., HBeam_size / 2., PS_housing_height / 2.);
		G4VSolid *HbeamV = new G4SubtractionSolid("HbeamV1", HbeamLong, beamCut, 0,
												  {HBeam_size + HBeam_thickness, 0, 0});
		HbeamV = new G4SubtractionSolid("HbeamV", HbeamV, beamCut, 0,
										{-HBeam_size - HBeam_thickness, 0, 0});
		G4LogicalVolume *HbeamV_LV = new G4LogicalVolume(HbeamV, _iron2, "HbeamV_LV");
		HbeamV_LV->SetVisAttributes(ironVis);

		// H-beam horizontal
		beamRotMtx->rotateX(90 * deg);
		G4Box *HbeamLong2 = new G4Box("HbeamLong2", HBeam_size / 2., HBeam_housingY / 2., HBeam_size / 2.);
		G4VSolid *HbeamH = new G4SubtractionSolid("HbeamH1", HbeamLong2, beamCut,
												  G4Transform3D(*beamRotMtx, {HBeam_size + HBeam_thickness, 0, 0}));
		HbeamH = new G4SubtractionSolid("HbeamH", HbeamH, beamCut,
										G4Transform3D(*beamRotMtx, {-HBeam_size + -HBeam_thickness, 0, 0}));
		G4LogicalVolume *HbeamH_LV = new G4LogicalVolume(HbeamH, _iron2, "HbeamH_LV");
		HbeamH_LV->SetVisAttributes(ironVis);

		// Small H-beams
		G4Box *HbeamShort[7];
		G4VSolid *HbeamSH[7];
		G4LogicalVolume *HbeamSH_LV[7];
		beamRotMtx->rotateZ(90 * deg);
		for (int i = 0; i < 7; i++)
		{
			HbeamShort[i] = new G4Box(Form("HbeamShort%d", i), dist[i] / 2., HBeam_size / 2., HBeam_size / 2.);
			HbeamSH[i] = new G4SubtractionSolid(Form("HbeamSH1%d", i), HbeamShort[i], beamCut,
												G4Transform3D(*beamRotMtx, {0, HBeam_size + HBeam_thickness, 0}));
			HbeamSH[i] = new G4SubtractionSolid(Form("HbeamSH1%d", i), HbeamSH[i], beamCut,
												G4Transform3D(*beamRotMtx, {0, -HBeam_size - HBeam_thickness, 0}));
			HbeamSH_LV[i] = new G4LogicalVolume(HbeamSH[i], _iron2, Form("HbeamSH%d_LV", i));
			HbeamSH_LV[i]->SetVisAttributes(ironVis);
		}

		////////////////////////////////////////////////////////////////
		//// H-beam, pit positioning
		////////////////////////////////////////////////////////////////
		// WC supporting H-beam
		G4ThreeVector vbeamPos = G4ThreeVector(HBeam_housingX / 2. - HBeam_size / 2., HBeam_housingY / 2. - HBeam_size / 2., -HBeam_size / 2.);
		G4ThreeVector hbeamPos = G4ThreeVector(HBeam_housingX / 2. - HBeam_size / 2., 0, HbeamHousingOut->GetZHalfLength() - HBeam_size / 2.);

		new G4PVPlacement(nullptr, vbeamPos, HbeamV_LV, "HbeamV_PV", HbeamHousingLV, false, 0, OverlapCheck);
		vbeamPos[1] *= -1;
		new G4PVPlacement(nullptr, vbeamPos, HbeamV_LV, "HbeamV_PV", HbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(nullptr, hbeamPos, HbeamH_LV, "HbeamH_PV", HbeamHousingLV, false, 0, OverlapCheck);
		for (int i = 0; i < 7; i++)
		{
			vbeamPos[0] -= dist[i] + HBeam_size;
			new G4PVPlacement(nullptr, vbeamPos, HbeamV_LV, "HbeamV_PV", HbeamHousingLV, false, 0, OverlapCheck);
			vbeamPos[1] *= -1;
			new G4PVPlacement(nullptr, vbeamPos, HbeamV_LV, "HbeamV_PV", HbeamHousingLV, false, 0, OverlapCheck);
			hbeamPos[0] -= dist[i] + HBeam_size;
			new G4PVPlacement(nullptr, hbeamPos, HbeamH_LV, "HbeamH_PV", HbeamHousingLV, false, 0, OverlapCheck);
		}
		for (int i = 0; i < 3; i++)
		{
			vbeamPos[1] -= len[i] + HBeam_size;
			new G4PVPlacement(nullptr, vbeamPos, HbeamV_LV, "HbeamV_PV", HbeamHousingLV, false, 0, OverlapCheck);
			vbeamPos[0] *= -1;
			new G4PVPlacement(nullptr, vbeamPos, HbeamV_LV, "HbeamV_PV", HbeamHousingLV, false, 0, OverlapCheck);
		}

		G4ThreeVector sbeamPos = G4ThreeVector(HBeam_housingX / 2., HBeam_housingY / 2. - HBeam_size / 2.,
											   HbeamHousingOut->GetZHalfLength() - HBeam_size / 2.);
		for (int i = 0; i < 7; i++)
		{
			sbeamPos[0] -= dist[i] / 2. + HBeam_size;
			new G4PVPlacement(nullptr, sbeamPos, HbeamSH_LV[i], "HbeamSH_PV", HbeamHousingLV, false, 0, OverlapCheck);
			for (int j = 0; j < 4; j++)
			{
				if ((i == 1 || i == 2 || i == 3) && (j == 0 || j == 3))
				{
					sbeamPos[1] -= len[j] / 2. + HBeam_size / 2.;
					new G4PVPlacement(nullptr, sbeamPos, HbeamSH_LV[i], "HbeamSH_PV", HbeamHousingLV, false, 0, OverlapCheck);
					sbeamPos[1] -= len[j] / 2. + HBeam_size / 2.;
					new G4PVPlacement(nullptr, sbeamPos, HbeamSH_LV[i], "HbeamSH_PV", HbeamHousingLV, false, 0, OverlapCheck);
				}
				else if (i == 2 && j == 1)
				{
					sbeamPos[1] -= len[j] + HBeam_size;
				}
				else
				{
					sbeamPos[1] -= len[j] + HBeam_size;
					new G4PVPlacement(nullptr, sbeamPos, HbeamSH_LV[i], "HbeamSH_PV", HbeamHousingLV, false, 0, OverlapCheck);
				}
			}
			sbeamPos[1] = HBeam_housingY / 2. - HBeam_size / 2.;
			sbeamPos[0] -= dist[i] / 2.;
		}

		//////////////////////////////////////////////////////////////
		G4PVPlacement *HbeamBotPV = nullptr;
		G4PVPlacement *pitPV = nullptr;
		G4PVPlacement *HbeamHousingPV = nullptr;

		switch (whichSimType )
		{
			case kRockGammaMode:
				break;
			case kNeutronMode:
				HbeamBotPV = new G4PVPlacement(nullptr, {0, 0, -plasticVetoHousing1Box->GetZHalfLength() * 2 - nShield_GapFromCeiling-HBeam_size /2.},
										HbeamBotLV, "HbeamBot_PV", logiCavern, false, 0, OverlapCheck);

			// bottomTlate = G4ThreeVector(0, 0, -plasticVetoHousing1Box->GetZHalfLength() * 2 - nShield_GapFromCeiling 
												// - HBeam_size - plasticVetoHousing2Box->GetZHalfLength());

				pitPV = new G4PVPlacement(nullptr, 
					bottomTlate + G4ThreeVector(-pitBox_x/2., 0, -pitBox_x/2.-HBeam_size),
				// {-pitBox_x/2., 0, -plasticVetoHousing1Box->GetZHalfLength()*2 - nShield_GapFromCeiling - pitBox_x/2.},
									PitLV, "pit_PV", logiCavern, false, 0, OverlapCheck);
				HbeamHousingPV = new G4PVPlacement(nullptr, {-HBeam_housingDist, 0, 
												-plasticVetoHousing1Box->GetZHalfLength()*2 - nShield_GapFromCeiling + HbeamHousingOut->GetZHalfLength()},
												HbeamHousingLV, "HbeamHousing_PV", logiCavern, false, 0, OverlapCheck);
				break;
			case kIdealMode:
				HbeamBotPV = new G4PVPlacement(nullptr, {0, 0, rock_floor_thickness - HBeam_size / 2.},
										HbeamBotLV, "HbeamBot_PV", logiFloor, false, 0, OverlapCheck);
				pitPV = new G4PVPlacement(nullptr, {-pitBox_x/2., 0, rock_floor_thickness - pitBox_z/2.},
									PitLV, "pit_PV", logiFloor, false, 0, OverlapCheck);	
				HbeamHousingPV = new G4PVPlacement(nullptr, {-HBeam_housingDist, 0, HbeamHousingOut->GetZHalfLength()},
												HbeamHousingLV, "HbeamHousing_PV", logiCavern, false, 0, OverlapCheck);
				break;
			case kRealMode:
				HbeamBotPV = new G4PVPlacement(nullptr, {0, 0, rock_floor_thickness - HBeam_size / 2.},
										HbeamBotLV, "HbeamBot_PV", logiFloor, false, 0, OverlapCheck);
				pitPV = new G4PVPlacement(nullptr, {-pitBox_x/2., 0, 0},
									PitLV, "pit_PV", logiFloor, false, 0, OverlapCheck);
				HbeamHousingPV = new G4PVPlacement(nullptr, {room_dist_x - HBeam_housingDist, -room_dist_y,
												-cavernLoafBox->GetZHalfLength() + HbeamHousingOut->GetZHalfLength()},
												HbeamHousingLV, "HbeamHousing_PV", logiCavern, false, 0, OverlapCheck);
				break;
			default:
				break;
		}

		if (fDbgMsgOn)
		{
			cout << " ======================== H-beams =========" << endl;
			cout << "  ( length: mm, weight: kg, global coordinate)" << endl;
			cout << " Bottom H-beam " << endl;
			cout << "     mass      : " << HbeamBotLV->GetMass(true, false) / kg << endl;
			if (fFloorPhysical != nullptr)
			{
				cout << "     coordinate: " << HbeamBotPV->GetTranslation() + fFloorPhysical->GetTranslation() << endl;
			}
			else if (whichSimType != kRockGammaMode){
				cout << "     coordinate: " << HbeamBotPV->GetTranslation() << endl;
				cout << " Pit " << endl;
				cout << "     mass      : " << PitLV->GetMass(true, false) / kg << endl;
			}
			if (fFloorPhysical != nullptr)
			{
				cout << "     coordinate: " << pitPV->GetTranslation() + fFloorPhysical->GetTranslation() << endl;
			}
			else if (whichSimType != kRockGammaMode){
				cout << "     coordinate: " << pitPV->GetTranslation() << endl;
				cout << " H-beam housing 1" << endl;
				cout << "     dimension                     : " << HbeamHousingOut->GetXHalfLength() * 2 << " x "
				 << HbeamHousingOut->GetYHalfLength() * 2 << " x "
				 << HbeamHousingOut->GetZHalfLength() * 2 << endl;
				cout << "     coordinaate                   : "
				 << HbeamHousingPV->GetTranslation() + physCavern->GetTranslation() << endl;
				cout << "     mass (HbeamV, HbeamH, HbeamSH): " << HbeamV_LV->GetMass(true, false) / kg << ", "
				 << HbeamH_LV->GetMass(true, false) / kg << ", "
				 << HbeamSH_LV[0]->GetMass(true, false) / kg << endl;
				cout << " H-beam housing 2" << endl;
				cout << "     dimension                 : " << shieldHatAirBox->GetXHalfLength() * 2 << " x "
				 << shieldHatAirBox->GetYHalfLength() * 2 << " x "
				 << shieldHatAirBox->GetZHalfLength() * 2 << endl;
				cout << "     coordinaate               : " << HatBeamHousingPV->GetTranslation() + physCavern->GetTranslation() << endl;
				cout << "     mass (Long1, Long2, Short): " << HatBeamLong1LV->GetMass(true, false) / kg << ", "
				 << HatBeamLong2LV->GetMass(true, false) / kg << ", "
				 << HatBeamShortLV->GetMass(true, false) / kg << endl;
				cout << " H-beam housing 3" << endl;
				cout << "     dimension                               : " << DetHbeamHousingBox->GetXHalfLength() * 2 << " x "
				 << DetHbeamHousingBox->GetYHalfLength() * 2 << " x "
				 << DetHbeamHousingBox->GetZHalfLength() * 2 << endl;
				cout << "     coordinaate                             : " << DetHbeamHousingPV->GetTranslation() + physCavern->GetTranslation() << endl;
				cout << "     mass (Long, Short, Horizontal, Vertical): " << DetHbeamLV->GetMass(true, false) / kg << ", "
				 << DetHbeamSLV->GetMass(true, false) / kg << ", "
				 << DetHbeamHLV->GetMass(true, false) / kg << ", "
				 << DetHbeamVLV->GetMass(true, false) / kg << endl;
			}
		}
	}
	else
	{ // Without H-beams
		switch(whichSimType)
		{
			case kRockGammaMode:
				new G4PVPlacement(nullptr, {0, 0, -plasticVetoHousing1Box->GetZHalfLength() * 2 - nShield_GapFromCeiling + HbeamHousingOut->GetZHalfLength()},
										HbeamHousingLV, "HbeamHousing_PV", logiCavern, false, 0, OverlapCheck);
				break;
			case kNeutronMode:
				new G4PVPlacement(nullptr, {-HBeam_housingDist, 0, 
										-plasticVetoHousing1Box->GetZHalfLength() * 2 - nShield_GapFromCeiling + HbeamHousingOut->GetZHalfLength()},
										HbeamHousingLV, "HbeamHousing_PV", logiCavern, false, 0, OverlapCheck);
				break;
			case kIdealMode:
				new G4PVPlacement(nullptr, {-HBeam_housingDist, 0, HbeamHousingOut->GetZHalfLength()},
										HbeamHousingLV, "HbeamHousing_PV", logiCavern, false, 0, OverlapCheck);
				break;
			case kRealMode:
				new G4PVPlacement(nullptr, {room_dist_x - HBeam_housingDist, -room_dist_y,
										-cavernLoafBox->GetZHalfLength() + HbeamHousingOut->GetZHalfLength()},
										HbeamHousingLV, "HbeamHousing_PV", logiCavern, false, 0, OverlapCheck);
				break;
			default:
				break;
		}
	}

	//////////////////////////////////////////
	// Additional PE
	//////////////////////////////////////////
	G4Box *addPEBox1 = nullptr;
	G4Box *addPEBox2 = nullptr;
	G4LogicalVolume *addPE1LV = nullptr;
	G4LogicalVolume *addPE2LV = nullptr;

	if (fAdditionalPE)
	{
		G4ThreeVector addPEPos = G4ThreeVector(
			HBeam_housingX / 2. - HBeam_size - dist[0] / 2.,
			HBeam_housingY / 2. - HBeam_size * 3 / 2 - len[0] / 2.,
			HbeamHousingOut->GetZHalfLength() - HBeam_size / 2.);

		G4PVPlacement *AddPE_PV[12] = {nullptr,};
		int pvid = 0;
		for (int i = 1; i < 5; i++)
		{
			addPEPos[0] -= dist[i - 1] / 2. + HBeam_size + dist[i] / 2.;
			for (int j = 0; j < 4; j++)
			{
				if (j == 0 || j == 3)
				{
					addPEBox1 = new G4Box("addPEBox",
										  dist[i] / 2. + HBeam_size / 2. - HBeam_thickness,
										  len[j] / 4. + HBeam_size / 4. - HBeam_thickness - 5 * solidBooleanTolBig,
										  HBeam_size / 2. - HBeam_thickness);

					addPE1LV = new G4LogicalVolume(addPEBox1, _polyethylene, "additionalPE_LV");
					addPE1LV->SetVisAttributes(shieldPE_VisAttr);

					if (j == 3)
					{
						addPEPos[1] -= (len[j] - HBeam_size * 3) / 4.;
					}
					addPEPos[1] -= (len[j] - HBeam_size) / 4.;
					AddPE_PV[pvid] =  new G4PVPlacement(nullptr, addPEPos, addPE1LV, "additionalPE1_PV", HbeamHousingLV, false, 0, OverlapCheck);
					pvid++;
					addPEPos[1] -= (len[j] + HBeam_size * 3) / 4.;
				}
				else
				{
					addPEBox2 = new G4Box("addPEBox",
										  dist[i] / 2. + HBeam_size / 2. - HBeam_thickness,
										  len[j] / 2. + HBeam_size / 2. - HBeam_thickness - solidBooleanTol,
										  HBeam_size / 2. - HBeam_thickness);

					addPE2LV = new G4LogicalVolume(addPEBox2, _polyethylene, "additionalPE_LV");
					addPE2LV->SetVisAttributes(shieldPE_VisAttr);

					addPEPos[1] -= len[j] / 2.;
					if (i != 2 && i != 3)
					{
						AddPE_PV[pvid] =  new G4PVPlacement(nullptr, addPEPos, addPE2LV, "additionalPE2_PV", HbeamHousingLV, false, 0, OverlapCheck);
						pvid++;
					}
					addPEPos[1] -= len[j] / 2. + HBeam_size;
				}
			}
			addPEPos[1] = HBeam_housingY / 2. - HBeam_size * 3 / 2 - len[0] / 2.;
		}
		if(whichSimType==kRockGammaMode){
			for (int ipv = 0; ipv < 12; ipv++){
				AddPE_PV[ipv]->SetTranslation(AddPE_PV[ipv]->GetTranslation() + G4ThreeVector(-HBeam_housingDist,0,0));
			}

		}
	}

	//////////////////////////////////////////
	// Additional PE Ver.2, su-yeon
	//////////////////////////////////////////
	/*
	//				Hbeam size reference
	G4double dist[7] = {2150,380,3800,565,465,1680,1680};
	G4double len[4] = {1192,2208,2218,882};
	*/
	// long one
	// G4double ParBoruBox3_x = dist[1]/2. + HBeam_size/2. - HBeam_thickness - 10./2.*mm; // ParallelPEOnly
	G4double ParBoruBox3_x = 80. / 2. * cm; // ParallelPEOnly_measuredsize
	G4double ParBoruBox3_y = len[1] + HBeam_size - HBeam_thickness - solidBooleanTol - 5. / 2. * cm;
	// G4double ParBoruBox3_z = HBeam_size/2. - HBeam_thickness; // ParallelPEOnly
	G4double ParBoruBox3_z = 20. / 2. * cm; // ParallelPEOnly_measuredsize

	// G4double ParPEBox3_x = dist[1]/2. + HBeam_size/2. - HBeam_thickness - 20./2.*mm; // ParallelPEOnly
	G4double ParPEBox3_x = (80. - 1.) / 2. * cm; // ParallelPEOnly_measuredsize
	G4double ParPEBox3_y = len[1] + HBeam_size - HBeam_thickness - solidBooleanTol - 5. / 2. * cm;
	// G4double ParPEBox3_z = HBeam_size/2. - HBeam_thickness - 10./2.*mm; // ParallelPEOnly
	G4double ParPEBox3_z = (20. - 1.) / 2. * cm; // ParallelPEOnly_measuredsize

	G4Box *ParBoruBox3 = new G4Box("ParBoruBox3", ParBoruBox3_x, ParBoruBox3_y, ParBoruBox3_z);
	G4Box *ParPEBox3 = new G4Box("ParPEBox3", ParPEBox3_x, ParPEBox3_y, ParPEBox3_z);
	G4LogicalVolume *ParBoru3LV = new G4LogicalVolume(ParBoruBox3, _BoricAcidRubber, "ParallelBoru3_LV");
	G4LogicalVolume *ParPE3LV = new G4LogicalVolume(ParPEBox3, _polyethylene, "ParallelPE3_LV");
	ParBoru3LV->SetVisAttributes(boricAcidVisAttr);
	ParPE3LV->SetVisAttributes(shieldPE_VisAttr);

	// short one
	// G4double ParBoruBox4_x = dist[4]/2. + HBeam_size/2. - HBeam_thickness - 10./2.*mm; // ParallelPEOnly
	G4double ParBoruBox4_x = 80. / 2. * cm; // ParallelPEOnly_measuredsize
	G4double ParBoruBox4_y = (len[1] + HBeam_size - HBeam_thickness - solidBooleanTol) / 2. - 50. / 2. * mm;
	// G4double ParBoruBox4_z = HBeam_size/2. - HBeam_thickness; // ParallelPEOnly
	G4double ParBoruBox4_z = 20. / 2. * cm; // ParallelPEOnly_measuredsize

	// G4double ParPEBox4_x = dist[4]/2. + HBeam_size/2. - HBeam_thickness - 10./2.*mm - 10./2.*mm; // ParallelPEOnly
	G4double ParPEBox4_x = (80. - 1.) / 2. * cm; // ParallelPEOnly_measuredsize
	G4double ParPEBox4_y = (len[1] + HBeam_size - HBeam_thickness - solidBooleanTol) / 2. - 50. / 2. * mm;
	// G4double ParPEBox4_z = HBeam_size/2. - HBeam_thickness - 10./2.*mm; // ParallelPEOnly
	G4double ParPEBox4_z = (20. - 1.) / 2. * cm; // ParallelPEOnly_measuredsize

	G4Box *ParBoruBox4 = new G4Box("ParBoruBox4", ParBoruBox4_x, ParBoruBox4_y, ParBoruBox4_z);
	G4Box *ParPEBox4 = new G4Box("ParPEBox4", ParPEBox4_x, ParPEBox4_y, ParPEBox4_z);

	G4LogicalVolume *ParBoru4LV = new G4LogicalVolume(ParBoruBox4, _BoricAcidRubber, "ParallelBoru4_LV");
	G4LogicalVolume *ParPE4LV = new G4LogicalVolume(ParPEBox4, _polyethylene, "ParallelePE4_LV");
	ParBoru4LV->SetVisAttributes(boricAcidVisAttr);
	ParPE4LV->SetVisAttributes(shieldPE_VisAttr);

	G4Box *ParBoruBox5 = new G4Box("ParBoruBox5", ParBoruBox4_x, ParBoruBox4_y, ParBoruBox4_z);
	G4Box *ParPEBox5 = new G4Box("ParPEBox5", ParPEBox4_x, ParPEBox4_y, ParPEBox4_z);

	G4LogicalVolume *ParBoru5LV = new G4LogicalVolume(ParBoruBox5, _BoricAcidRubber, "ParallelBoru5_LV");
	G4LogicalVolume *ParPE5LV = new G4LogicalVolume(ParPEBox5, _polyethylene, "ParallelePE5_LV");
	ParBoru5LV->SetVisAttributes(boricAcidVisAttr);
	ParPE5LV->SetVisAttributes(shieldPE_VisAttr);
	/*
	   cout << "ParBoruBox3_x = " << ParBoruBox3_x << endl;
	   cout << "ParBoruBox3_y = " << ParBoruBox3_y << endl;
	   cout << "ParBoruBox3_z = " << ParBoruBox3_z << endl;
	   cout << "ParBoruBox4_x = " << ParBoruBox4_x << endl;
	   cout << "ParBoruBox4_y = " << ParBoruBox4_y << endl;
	   cout << "ParBoruBox4_z = " << ParBoruBox4_z << endl;
	   */

	// Additional PE V3. (Adding PE to v2)
	G4double FNBBoruBox_x = 200. / 2. * cm;
	G4double FrontBoruBox_y = 60. / 2. * cm;
	G4double BackBoruBox_y = 80. / 2. * cm;
	G4double FNBBoruBox_z = 15. / 2. * cm;

	G4double FNBPEBox_x = 200. / 2. * cm;
	G4double FrontPEBox_y = (60. - 1.) / 2. * cm;
	G4double BackPEBox_y = (80. - 1.) / 2. * cm;
	G4double FNBPEBox_z = (15. - 1.) / 2. * cm;

	G4Box *FrontBoruBox = new G4Box("FrontBoruBox", FNBBoruBox_x, FrontBoruBox_y, FNBBoruBox_z);
	G4Box *FrontPEBox = new G4Box("FrontPEBox", FNBPEBox_x, FrontPEBox_y, FNBPEBox_z);

	G4Box *BackBoruBox = new G4Box("BackBoruBox", FNBBoruBox_x, BackBoruBox_y, FNBBoruBox_z);
	G4Box *BackPEBox = new G4Box("BackPEBox", FNBPEBox_x, BackPEBox_y, FNBPEBox_z);

	G4LogicalVolume *FrontBoruLV = new G4LogicalVolume(FrontBoruBox, _BoricAcidRubber, "FrontBoru_LV");
	G4LogicalVolume *FrontPELV = new G4LogicalVolume(FrontPEBox, _polyethylene, "FrontPE_LV");
	FrontBoruLV->SetVisAttributes(boricAcidVisAttr);
	FrontPELV->SetVisAttributes(shieldPE_VisAttr);

	G4LogicalVolume *BackBoruLV = new G4LogicalVolume(BackBoruBox, _BoricAcidRubber, "BackBoru_LV");
	G4LogicalVolume *BackPELV = new G4LogicalVolume(BackPEBox, _polyethylene, "BackPE_LV");
	BackBoruLV->SetVisAttributes(boricAcidVisAttr);
	BackPELV->SetVisAttributes(shieldPE_VisAttr);
	
	if (fAdditionalPE)
	{
		G4PVPlacement *AddPE1_PV =  new G4PVPlacement(nullptr,
						  G4ThreeVector(
							  //			2860 + 40 + 10, // ParallelPEOnly
							  2860, // ParallelPEOnly_measuredsize
							  1104 - (-len[2] / 2. + HBeam_size / 2. - HBeam_thickness - solidBooleanTol + ParPEBox3_y) + 250,
							  HbeamHousingOut->GetZHalfLength() - HBeam_size / 2.),
						  ParBoru3LV, "ParallelBoru3_PV", HbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(nullptr, G4ThreeVector(ParBoruBox3_x - ParPEBox3_x, 0, -ParBoruBox3_z + ParPEBox3_z), ParPE3LV, "ParallelPE3_PV", ParBoru3LV, false, 0, OverlapCheck);

		G4PVPlacement *AddPE2_PV =  new G4PVPlacement(nullptr,
						  G4ThreeVector(
							  //			-1220 - 5,  // ParallelPEOnly
							  -1220 - 30, // ParallelPEOnly_measuredsize
							  1104,
							  HbeamHousingOut->GetZHalfLength() - HBeam_size / 2.),
						  ParBoru4LV, "ParallelBoru4_PV", HbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(nullptr, G4ThreeVector(-ParBoruBox4_x + ParPEBox4_x, 0, -ParBoruBox4_z + ParPEBox4_z), ParPE4LV, "ParallelPE4_PV", ParBoru4LV, false, 0, OverlapCheck);

		G4PVPlacement *AddPE3_PV =  new G4PVPlacement(nullptr,
						  G4ThreeVector(
							  //			- 1220 - 5, // ParallelPEOnly
							  -1220 - 30, // ParallelPEOnly_measuredsize
							  -1104 - 300,
							  HbeamHousingOut->GetZHalfLength() - HBeam_size / 2.),
						  ParBoru5LV, "ParallelBoru5_PV", HbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(nullptr, G4ThreeVector(-ParBoruBox4_x + ParPEBox4_x, 0, -ParBoruBox4_z + ParPEBox4_z), ParPE5LV, "ParallelPE5_PV", ParBoru5LV, false, 0, OverlapCheck);

		G4PVPlacement *AddPE4_PV =  new G4PVPlacement(nullptr,
						  G4ThreeVector(740, 2000, HbeamHousingOut->GetZHalfLength() - HBeam_size / 2.),
						  FrontBoruLV, "FrontBoru_PV", HbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(nullptr, G4ThreeVector(0, FrontBoruBox_y - FrontPEBox_y, -FNBBoruBox_z + FNBPEBox_z), FrontPELV, "FrontPE_PV", FrontBoruLV, false, 0, OverlapCheck);

		G4PVPlacement *AddPE5_PV =  new G4PVPlacement(nullptr,
						  G4ThreeVector(740, -2200, HbeamHousingOut->GetZHalfLength() - HBeam_size / 2.),
						  BackBoruLV, "BackBoru_PV", HbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(nullptr, G4ThreeVector(0, -BackBoruBox_y + BackPEBox_y, -FNBBoruBox_z + FNBPEBox_z), BackPELV, "BackPE_PV", BackBoruLV, false, 0, OverlapCheck);
		
		if(whichSimType==kRockGammaMode){
			AddPE1_PV->SetTranslation(AddPE1_PV->GetTranslation() + G4ThreeVector(-HBeam_housingDist,0,0));
			AddPE2_PV->SetTranslation(AddPE2_PV->GetTranslation() + G4ThreeVector(-HBeam_housingDist,0,0));
			AddPE3_PV->SetTranslation(AddPE3_PV->GetTranslation() + G4ThreeVector(-HBeam_housingDist,0,0));
			AddPE4_PV->SetTranslation(AddPE4_PV->GetTranslation() + G4ThreeVector(-HBeam_housingDist,0,0));
			AddPE5_PV->SetTranslation(AddPE5_PV->GetTranslation() + G4ThreeVector(-HBeam_housingDist,0,0));
		}
	}

	//////////////////////////////////////////////////////
	if (fDbgMsgOn)
	{
		cout << "\n ================================ Outer shields ======" << endl;
		cout << " (length: mm, weight: kg, global coordinate)\n"
			 << endl;
		cout << " World dimension(box): "
			 << worldSolid->GetXHalfLength() * 2 << " x "
			 << worldSolid->GetYHalfLength() * 2 << " x "
			 << worldSolid->GetZHalfLength() * 2 << endl;
		switch(whichSimType)
		{
			case kRockGammaMode:
				cout << " Rock for rockgamma simulation" << endl;
				cout << "     shell mass  : " << logiRockShell->GetMass(true, false) / kg << endl;
				break;
			case kNeutronMode:
				cout << " Cavern" << endl;
				cout << "     dimension (r): " << neutronmodeCavern->GetOuterRadius() << endl;
				cout << "     coordinate   : " << physCavern->GetTranslation() << endl;
				break;
			case kIdealMode:
			case kRealMode:
				cout << " Rock floor (tube)" << endl;
				cout << "     mass             : " << logiFloor->GetMass(true, false) / kg << endl;
				cout << "     dimension (r x h): "
				 	<< floorSolid->GetOuterRadius() << " x " << floorSolid->GetZHalfLength() * 2 << endl;
				cout << "     coordinate       : " << fFloorPhysical->GetTranslation() << endl;
				cout << " Rock body (hemisphere)" << endl;
				cout << "     mass         : " << logiRock->GetMass(true, false) / kg << endl;
				cout << "     dimension (r): " << solidRock->GetOuterRadius() << endl;
				cout << "     coordinate   : " << fRockPhysical->GetTranslation() << endl;
				cout << " Cavern ( '" << AmoreDetectorConstruction::GetSimTypeName(whichSimType)
				 	<< "' option selected )" << endl;
				if (whichSimType == kIdealMode)	{
					cout << "     dimension (r): " << hemiSphereCavern->GetOuterRadius() << endl;
					cout << "     coordinate   : " << fRockPhysical->GetTranslation() + physCavern->GetTranslation() << endl;
				} else{
					cout << "     dimension (box part): "
						<< cavernLoafBox->GetXHalfLength() * 2 << " x "
						<< cavernLoafBox->GetYHalfLength() * 2 << " x "
						<< cavernLoafBox->GetZHalfLength() * 2 << endl;
				}
				break;
			default:
				break;
		}	

		// if (whichSimType!=kRockGammaMode)
		// {
			/*
			cout << " ============================ Muon Veto ======" << endl;
			cout << " PS veto housing " << endl; //.......................
			cout << "     dimension : "
				<< plasticVetoHousing1Box->GetXHalfLength()*2 << " x "
				<< plasticVetoHousing1Box->GetYHalfLength()*2 << " x "
				<< plasticVetoHousing1Box->GetZHalfLength()*2 << endl;
			cout << "     thickness : "
				<< plasticVetoHousing1Box->GetXHalfLength() - shieldPEBox->GetXHalfLength() << endl;
			cout << "     coordinate: "
				<< physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() << endl;
			cout << " PS supporter (Aluminium profile) mass: \n"
				<< "     " << plasticVetoSupporterV1LV->GetMass(true,false)/kg << "("
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
				<< plasticVetoSupporterH2Box->GetZHalfLength()*2. << ")" << endl;
			cout << " PS veto stainless flame" << endl;
			cout << "     mass       : " << plasticVetoLV[0]->GetMass(true,false)/kg << endl;
			cout << "     dimension  : "
				<< plasticVetoBox[0]->GetXHalfLength()*2 << " x "
				<< plasticVetoBox[0]->GetYHalfLength()*2 << " x "
				<< plasticVetoBox[0]->GetZHalfLength()*2 << endl;
			cout << "     coordinate : " <<  physCavern->GetTranslation() + vetoHousing1PV->GetTranslation()
				+ GetPhysicalVolumeByName("PlasticVeto0_PV")->GetTranslation() << endl;
			cout << " PS veto Aluminium plate " << endl;
			cout << "     mass      : " << aluminiumHolderLV[0]->GetMass(true,false)/kg << endl;
			cout << "     dimension : " <<  plasticScintHolderBox[0]->GetXHalfLength()*2 << " x "
				<<  plasticScintHolderBox[0]->GetYHalfLength()*2 << " x "
				<<  plasticScintHolderBox[0]->GetZHalfLength()*2 << endl;
			cout << " Plastic scintillator" << endl;
			cout << "     dimension: "
				<< plasticScintBox[0]->GetXHalfLength()*2 << " x "
				<< plasticScintBox[0]->GetYHalfLength()*2 << " x "
				<< plasticScintBox[0]->GetZHalfLength()*2 << endl;
			cout << "     mass     : " << plasticScintOLV[0]->GetMass(true,false)/kg << endl;
			cout << " PS veto detector ( # of veto: " << f200_VetoTotCNum
				<< " =>  barrel: " << nVetoZ*8 << ", bottom: " << nVetoB*2 << ")" << endl;
				*/
		cout << " ============================= Bottom shields ======" << endl;
		cout << " PE (box) " << endl; //..............................
		cout << "     mass      : " << shieldPELV->GetMass(true, false) / kg << endl;
		cout << "     dimension : "
				<< shieldPEBox->GetXHalfLength() * 2 << " x "
				<< shieldPEBox->GetYHalfLength() * 2 << " x "
				<< shieldPEBox->GetZHalfLength() * 2 << endl;
		cout << "     thickness : "
				<< shieldPEBox->GetXHalfLength() - shieldBoricAcidBox->GetXHalfLength() << endl;
		cout << "     coordinate: "
				<< physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() + shieldPEPV->GetTranslation() << endl;
		cout << " BoricAcid " << endl; //..............................
		cout << "     mass      : " << shieldBoricAcidLV->GetMass(true, false) / kg << endl;
		cout << "     dimension : "
				<< shieldBoricAcidBox->GetXHalfLength() * 2 << " x "
				<< shieldBoricAcidBox->GetYHalfLength() * 2 << " x "
				<< shieldBoricAcidBox->GetZHalfLength() * 2 << endl;
		cout << "     thickness : "
				<< shieldBoricAcidBox->GetXHalfLength() - leadShieldBox->GetXHalfLength() << endl;
		cout << "     coordinate: "
				<< physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() + shieldPEPV->GetTranslation() + shieldBoricAcidPV->GetTranslation() << endl;
		cout << " Lead shield " << endl; //..............................
		cout << "     mass      : " << leadShieldLV->GetMass(true, false) / kg << endl;
		cout << "     dimension : "
				<< leadShieldBox->GetXHalfLength() * 2 << " x "
				<< leadShieldBox->GetYHalfLength() * 2 << " x "
				<< leadShieldBox->GetZHalfLength() * 2 << endl;
		cout << "     thickness : "
				<< leadShieldBox->GetXHalfLength() - ThinLeadShieldBox->GetXHalfLength() << endl;
		cout << "     coordinate: "
				<< physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() + shieldPEPV->GetTranslation() + shieldBoricAcidPV->GetTranslation() + f200_VetoMaterialPhysical->GetTranslation() + BufferMother->GetTranslation() << endl;
		cout << " Thin lead shield " << endl; //..............................
		cout << "     mass      : " << ThinLeadShieldLV->GetMass(true, false) / kg << endl;
		cout << "     dimension : "
				<< ThinLeadShieldBox->GetXHalfLength() * 2 << " x "
				<< ThinLeadShieldBox->GetYHalfLength() * 2 << " x "
				<< ThinLeadShieldBox->GetZHalfLength() * 2 << endl;
		cout << "     thickness : "
				<< ThinLeadShieldBox->GetXHalfLength() - boricAcidBox->GetXHalfLength() << endl;
		cout << "     coordinate: "
				<< physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() + shieldPEPV->GetTranslation() + shieldBoricAcidPV->GetTranslation() + f200_VetoMaterialPhysical->GetTranslation() + BufferMother->GetTranslation() + GetPhysicalVolumeByName("ThinLeadShield_PV")->GetTranslation() << endl;
				 /*
			cout << " Inner boric acid " << endl; //..............................
			cout << "     mass      : " << boricAcidLV->GetMass(true, false) / kg << endl;
			cout << "     dimension : "
				 << boricAcidBox->GetXHalfLength() * 2 << " x "
				 << boricAcidBox->GetYHalfLength() * 2 << " x "
				 << boricAcidBox->GetZHalfLength() * 2 << endl;
			cout << "     thickness : "
				 << boricAcidBox->GetXHalfLength() - IDspaceBox->GetXHalfLength() << endl;
			cout << "     coordinate: "
				 << physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() + shieldPEPV->GetTranslation() + shieldBoricAcidPV->GetTranslation() + f200_VetoMaterialPhysical->GetTranslation() + BufferMother->GetTranslation() + GetPhysicalVolumeByName("ThinLeadShield_PV")->GetTranslation() + GetPhysicalVolumeByName("InnerBoricAcid_PV")->GetTranslation() << endl;
				 */
		cout << " Inner detector space " << endl; //..............................
		cout << "     dimension : "
				<< IDspaceBox->GetXHalfLength() * 2 << " x "
				<< IDspaceBox->GetYHalfLength() * 2 << " x "
				<< IDspaceBox->GetZHalfLength() * 2 << endl;
		cout << " ============================= Upper shields ======" << endl;
		if (fAdditionalPE)
		{
			cout << " Additional PE shield" << endl;
			cout << "     mass        : " << addPE1LV->GetMass(true, false) / kg
					<< ", " << addPE2LV->GetMass(true, false) / kg << endl;
			cout << "     dimension(1): " << addPEBox1->GetXHalfLength() * 2 << " x "
					<< addPEBox1->GetYHalfLength() * 2 << " x "
					<< addPEBox1->GetZHalfLength() * 2 << endl;
			cout << "     dimension(2): " << addPEBox2->GetXHalfLength() * 2 << " x "
					<< addPEBox2->GetYHalfLength() * 2 << " x "
					<< addPEBox2->GetZHalfLength() * 2 << endl;
			cout << "     coordinate  : " << GetPhysicalVolumeByName("additionalPE1_PV")->GetTranslation() + physCavern->GetTranslation() + GetPhysicalVolumeByName("HbeamHousing_PV")->GetTranslation() << ", "
					<< GetPhysicalVolumeByName("additionalPE2_PV")->GetTranslation() + physCavern->GetTranslation() + GetPhysicalVolumeByName("HbeamHousing_PV")->GetTranslation() << endl;
		}
		cout << " WC veto housing " << endl;
		cout << "     mass              : " << shieldWaterHousingLV->GetMass(true, false) / kg << endl;
		cout << "     dimension (out)   : "
				<< shieldWaterHousingBox->GetXHalfLength() * 2 << " x "
				<< shieldWaterHousingBox->GetYHalfLength() * 2 << " x "
				<< shieldWaterHousingBox->GetZHalfLength() * 2 << endl;
		cout << "     dimension (inner) : "
				<< shieldHatAirBox->GetXHalfLength() * 2 << " x "
				<< shieldHatAirBox->GetYHalfLength() * 2 << " x "
				<< shieldHatAirBox->GetZHalfLength() * 2 << endl;
		cout << "     thickness         : " << waterhousing_thickness << endl;
		cout << "     coordinate        : " << shieldWaterHousingPV->GetTranslation() + physCavern->GetTranslation() << endl;
		cout << " WC water" << endl;
		cout << "     mass      : " << shieldWaterTankLV->GetMass(true, false) / kg << endl;
		cout << "     dimension : " << shieldWaterTankBox->GetXHalfLength() * 2 << " x "
				<< shieldWaterTankBox->GetYHalfLength() * 2 << " x "
				<< shieldWaterTankBox->GetZHalfLength() * 2 << endl;
		cout << "     thickness : " << watertank_thickness << "(top thickness: " << watertank_top_thickness << ")" << endl;
		cout << "     coordinate: " << shieldWaterTankPV->GetTranslation() + shieldWaterHousingPV->GetTranslation() + physCavern->GetTranslation() << endl;
		cout << " Hat Aluminium plate " << endl;
		cout << "     mass      : " << HatAlPlateLV->GetMass(true, false) / kg << endl;
		cout << "     dimension : " << HatBeamHousingInBox->GetXHalfLength() * 2 << " x "
				<< HatBeamHousingInBox->GetYHalfLength() * 2 << " x "
				<< HatBeamHousingInBox->GetZHalfLength() * 2 << endl;
		cout << "     thickness : " << HatAlPlate_thickness << endl;
		cout << "     coordinate: " << HatAlPlatePV->GetTranslation() + physCavern->GetTranslation() << endl;
		cout << " Hat boric acid " << endl;
		cout << "     mass      : " << shieldHatBoricLV->GetMass(true, false) / kg << endl;
		cout << "     dimension : " << HatAlPlateInBox->GetXHalfLength() * 2 << " x "
				<< HatAlPlateInBox->GetYHalfLength() * 2 << " x "
				<< HatAlPlateInBox->GetZHalfLength() * 2 << endl;
		cout << "     thickness : " << boricacid_thickness << endl;
		cout << "     coordinate: " << shieldHatBoricPV->GetTranslation() + physCavern->GetTranslation() << endl;
	}
	// }

	cout << " end " << endl;

	retvalLV = logiCavern;
	return retvalLV;
}

////////////////////////////////////////////////////////////////
/// H-beam shape is not H shape but square pillar.
static G4LogicalVolume *MakeHBeam1(G4Material *mat)
{
	G4double beam1Boxsize_x = 8800 / 2.;
	G4double beam1Boxsize_y = 4975 / 2.;
	G4double beam1Boxsize_z = 300;
	G4double ssize_x = 300 / 2.;
	G4double ssize_y = 675 / 2.;
	G4double ssize_z = 600;
	G4double lsize_x = 3600 / 2.;
	G4double lsize_y = 3375 / 2.;
	G4double lsize_z = 600;

	G4ThreeVector lposHole;
	G4ThreeVector rposHole;
	G4SubtractionSolid *HbeamBot1Solid[20];
	G4SubtractionSolid *HbeamBot2Solid[20];

	G4Box *Hbeam1Box = new G4Box("Hbeam1_Box", beam1Boxsize_x, beam1Boxsize_y, beam1Boxsize_z / 2.);
	G4Box *smallsquarehole = new G4Box("smallsquarehole", ssize_x, ssize_y, ssize_z / 2.);
	G4Box *largesquarehole = new G4Box("largesquarehole", lsize_x, lsize_y, lsize_z / 2.);

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			lposHole = G4ThreeVector(
				-beam1Boxsize_x + 400 + ssize_x + ssize_x * 2 * i + 200 * i,
				-beam1Boxsize_y + 400 + ssize_y + ssize_y * 2 * j + 200 * j,
				0);
			if (j == 0 && i == 0)
			{
				HbeamBot1Solid[j + i * 5] = new G4SubtractionSolid("HbeamBot_Solid",
																   Hbeam1Box, smallsquarehole, 0, lposHole);
			}
			else
			{
				HbeamBot1Solid[j + i * 5] = new G4SubtractionSolid("HbeamBot_Solid",
																   HbeamBot1Solid[j - 1 + i * 5], smallsquarehole, 0, lposHole);
			}
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			rposHole = G4ThreeVector(
				beam1Boxsize_x - 400 - ssize_x - ssize_x * 2 * i - 200 * i,
				-beam1Boxsize_y + 400 + ssize_y + ssize_y * 2 * j + 200 * j,
				0);
			if (j == 0 && i == 0)
			{
				HbeamBot2Solid[j + i * 5] = new G4SubtractionSolid("HbeamBot_Solid",
																   HbeamBot1Solid[19], smallsquarehole, 0, rposHole);
			}
			else
			{
				HbeamBot2Solid[j + i * 5] = new G4SubtractionSolid("HbeamBot_Solid",
																   HbeamBot2Solid[j - 1 + i * 5], smallsquarehole, 0, rposHole);
			}
		}
	}
	G4SubtractionSolid *HbeamBotSolid = new G4SubtractionSolid("HbeamBot_Solid",
															   HbeamBot2Solid[19], largesquarehole, 0, G4ThreeVector(0, 0, 0));
	G4LogicalVolume *HBeam1LV = new G4LogicalVolume(HbeamBotSolid, mat, "HbeamBot_LV");

	return HBeam1LV;
}

static G4LogicalVolume *MakePit(G4double pitBox_x, G4double pitBox_z, G4Material *mat, G4double HatBarrel_gap)
{
	// G4double pitBox_x = 9600/2.;
	// G4double pitBox_z = 6900;
	G4double pitBox_y = 3600 / 2.;
	G4double inBox1_y = (3600 - 600) / 2;
	G4double inBox1_x = 4200 / 2;
	G4double inBox2_x = 6000 / 2;
	G4double inBox3_y = 5000 / 2.;
	G4double beam1Boxsize_x = 8800 / 2.;
	G4double beam1Boxsize_y = 4975 / 2.;
	// G4double beam1Boxsize_z = 300;
	G4double beam1Boxsize_z = 460;

	G4Box *Hbeam1Box = new G4Box("Hbeam1_Box", beam1Boxsize_x, beam1Boxsize_y, beam1Boxsize_z / 2.);
	G4Box *pitMotherBox = new G4Box("pitMother_Box", pitBox_x, pitBox_y, pitBox_z / 2.);
	G4Box *pitInBox1 = new G4Box("pitInner1_Box", inBox1_x, inBox1_y, pitBox_z / 2.);
	G4Box *pitInBox2 = new G4Box("pitInner2_Box", inBox2_x, inBox1_y, pitBox_z / 2.);
	G4Box *pitInBox3 = new G4Box("pitInner2_Box", inBox2_x, inBox3_y, pitBox_z / 2.);

	G4VSolid *pit1Solid = new G4SubtractionSolid("pitSolid", pitMotherBox, Hbeam1Box, 0,
												 G4ThreeVector(pitBox_x - 2400, 0, pitBox_z / 2. - beam1Boxsize_z / 2.));
	G4VSolid *pit2Solid = new G4SubtractionSolid("pitSolid", pit1Solid, pitInBox1, 0,
												 G4ThreeVector(pitBox_x - inBox1_x - HatBarrel_gap, 0., HatBarrel_gap));
	G4VSolid *pit3Solid = new G4SubtractionSolid("pitSolid", pit2Solid, pitInBox2, 0,
												 G4ThreeVector(-pitBox_x + inBox2_x + HatBarrel_gap, 0, 2600 + HatBarrel_gap));
	G4VSolid *pitSolid = new G4SubtractionSolid("pitSolid", pit3Solid, pitInBox3, 0,
												G4ThreeVector(-inBox2_x, 0., -4300));

	G4LogicalVolume *pitLV = new G4LogicalVolume(pitSolid, mat, "pit_LV");
	return pitLV;
}
