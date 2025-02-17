#include "AmoreSim/AmoreParallelWorldConstruction.hh"
#include "AmoreSim/AmoreDetectorConstruction.hh"

#include "CupSim/CupPMTSD.hh"
#include "CupSim/CupParam.hh"
#include "CupSim/CupScintSD.hh"
#include "CupSim/Cup_PMT_LogicalVolume.hh"
#include "CupSim/CupPMTOpticalModel.hh"
#include "CupSim/CupDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4Region.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"

using namespace CLHEP;
using namespace std;

////////////////////////////////////////////////////////////////
// declaration of "private" static utility functions that we
// don't need in class definition
/*
static void MakeID_PMT_Support(Cup_PMT_LogicalVolume *that,
							   G4Material *SupportMat,
							   G4Material *ExteriorMat);
							   */

////////////////////////////////////////////////////////////////////
// CONSTRUCT PARALLAL WORLD for AMoRE200 WCMD's PMTs
//
void AmoreParallelWorldConstruction::ConstructAMoRE200_ParallelWorldForWCMD()
{
	// -- database
	CupParam &db(CupParam::GetDB());

	// Create optical surface for photocathode
	G4OpticalSurface *fPhotocathode_opsurf = new G4OpticalSurface("Photocathode_opsurf");
	fPhotocathode_opsurf->SetType(dielectric_metal); // ignored if RINDEX defined
	fPhotocathode_opsurf->SetMaterialPropertiesTable(
		G4Material::GetMaterial("photocathode")->GetMaterialPropertiesTable());

	// -- get pointer to physCavern, so we can put stuff in the pit
	G4VPhysicalVolume *fParaPhys = CupDetectorConstruction::GetPhysicalVolumeByName("paraPhysical");
	if (fParaPhys == NULL)
	{
		G4Exception(" ", " ", JustWarning, "Could not find fParaPhys!  Cannot build WCMD PMTs.");
	}
	else
		G4cout << "found fParaPhys!" << G4endl;

	// -- open pmt coordinates file, using AmoreDATA environment variable if set
	// G4std::ifstream whereWCPMT;
	ifstream whereWCPMT;
	const char *basic_fn = "pmtpos_amore200.dat";
	if (getenv("AmoreDATA") != NULL)
		whereWCPMT.open((G4String(getenv("AmoreDATA")) + "/" + G4String(basic_fn)).c_str());
	// print error message on failure of file open
	if (whereWCPMT.fail())
	{
		G4cerr << "Error, " << basic_fn << " could not be opened.\n";
		if (getenv("AmoreDATA") == NULL)
			G4cerr << "AmoreDATA environment variable is not set, so I was looking"
					  " for "
				   << basic_fn << " in the current directory." << G4endl;
		else
			G4cerr << "I was looking for it in the AmoreDATA directory, "
				   << getenv("AmoreDATA") << G4endl;
		G4Exception(" ", " ", JustWarning, "Error, pmt coordinates file could not be opened.\n");
	}
	// read max number of pmts
	// int maxWCPMTNo;
	whereWCPMT >> maxWCPMTNo;
	if (whereWCPMT.fail())
	{
		G4cerr << "Error, initial integer could not be read from pmt coordinates file.\n";
		G4Exception(" ", " ", JustWarning, "Error, initial integer could not be read from pmt coordinates file.\n");
	}

	// --- make the fundamental inner  PMT assembly
	Cup_PMT_LogicalVolume *_logiWCPMT10 = new Cup_10inch_LogicalVolume("WCPMT10",
																	   fWater,
																	   fGlass,
																	   fPhotocathode_opsurf,
																	   fPMT_Vac,
																	   fStainless, // dynode material
																				   /*( db["omit_pmt_masks"] != 0.0 ?
																					 NULL :       // no mask
																					 fBlackAcryl // physical mask on tubes to block non-sensitive areas
																					 )*/
																	   0,
																	   pmtSDWC, // sensitive detector hook
																	   // whichPmtStyle);
																	   kPmtStyle_TorusStack);

	Cup_PMT_LogicalVolume *_logiWCPMT8 = new Cup_8inch_LogicalVolume("WCPMT8",
																	   fWater,
																	   fGlass,
																	   fPhotocathode_opsurf,
																	   fPMT_Vac,
																	   fStainless, // dynode material
																				   /*( db["omit_pmt_masks"] != 0.0 ?
																					 NULL :       // no mask
																					 fBlackAcryl // physical mask on tubes to block non-sensitive areas
																					 )*/
																	   0,
																	   pmtSDWC, // sensitive detector hook
																	   // whichPmtStyle);
																	   kPmtStyle_TorusStack);
	/*
	MakeID_PMT_Support(_logiWCPMT10,
					   fWater,     // support material
					   fWater);   // external material
					   */

	// --- make the inner PMTs
	G4RotationMatrix *_rotWCPMT;

	G4int region_a, region_b, region_c;
	G4double cood_x, cood_y, cood_z, pmt_size;
	G4int WCPMTno;
	G4double angle_x, angle_z;
	char PMTname[64];
	G4ThreeVector new_x, new_y, new_z;

	G4double PMT_rotation_factor_for_cylindrical_ID =
		db.GetWithDefault("PMT_rotation_factor_for_cylindrical_ID", 0.0);

	// G4double angle_tank_corner = atan2(buffer_tank_height/2,buffer_tank_radius);

	// --- loop reading coordinates and making PMTs
	while (whereWCPMT.good())
	{
		// get a line from the file
		char linebuffer[128];
		whereWCPMT.getline(linebuffer, sizeof(linebuffer) - 1);
		if (whereWCPMT.fail())
			break;

		// skip blank lines and lines beginning with '#'
		if (linebuffer[0] == '#' || linebuffer[0] == '\0')
			continue;

		// put the line in an istrstream for convenient parsing
		// G4std::istrstream lineStream(linebuffer);
		istringstream lineStream(linebuffer);

		// parse out region, coordinates,
		region_a = region_b = region_c = -1;
		cood_x = cood_y = cood_z = 0.0;
		pmt_size = -911.;
		lineStream >> region_a >> region_b >> region_c >> cood_x >> cood_y >> cood_z >> pmt_size;

		// check for bad data
		if (lineStream.fail() || region_a < 0 || region_b < 0 || region_c < 0 || (cood_x == 0. && cood_y == 0. && cood_z == 0.) || pmt_size == -911.)
		{
			G4cerr << "BAD DATA in PMT file:  line=\"" << linebuffer << "\"\n";
			G4cerr.flush();
			continue;
		}

		// skip if pmt_size == 0
		if (pmt_size <= 0.0)
			continue;

		// a possible DAQ numbering
		WCPMTno = region_c;

		// name this PMT
		sprintf(PMTname, "physWCPMT%d", WCPMTno);

		// calculate angles and positions
		// we want PMTS to point normal to the surface they're mounted on
		// if rotation_factor==0.0, we want them to point to the center if
		// rotation_factor==1.0, and we want something in between for in-between
		// values.
		G4double r = sqrt(cood_x * cood_x + cood_y * cood_y + cood_z * cood_z);
		G4double dx = -cood_x / r;
		G4double dy = -cood_y / r;
		G4double dz = -cood_z / r;

		angle_z = atan2(dx, dy);
		angle_x = atan2(dz, sqrt(dx * dx + dy * dy));

		PMT_rotation_factor_for_cylindrical_ID = 0.0; 

		double normal_angle = -M_PI / 2; // only top PMTs
		angle_x = normal_angle + (angle_x - normal_angle) * PMT_rotation_factor_for_cylindrical_ID;

		_rotWCPMT = new G4RotationMatrix();
		_rotWCPMT->rotateZ(angle_z);
		_rotWCPMT->rotateX(M_PI / 2.0 - angle_x);

		cood_x = cood_x;
		cood_y = cood_y;
		cood_z = cood_z;

		G4ThreeVector pmtpos(cood_x, cood_y, cood_z);

		// creating the new PVPlacement automatically adds it to the
		//  singleton G4PhysicalVolumeStore()
		// ****************************************************************
		// * Use the constructor that specifies the PHYSICAL mother, since
		// * each PMT occurs only once in one physical volume.  This saves
		// * the GeometryManager some work. -GHS.
		// ****************************************************************
		if (db["omit_id_pmts"] == 0.0)
		{
			if (pmt_size == 10.0)
			{
				new G4PVPlacement(_rotWCPMT,
								  pmtpos,
								  PMTname,
								  _logiWCPMT10,
								  // BufferInteriorPhys,    // physical parent
								  fParaPhys,
								  false,
								  WCPMTno);
			}
			else if (pmt_size == 8.0)
			{
				new G4PVPlacement(_rotWCPMT,
								  pmtpos,
								  PMTname,
								  _logiWCPMT8,
								  // BufferInteriorPhys,    // physical parent
								  fParaPhys,
								  false,
								  WCPMTno);
			}
		}
	}
	whereWCPMT.close();
	// (all done setting inner PMTs from file)
}

////////////////////////////////////////////////////////////////
// definition of "private" static utility functions that we
// don't need in class definition
/*
static void MakeID_PMT_Support(Cup_PMT_LogicalVolume *that,
							   G4Material *SupportMat,
							   G4Material *ExteriorMat)
{
  // ... this should build the PMT support geometry
  // ... empty function for now
}
*/
