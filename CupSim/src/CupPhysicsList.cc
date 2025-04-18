//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo
//                    (27th November 2001)
//
// PhysicsList program
//
// Modified:
//
// 14-02-03 Fix bugs in msc and hIon instanciation + cut per region
//
// 05-02-05 AH - changes to G4Decay - added is not short lived protection
//          and redefined particles to allow non-static creation
//          i.e. changed construction to G4MesonConstructor, G4BaryonConstructor
//
// 23-10-09 LP - migrated EM physics from the LowEnergy processes (not supported) to
//          the new G4Livermore model implementation. Results unchanged.
//
// --------------------------------------------------------------

#include "G4Version.hh"

//#if G4VERSION_NUMBER >= 1000
#include <iomanip>

#include "CupSim/CupParam.hh"
#include "CupSim/CupPhysicsList.hh"
#include "CupSim/CupPhysicsListMessenger.hh"
#include "CupSim/CupDeferTrackProc.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4ios.hh"
#include "G4UserLimits.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4StepLimiter.hh"

// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4BetheHeitler5DModel.hh"

#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"

// e-
#include "G4eMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4UniversalFluctuation.hh"

// e+
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

// alpha and GenericIon and deuterons, triton, He3:
//muon:
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuMultipleScattering.hh"
#include "G4MuonMinusCapture.hh"

//OTHERS:
#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"
#include "G4hBremsstrahlung.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"

//em process options to allow msc step-limitation to be switched off
#include "G4EmParameters.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"

#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
//#include "G4OpBoundaryProcess.hh" // EJ: replaced by CupOpBoundaryProcess.hh
#include "G4OpticalParameters.hh"

// Elastic processes:
#include "G4HadronElasticProcess.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"

// Inelastic processes:
#include "G4HadronInelasticProcess.hh"

// High energy FTFP model and Bertini cascade
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4CascadeInterface.hh"

// Cross sections
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CrossSectionElastic.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4AntiNuclElastic.hh"

#include "G4CrossSectionInelastic.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGNucleonElasticXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronElasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"

#include "G4HadronElastic.hh"
#include "G4NeutronCaptureProcess.hh"

// Neutron high-precision models: <20 MeV
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPCapture.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPInelasticData.hh"

// Stopping processes
#include "G4HadronStoppingProcess.hh"
#include "G4HadronicAbsorptionBertini.hh"
#include "G4HadronicAbsorptionFritiof.hh"

#include "G4HadronicParameters.hh"

#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4PhysicsListHelper.hh"
#include "G4NuclideTable.hh"
#include "G4NuclearLevelData.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4EmStandardPhysics.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4ios.hh"

G4ProductionCuts *CupPhysicsList::DetectorCuts = nullptr;
G4double CupPhysicsList::cutForGamma           = 0;
G4double CupPhysicsList::cutForElectron        = 0;
G4double CupPhysicsList::cutForPositron        = 0;
G4double CupPhysicsList::cutForProton          = 0;
G4double CupPhysicsList::cutForAlpha           = 0;
G4double CupPhysicsList::cutForGenericIon      = 0;
G4bool CupPhysicsList::enableCrystalRegion     = 0;
G4bool CupPhysicsList::omitHadronicProc        = 0;
G4bool CupPhysicsList::omitNeutHP              = 0;
G4int CupPhysicsList::VerboseLevel             = 0;
G4int CupPhysicsList::OpVerbLevel              = 0;

// Constructor /////////////////////////////////////////////////////////////
CupPhysicsList::CupPhysicsList() : G4VModularPhysicsList() {
	defaultCutValue = 1. * mm; //
	cutForGamma     = defaultCutValue;
	cutForElectron  = defaultCutValue;
	cutForPositron  = defaultCutValue;

	DetectorCuts = 0;

	VerboseLevel = 1;
	OpVerbLevel  = 0;

	SetVerboseLevel(VerboseLevel);

	// EJ Messenger
	pMessenger       = new CupPhysicsListMessenger(this);
	CupParam &db     = CupParam::GetDB();
	omitHadronicProc = (db["omit_hadronic_processes"] != 0.0);
	//G4cout << "EJ: omitHadronicProc = " << omitHadronicProc << G4endl;
	if (omitHadronicProc) {
		G4cerr << "Warning, Hadronic processes omitted.\n";
		//return;
	}
	omitNeutHP = (db["omit_neutron_hp"] != 0.0);
	if (omitNeutHP) {
		G4cerr << "WARNING: --- OMITTING neutron_hp model! ---" << G4endl;
	} else {
		G4cerr << "Note: +++ INCLUDING neutron_hp model. +++" << G4endl;
	}

	// Hadronic physics extra configuration
  	G4double hlThreshold = 1000 * picosecond;
    	G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(hlThreshold);

  	//G4EmParameters::Instance()->AddPhysics("World","G4RadioactiveDecay");
  	G4EmParameters::Instance()->AddPhysics("World","G4Radioactivation");
  	G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
  	deex->SetStoreICLevelData(true);
  	deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()
                       /std::log(2.));
	// radioactivation of materials
  	deex->SetIsomerProduction(true);
  	deex->SetCorrelatedGamma(true);
  	SetVerboseLevel(VerboseLevel);
}

// Destructor //////////////////////////////////////////////////////////////
CupPhysicsList::~CupPhysicsList() { delete pMessenger; }

// Construct Particles /////////////////////////////////////////////////////
void CupPhysicsList::ConstructParticle() {

	// In this method, static member functions should be called
	// for all particles which yodu want to use.
	// This ensures that objects of these particle types will be
	// created in the program.

	ConstructMyBosons();
	ConstructMyLeptons();
	ConstructMyHadrons();
	ConstructMyShortLiveds();
}

// construct Bosons://///////////////////////////////////////////////////
void CupPhysicsList::ConstructMyBosons() {
	// pseudo-particles
	G4Geantino::GeantinoDefinition();
	G4ChargedGeantino::ChargedGeantinoDefinition();

	// gamma
	G4Gamma::GammaDefinition();

	// OpticalPhotons
	G4OpticalPhoton::OpticalPhotonDefinition();
}

// construct Leptons://///////////////////////////////////////////////////
void CupPhysicsList::ConstructMyLeptons() {
	// leptons
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();
	G4MuonPlus::MuonPlusDefinition();
	G4MuonMinus::MuonMinusDefinition();

	G4NeutrinoE::NeutrinoEDefinition();
	G4AntiNeutrinoE::AntiNeutrinoEDefinition();
	G4NeutrinoMu::NeutrinoMuDefinition();
	G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4MesonConstructor.hh"

// construct Hadrons://///////////////////////////////////////////////////
void CupPhysicsList::ConstructMyHadrons() {
	//  mesons
	G4MesonConstructor mConstructor;
	mConstructor.ConstructParticle();

	//  baryons
	G4BaryonConstructor bConstructor;
	bConstructor.ConstructParticle();

	//  ions
	G4IonConstructor iConstructor;
	iConstructor.ConstructParticle();
}

#include "G4ShortLivedConstructor.hh"
// construct Shortliveds://///////////////////////////////////////////////////
void CupPhysicsList::ConstructMyShortLiveds() {
	// ShortLiveds
	G4ShortLivedConstructor slConstructor;
	slConstructor.ConstructParticle();
}

// Construct Processes //////////////////////////////////////////////////////
void CupPhysicsList::ConstructProcess() {

	AddTransportation();

	AddParameterisation();
	// ConstructEM();
	if (emName == "livermore") {
		auto a = new G4EmLivermorePhysics;
		a->ConstructProcess();
	} else {
		ConstructEM();
		G4cout << "EM Physics is ConstructEM(): default!" << G4endl;
	}
	// ConstructOp();
	ConstructOp();

	ConstructHad();

	ConstructGeneral();
}

// G4FastSimulation Processes //////////////////////////////////////////////////////
#include "G4FastSimulationManagerProcess.hh"

void CupPhysicsList::AddParameterisation() {
	auto theParticleIterator = GetParticleIterator();

	G4FastSimulationManagerProcess *theFastSimulationManagerProcess =
		new G4FastSimulationManagerProcess();
	theParticleIterator->reset();
	while ((*theParticleIterator)()) {
		G4ParticleDefinition *particle = theParticleIterator->value();
		G4ProcessManager *pmanager     = particle->GetProcessManager();
		// both postStep and alongStep action are required if the detector
		// makes use of ghost volumes. If no ghost, the postStep
		// is sufficient (and faster?).
#define Cup_USES_GHOST_VOLUMES 0
#if Cup_USES_GHOST_VOLUMES
		pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);
#else
		pmanager->AddProcess(theFastSimulationManagerProcess, -1, -1, 1);
#endif
	}
}

// Electromagnetic Processes ////////////////////////////////////////////////
// all charged particles

// gamma
#include "G4LivermorePhotoElectricModel.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"

#include "G4LivermoreRayleighModel.hh"
#include "G4RayleighScattering.hh"

// e-
#include "G4eMultipleScattering.hh"

#include "G4LivermoreIonisationModel.hh"
#include "G4eIonisation.hh"

#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4eBremsstrahlung.hh"

// e+
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eplusAnnihilation.hh"

// alpha and GenericIon and deuterons, triton, He3:
#include "G4EnergyLossTables.hh"

// muon:
#include "G4MuBremsstrahlung.hh"
#include "G4MuIonisation.hh"
#include "G4MuPairProduction.hh"
//#include "G4MuonMinusCaptureAtRest.hh"
#include "G4MuonMinusCapture.hh"

// OTHERS:
#include "G4IonParametrisedLossModel.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"
#include "G4ionIonisation.hh"

#include "G4LossTableManager.hh"

void CupPhysicsList::ConstructEM() {

	// set a finer grid of the physic tables in order to improve precision
	// former LowEnergy models have 200 bins up to 100 GeV

	auto theParticleIterator = GetParticleIterator();

	theParticleIterator->reset();
	while ((*theParticleIterator)()) {
		G4ParticleDefinition *particle = theParticleIterator->value();
		G4ProcessManager *pmanager     = particle->GetProcessManager();
		G4String particleName          = particle->GetParticleName();
		G4String particleType          = particle->GetParticleType();
		G4double charge                = particle->GetPDGCharge();

		if (particleName == "gamma") {
			// gamma
			G4RayleighScattering *theRayleigh = new G4RayleighScattering();
			theRayleigh->SetEmModel(new G4LivermoreRayleighModel()); // not strictly necessary
			pmanager->AddDiscreteProcess(theRayleigh);

			G4PhotoElectricEffect *thePhotoElectricEffect = new G4PhotoElectricEffect();
			thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
			pmanager->AddDiscreteProcess(thePhotoElectricEffect);

			G4ComptonScattering *theComptonScattering = new G4ComptonScattering();
			theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
			pmanager->AddDiscreteProcess(theComptonScattering);

			G4GammaConversion *theGammaConversion = new G4GammaConversion();
			theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel());
			pmanager->AddDiscreteProcess(theGammaConversion);

		} else if (particleName == "e-") {
			// electron
			// process ordering: AddProcess(name, at rest, along step, post step)
			// Multiple scattering
			G4eMultipleScattering *msc = new G4eMultipleScattering();
			pmanager->AddProcess(msc, -1, 1, 1);

			// Ionisation
			G4eIonisation *eIonisation = new G4eIonisation();
			eIonisation->SetEmModel(new G4LivermoreIonisationModel());
			eIonisation->SetStepFunction(0.2, 100 * um); // improved precision in tracking
			pmanager->AddProcess(eIonisation, -1, 2, 2);

			// Bremsstrahlung
			G4eBremsstrahlung *eBremsstrahlung = new G4eBremsstrahlung();
			eBremsstrahlung->SetEmModel(new G4LivermoreBremsstrahlungModel());
			pmanager->AddProcess(eBremsstrahlung, -1, -3, 3);
		} else if (particleName == "e+") {
			// positron
			G4eMultipleScattering *msc = new G4eMultipleScattering();
			pmanager->AddProcess(msc, -1, 1, 1);

			// Ionisation
			G4eIonisation *eIonisation = new G4eIonisation();
			eIonisation->SetStepFunction(0.2, 100 * um); //
			pmanager->AddProcess(eIonisation, -1, 2, 2);

			// Bremsstrahlung (use default, no low-energy available)
			pmanager->AddProcess(new G4eBremsstrahlung(), -1, -1, 3);

			// Annihilation
			pmanager->AddProcess(new G4eplusAnnihilation(), 0, -1, 4);
		} else if (particleName == "mu+" || particleName == "mu-") {
			// muon
			pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
			pmanager->AddProcess(new G4MuIonisation(), -1, 2, 2);
			pmanager->AddProcess(new G4MuBremsstrahlung(), -1, -1, 3);
			pmanager->AddProcess(new G4MuPairProduction(), -1, -1, 4);
			if (particleName == "mu-")
				// pmanager->AddProcess(new G4MuonMinusCaptureAtRest(), 0,-1,-1);
				pmanager->AddProcess(new G4MuonMinusCapture(), 0, -1, -1);
		} else if (particleName == "proton" || particleName == "pi+" || particleName == "pi-") {
			// multiple scattering
			pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

			// ionisation
			G4hIonisation *hIonisation = new G4hIonisation();
			hIonisation->SetStepFunction(0.2, 50 * um);
			pmanager->AddProcess(hIonisation, -1, 2, 2);

			// bremmstrahlung
			pmanager->AddProcess(new G4hBremsstrahlung, -1, -3, 3);
		} else if (particleName == "alpha" || particleName == "deuteron" ||
				particleName == "triton" || particleName == "He3") {
			// multiple scattering
			pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

			// ionisation
			G4ionIonisation *ionIoni = new G4ionIonisation();
			ionIoni->SetStepFunction(0.1, 20 * um);
			pmanager->AddProcess(ionIoni, -1, 2, 2);
		} else if (particleName == "GenericIon") {
			// OBJECT may be dynamically created as either a GenericIon or nucleus
			// G4Nucleus exists and therefore has particle type nucleus
			// genericIon:

			// multiple scattering
			pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

			// ionisation
			G4ionIonisation *ionIoni = new G4ionIonisation();
			ionIoni->SetEmModel(new G4IonParametrisedLossModel());
			ionIoni->SetStepFunction(0.1, 20 * um);
			pmanager->AddProcess(ionIoni, -1, 2, 2);
		}

		else if ((!particle->IsShortLived()) && (charge != 0.0) &&
				(particle->GetParticleName() != "chargedgeantino")) {
			// all others charged particles except geantino
			G4hMultipleScattering *aMultipleScattering = new G4hMultipleScattering();
			G4hIonisation *ahadronIon                  = new G4hIonisation();

			// multiple scattering
			pmanager->AddProcess(aMultipleScattering, -1, 1, 1);

			// ionisation
			pmanager->AddProcess(ahadronIon, -1, 2, 2);
		}
	}

}

// Optical Processes ////////////////////////////////////////////////////////
// EJ: start
#include "CupSim/CupOpAttenuation.hh"
#include "CupSim/CupOpBoundaryProcess.hh"
#include "CupSim/CupScintillation.hh"
#include "G4Cerenkov.hh"
#include "G4EmSaturation.hh"
// EJ: end

void CupPhysicsList::ConstructOp() {
	// EJ: start
	// scintillation process
	CupScintillation *theScintProcessDef = new CupScintillation("Scintillation");
	// theScintProcessDef->DumpPhysicsTable();
	theScintProcessDef->SetTrackSecondariesFirst(true);
	theScintProcessDef->SetScintillationYieldFactor(1.0);     //
	theScintProcessDef->SetScintillationExcitationRatio(0.0); //
	theScintProcessDef->SetVerboseLevel(OpVerbLevel);

	G4EmSaturation *emSaturation = G4LossTableManager::Instance()->EmSaturation();
	theScintProcessDef->AddSaturation(emSaturation);

	// scintillation process for alpha:
	CupScintillation *theScintProcessAlpha = new CupScintillation("Scintillation");
	// theScintProcessNuc->DumpPhysicsTable();
	theScintProcessAlpha->SetTrackSecondariesFirst(true);
	theScintProcessAlpha->SetScintillationYieldFactor(1.1);
	theScintProcessAlpha->SetScintillationExcitationRatio(1.0);
	theScintProcessAlpha->SetVerboseLevel(OpVerbLevel);

	theScintProcessAlpha->AddSaturation(emSaturation);

	// scintillation process for heavy nuclei
	CupScintillation *theScintProcessNuc = new CupScintillation("Scintillation");
	// theScintProcessNuc->DumpPhysicsTable();
	theScintProcessNuc->SetTrackSecondariesFirst(true);
	theScintProcessNuc->SetScintillationYieldFactor(0.2);
	theScintProcessNuc->SetScintillationExcitationRatio(1.0);
	theScintProcessNuc->SetVerboseLevel(OpVerbLevel);

	theScintProcessNuc->AddSaturation(emSaturation);

	// optical processes
	CupOpAttenuation *theAttenuationProcess = new CupOpAttenuation();
	theAttenuationProcess->UseTimeProfile("exponential");
	theAttenuationProcess->SetVerboseLevel(OpVerbLevel);

	//G4OpBoundaryProcess *theBoundaryProcess = new G4OpBoundaryProcess();
	CupOpBoundaryProcess *theBoundaryProcess = new CupOpBoundaryProcess();
	theBoundaryProcess->SetVerboseLevel(OpVerbLevel);

	// Cerenkov
	G4Cerenkov *theCerenkovProcess = new G4Cerenkov();
	theCerenkovProcess->SetTrackSecondariesFirst(true);

	auto theParticleIterator = GetParticleIterator();

	theParticleIterator->reset();
	while ((*theParticleIterator)()) {
		G4ParticleDefinition *particle = theParticleIterator->value();
		G4ProcessManager *pmanager     = particle->GetProcessManager();
		G4String particleName          = particle->GetParticleName();
		if (theScintProcessDef->IsApplicable(*particle)) {
			//      if(particle->GetPDGMass() > 5.0*GeV)
			if (particle->GetParticleName() == "GenericIon") {
				pmanager->AddProcess(theScintProcessNuc); // AtRestDiscrete
				pmanager->SetProcessOrderingToLast(theScintProcessNuc, idxAtRest);
				pmanager->SetProcessOrderingToLast(theScintProcessNuc, idxPostStep);
			} else if (particle->GetParticleName() == "alpha") {
				pmanager->AddProcess(theScintProcessAlpha);
				pmanager->SetProcessOrderingToLast(theScintProcessAlpha, idxAtRest);
				pmanager->SetProcessOrderingToLast(theScintProcessAlpha, idxPostStep);
			} else {
				pmanager->AddProcess(theScintProcessDef);
				pmanager->SetProcessOrderingToLast(theScintProcessDef, idxAtRest);
				pmanager->SetProcessOrderingToLast(theScintProcessDef, idxPostStep);
				pmanager->AddProcess(theCerenkovProcess);
				pmanager->SetProcessOrdering(theCerenkovProcess, idxPostStep);
			}
		}

		if (particleName == "opticalphoton") {
			pmanager->AddDiscreteProcess(theAttenuationProcess);
			pmanager->AddDiscreteProcess(theBoundaryProcess);
		}
	}
}

// Hadronic processes ////////////////////////////////////////////////////////

// Elastic processes:
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "G4HadronElasticProcess.hh"

// Inelastic processes:
#include "G4HadronInelasticProcess.hh"

// High energy FTFP model and Bertini cascade
#include "G4CascadeInterface.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4FTFModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4LundStringFragmentation.hh"
#include "G4PreCompoundModel.hh"
#include "G4TheoFSGenerator.hh"

// Cross sections
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CrossSectionElastic.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4AntiNuclElastic.hh"

#include "G4CrossSectionInelastic.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGNucleonElasticXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronElasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"

#include "G4HadronElastic.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4NeutronFissionProcess.hh"

// Neutron high-precision models: <20 MeV
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPCapture.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPFissionData.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPFission.hh"

// Stopping processes
#include "G4HadronStoppingProcess.hh"
#include "G4HadronicAbsorptionBertini.hh"
#include "G4HadronicAbsorptionFritiof.hh"

#include "G4HadronicParameters.hh"

#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4PhysicsListHelper.hh"
#include "G4NuclideTable.hh"
#include "G4NuclearLevelData.hh"

// EJ
// Muon Nuclear processes
#include "G4MuonNuclearProcess.hh"
#include "G4MuonVDNuclearModel.hh"

void CupPhysicsList::ConstructHad() {
	if(omitHadronicProc) return;

  //Elastic models
  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE();

  // Inelastic scattering
  const G4double theFTFMin0 =    0.0*GeV;
  const G4double theFTFMin1 =    3.0*GeV;
  const G4double theFTFMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  const G4double theBERTMin0 =   0.0*GeV;
  const G4double theBERTMin1 =  19.0*MeV;
  const G4double theBERTMax =    6.0*GeV;
  const G4double theHPMin =      0.0*GeV;
  const G4double theHPMax =     20.0*MeV;

	G4FTFModel *theStringModel           = new G4FTFModel;
	G4ExcitedStringDecay *theStringDecay = new G4ExcitedStringDecay(new G4LundStringFragmentation);
	theStringModel->SetFragmentationModel(theStringDecay);
	G4PreCompoundModel *thePreEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
	G4GeneratorPrecompoundInterface *theCascade =
		new G4GeneratorPrecompoundInterface(thePreEquilib);

	G4TheoFSGenerator *theFTFModel0 = new G4TheoFSGenerator("FTFP");
	theFTFModel0->SetHighEnergyGenerator(theStringModel);
	theFTFModel0->SetTransport(theCascade);
	theFTFModel0->SetMinEnergy(theFTFMin0);
	theFTFModel0->SetMaxEnergy(theFTFMax);

	G4TheoFSGenerator *theFTFModel1 = new G4TheoFSGenerator("FTFP");
	theFTFModel1->SetHighEnergyGenerator(theStringModel);
	theFTFModel1->SetTransport(theCascade);
	theFTFModel1->SetMinEnergy(theFTFMin1);
	theFTFModel1->SetMaxEnergy(theFTFMax);

	G4CascadeInterface *theBERTModel0 = new G4CascadeInterface;
	theBERTModel0->SetMinEnergy(theBERTMin0);
	theBERTModel0->SetMaxEnergy(theBERTMax);

	G4CascadeInterface *theBERTModel1 = new G4CascadeInterface;
	theBERTModel1->SetMinEnergy(theBERTMin1);
	theBERTModel1->SetMaxEnergy(theBERTMax);

  G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  G4ComponentGGNuclNuclXsc * ggNuclNuclXsec = new G4ComponentGGNuclNuclXsc();
  G4VCrossSectionDataSet * theGGNuclNuclData = new G4CrossSectionInelastic(ggNuclNuclXsec);
  G4VCrossSectionDataSet * theGGNNEl = new G4CrossSectionElastic(ggNuclNuclXsec);
  G4ComponentGGHadronNucleusXsc * ggHNXsec = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet * theGGHNEl = new G4CrossSectionElastic(ggHNXsec);
  G4VCrossSectionDataSet * theGGHNInel = new G4CrossSectionInelastic(ggHNXsec);

	auto theParticleIterator = GetParticleIterator();

	theParticleIterator->reset();
	while ((*theParticleIterator)()) {
		G4ParticleDefinition *particle = theParticleIterator->value();
		G4ProcessManager *pmanager     = particle->GetProcessManager();
		G4String particleName          = particle->GetParticleName();

		if (particleName == "pi+") {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_he );
          pmanager->AddDiscreteProcess( theElasticProcess );
          //Inelastic scattering
          G4HadronInelasticProcess* theInelasticProcess = 
            new G4HadronInelasticProcess( "inelastic", G4PionPlus::Definition() );
          theInelasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 ); 
		  pmanager->AddDiscreteProcess( theInelasticProcess );
		}
		else if (particleName == "pi-") {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_he );
          pmanager->AddDiscreteProcess( theElasticProcess );
          //Inelastic scattering
          G4HadronInelasticProcess* theInelasticProcess =
            new G4HadronInelasticProcess( "inelastic", G4PionMinus::Definition() );
          theInelasticProcess->AddDataSet( new G4BGGPionInelasticXS( particle ) );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
          //Absorption
          pmanager->AddRestProcess(new G4HadronicAbsorptionBertini(G4PionMinus::Definition()), ordDefault);
		}
		else if (particleName == "kaon+") {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( theGGHNEl );
          theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering       
          G4HadronInelasticProcess* theInelasticProcess =
            new G4HadronInelasticProcess( "inelastic", G4KaonPlus::Definition() );
          theInelasticProcess->AddDataSet( theGGHNInel );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
		}
		else if (particleName == "kaon0S") {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( theGGHNEl );
          theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering        
          G4HadronInelasticProcess* theInelasticProcess =
            new G4HadronInelasticProcess( "inelastic", G4KaonZeroShort::Definition() );
          theInelasticProcess->AddDataSet( theGGHNInel );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
		}
		else if (particleName == "kaon0L") {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( theGGHNEl );
          theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4HadronInelasticProcess* theInelasticProcess =
            new G4HadronInelasticProcess( "inelastic", G4KaonZeroLong::Definition() );
          theInelasticProcess->AddDataSet( theGGHNInel );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
		}

		else if (particleName == "kaon-") {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( theGGHNEl );
          theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4HadronInelasticProcess* theInelasticProcess =
            new G4HadronInelasticProcess( "inelastic", G4KaonMinus::Definition() );
          theInelasticProcess->AddDataSet( theGGHNInel );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
          pmanager->AddRestProcess(new G4HadronicAbsorptionBertini(G4KaonMinus::Definition()), ordDefault);
		}

		else if (particleName == "proton") {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGNucleonElasticXS( G4Proton::Proton() ) );
          theElasticProcess->RegisterMe( elastic_chip );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4HadronInelasticProcess* theInelasticProcess =
            new G4HadronInelasticProcess( "inelastic", G4Proton::Definition() );
          theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
		}

		else if (particleName == "anti_proton") {
          // Elastic scattering
          const G4double elastic_elimitAntiNuc = 100.0*MeV;
          G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
          elastic_anuc->SetMinEnergy( elastic_elimitAntiNuc );
          G4CrossSectionElastic* elastic_anucxs = new G4CrossSectionElastic( elastic_anuc->GetComponentCrossSection() );    
          G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
          elastic_lhep2->SetMaxEnergy( elastic_elimitAntiNuc );
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( elastic_anucxs );
          theElasticProcess->RegisterMe( elastic_lhep2 );
          theElasticProcess->RegisterMe( elastic_anuc );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4HadronInelasticProcess* theInelasticProcess = 
            new G4HadronInelasticProcess( "inelastic", G4AntiProton::Definition() );
          theInelasticProcess->AddDataSet( theAntiNucleonData );
          theInelasticProcess->RegisterMe( theFTFModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
          // Absorption
          pmanager->AddRestProcess(new G4HadronicAbsorptionFritiof(G4AntiProton::Definition()), ordDefault);
		}

		else if (particleName == "neutron") {
			// elastic scattering
			G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
        		theElasticProcess->AddDataSet(new G4NeutronElasticXS());
        	G4HadronElastic* elastic_neutronChipsModel = new G4ChipsElasticModel();
			if (omitNeutHP) {
				elastic_neutronChipsModel->SetMinEnergy( 19.0*MeV );
        		theElasticProcess->RegisterMe( elastic_neutronChipsModel );
			} else {
				elastic_neutronChipsModel->SetMinEnergy( 19.0*MeV );
        		theElasticProcess->RegisterMe( elastic_neutronChipsModel );
				G4ParticleHPElastic * theElasticNeutronHP = new G4ParticleHPElastic;
        		theElasticNeutronHP->SetMinEnergy( theHPMin );
        		theElasticNeutronHP->SetMaxEnergy( theHPMax );
				theElasticProcess->RegisterMe( theElasticNeutronHP );
				theElasticProcess->AddDataSet( new G4ParticleHPElasticData );
			}
			pmanager->AddDiscreteProcess( theElasticProcess );
			// inelastic scattering
			G4HadronInelasticProcess* theInelasticProcess =
	  		  new G4HadronInelasticProcess( "inelastic", G4Neutron::Definition() );
				theInelasticProcess->AddDataSet( new G4NeutronInelasticXS() );
			if (omitNeutHP) {
				theInelasticProcess->RegisterMe(theFTFModel1);
				theInelasticProcess->RegisterMe(theBERTModel1);
			} else {
				theInelasticProcess->RegisterMe(theFTFModel1);
				theInelasticProcess->RegisterMe(theBERTModel1);
				G4ParticleHPInelastic * theNeutronInelasticHPModel = new G4ParticleHPInelastic;
				theNeutronInelasticHPModel->SetMinEnergy(theHPMin);
				theNeutronInelasticHPModel->SetMaxEnergy(theHPMax);
				theInelasticProcess->RegisterMe(theNeutronInelasticHPModel);
				theInelasticProcess->AddDataSet(new G4ParticleHPInelasticData);
			}
			pmanager->AddDiscreteProcess(theInelasticProcess);
			// capture
			G4NeutronCaptureProcess* theCaptureProcess =
	  		  new G4NeutronCaptureProcess;
			G4ParticleHPCapture * theLENeutronCaptureModel = new G4ParticleHPCapture;
			if (omitNeutHP) {
				theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
			} else {
				theLENeutronCaptureModel->SetMinEnergy(theHPMin);
				theLENeutronCaptureModel->SetMaxEnergy(theHPMax);
				theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
				theCaptureProcess->AddDataSet(new G4ParticleHPCaptureData);
			}
			pmanager->AddDiscreteProcess(theCaptureProcess);
			// fission
			G4NeutronFissionProcess *theNeutronFissionProcess = new G4NeutronFissionProcess();
			G4ParticleHPFission *theNeutronLFission          = new G4ParticleHPFission();
			if (omitNeutHP) {
				theNeutronFissionProcess->RegisterMe(theNeutronLFission);
			} else {
				G4ParticleHPFission *theNeutronHPFission = new G4ParticleHPFission();
				theNeutronHPFission->SetMaxEnergy(20. * MeV);
				theNeutronHPFission->SetMinEnergy(20. * MeV);
				theNeutronFissionProcess->RegisterMe(theNeutronHPFission);
				theNeutronFissionProcess->RegisterMe(theNeutronLFission);
				theNeutronFissionProcess->AddDataSet(new G4ParticleHPFissionData);
			}
			pmanager->AddDiscreteProcess(theNeutronFissionProcess);
		
		} else if (particleName == "anti_neutron") {
			// Elastic scattering
			G4HadronElasticProcess *theElasticProcess = new G4HadronElasticProcess;
			theElasticProcess->AddDataSet(theAntiNucleonData);	
			theElasticProcess->RegisterMe(elastic_lhep0);
			pmanager->AddDiscreteProcess(theElasticProcess);
			// Inelastic scattering (include annihilation on-fly)
			G4HadronInelasticProcess* theInelasticProcess = 
	    	  new G4HadronInelasticProcess( "inelastic", G4AntiNeutron::Definition() );
			theInelasticProcess->AddDataSet(theAntiNucleonData);
			theInelasticProcess->RegisterMe(theFTFModel0);
			pmanager->AddDiscreteProcess(theInelasticProcess);
		}

		else if (particleName == "deuteron") {
			// Elastic scattering
			G4HadronElasticProcess *theElasticProcess = new G4HadronElasticProcess;
			theElasticProcess->AddDataSet(theGGNuclNuclData);
			theElasticProcess->RegisterMe(elastic_lhep0);
			pmanager->AddDiscreteProcess(theElasticProcess);
			// Inelastic scattering
			G4HadronInelasticProcess* theInelasticProcess = 
	    	  new G4HadronInelasticProcess( "inelastic", G4Deuteron::Definition() );
			theInelasticProcess->AddDataSet(theGGNuclNuclData);
			theInelasticProcess->RegisterMe(theFTFModel1);
			theInelasticProcess->RegisterMe(theBERTModel0);
			pmanager->AddDiscreteProcess(theInelasticProcess);
		}

		else if (particleName == "triton") {
			// Elastic scattering
			G4HadronElasticProcess *theElasticProcess = new G4HadronElasticProcess;
			theElasticProcess->AddDataSet(theGGNuclNuclData);
			theElasticProcess->RegisterMe(elastic_lhep0);
			pmanager->AddDiscreteProcess(theElasticProcess);
			// Inelastic scattering
			G4HadronInelasticProcess* theInelasticProcess = 
	    	  new G4HadronInelasticProcess( "inelastic", G4Triton::Definition() );
			theInelasticProcess->AddDataSet(theGGNuclNuclData);
			theInelasticProcess->RegisterMe(theFTFModel1);
			theInelasticProcess->RegisterMe(theBERTModel0);
			pmanager->AddDiscreteProcess(theInelasticProcess);
		}

		else if (particleName == "alpha") {
			// Elastic scattering
			G4HadronElasticProcess *theElasticProcess = new G4HadronElasticProcess;
			theElasticProcess->AddDataSet(theGGNuclNuclData);
			theElasticProcess->RegisterMe(elastic_lhep0);
			pmanager->AddDiscreteProcess(theElasticProcess);
			// Inelastic scattering
			G4HadronInelasticProcess* theInelasticProcess = 
	    	  new G4HadronInelasticProcess( "inelastic", G4Alpha::Definition() );	
			theInelasticProcess->AddDataSet(theGGNuclNuclData);
			theInelasticProcess->RegisterMe(theFTFModel1);
			theInelasticProcess->RegisterMe(theBERTModel0);
			pmanager->AddDiscreteProcess(theInelasticProcess);
		}

		// EJ: for muon nuclear processes
		else if (particleName == "mu-") {
			G4MuonNuclearProcess *muNuclearProcess = new G4MuonNuclearProcess();
			G4MuonVDNuclearModel *muNuclearModel   = new G4MuonVDNuclearModel();
			muNuclearProcess->RegisterMe(muNuclearModel);
			pmanager->AddDiscreteProcess(muNuclearProcess);
		}
	}
}

// Decays ///////////////////////////////////////////////////////////////////
#include "G4Decay.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4RadioactiveDecay.hh"
#include "G4Radioactivation.hh"

#include "G4NuclideTable.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"

void CupPhysicsList::ConstructGeneral() {

	auto theParticleIterator = GetParticleIterator();

	// Add Decay Process
	G4Decay *theDecayProcess           = new G4Decay();
	CupDeferTrackProc *theDeferProcess = new CupDeferTrackProc();
	theParticleIterator->reset();
	while ((*theParticleIterator)()) {
		G4ParticleDefinition *particle = theParticleIterator->value();
		G4ProcessManager *pmanager     = particle->GetProcessManager();

		if (theDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()) {
			pmanager->AddProcess(theDecayProcess);
			// set ordering for PostStepDoIt and AtRestDoIt
			pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
			pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
		}
		if (!particle->IsShortLived()) pmanager->AddDiscreteProcess(theDeferProcess);
	}

	//G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();
	auto theRadioactiveDecay  = new G4Radioactivation();
  	theRadioactiveDecay->SetARM(true);    

	//set a finer grid of the physic tables in order to improve precision
  	//former LowEnergy models have 200 bins up to 100 GeV
  	G4EmParameters* emParameters = G4EmParameters::Instance();
  	//emParameters->SetVerbose(2);
  	//emParameters->SetNumberOfBinsPerDecade(20);
  	//emParameters->SetMscStepLimitType(fMinimal);
	// EM physics extra configuration
  	// this physics constructor should be defined after EM constructor
	//emParameters->SetMinEnergy(100.0 * eV);
	//emParameters->SetMaxEnergy(100.0 * GeV);
	emParameters->SetAuger(true);  // Enable Auger electron production
	emParameters->SetFluo(true);   // Enable fluorescence
	emParameters->SetPixe(true);   // Enable PIXE
	emParameters->SetDeexcitationIgnoreCut(true);

	// EM physics constructor is not used in this example, so
  	// it is needed to instantiate and to initialize atomic deexcitation
  	G4LossTableManager* man = G4LossTableManager::Instance();
  	G4VAtomDeexcitation* ad = man->AtomDeexcitation();
  	if(!ad) {
    	G4EmParameters::Instance()->SetAugerCascade(true);
    	ad = new G4UAtomicDeexcitation();
    	man->SetAtomDeexcitation(ad);
    	ad->InitialiseAtomicDeexcitation();
	}

	// Declare radioactive decay to the GenericIon in the IonTable.
	// register radioactiveDecay
  	G4PhysicsListHelper::GetPhysicsListHelper()->
    	RegisterProcess(theRadioactiveDecay, G4GenericIon::GenericIon());

	G4ParticleDefinition* triton = G4Triton::Definition();
    	triton->SetPDGStable(false);
	G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(theRadioactiveDecay,
                                                               G4Triton::Triton());
	
/*
	// Declare radioactive decay to the GenericIon in the IonTable.
	G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();
	//auto theRadioactiveDecay  = new G4Radioactivation();

	const G4IonTable *theIonTable           = G4ParticleTable::GetParticleTable()->GetIonTable();

	G4ParticleDefinition* triton = G4Triton::Definition();
        triton->SetPDGStable(false);

	for (G4int i = 0; i < theIonTable->Entries(); i++) {
		G4String particleName = theIonTable->GetParticle(i)->GetParticleName();

		if (particleName == "GenericIon") {
			G4ProcessManager *pmanager = theIonTable->GetParticle(i)->GetProcessManager();
			pmanager->SetVerboseLevel(VerboseLevel);
			pmanager->AddProcess(theRadioactiveDecay);
			pmanager->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
			pmanager->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
		}

		if (particleName == "triton") {
			G4ProcessManager *pmanager = theIonTable->GetParticle(i)->GetProcessManager();
			pmanager->SetVerboseLevel(VerboseLevel);
			pmanager->AddProcess(theRadioactiveDecay);
			pmanager->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
			pmanager->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
		}
	}

	// Configure atomic deexcitation parameters
    //G4EmParameters* emParameters = G4EmParameters::Instance();
    auto emParameters = G4EmParameters::Instance();
    emParameters->SetVerbose(1);
    emParameters->SetMinEnergy(100.0 * eV);
    emParameters->SetMaxEnergy(100.0 * GeV);
    emParameters->SetAuger(true);  // Enable Auger electron production
    emParameters->SetFluo(true);   // Enable fluorescence
    emParameters->SetPixe(true);   // Enable PIXE
    emParameters->SetDeexcitationIgnoreCut(true); // Ensure deexcitation processes are not cut off
	*/
}

// Cuts /////////////////////////////////////////////////////////////////////
void CupPhysicsList::SetCuts() {

	if (verboseLevel > 1) G4cout << "CupSim/CupPhysicsList::SetCuts:";

	if (verboseLevel > 0) {
		G4cout << "CupSim/CupPhysicsList::SetCuts:";
		G4cout << "CutLength : " << G4BestUnit(defaultCutValue, "Length") << G4endl;
	}

	// special for low energy physics
	// G4double lowlimit=250*eV;
	// G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,100.*GeV);

	// set cut values for gamma at first and for e- second and next for e+,
	// because some processes for e+/e- need cut values for gamma
	SetCutValue(cutForGamma, "gamma");
	SetCutValue(cutForElectron, "e-");
	SetCutValue(cutForPositron, "e+");

	// setcut for region

	G4cout << "PhysicsList: setcut for region  \n";

	if (!DetectorCuts) SetDetectorCut(cutForElectron);
	G4Region *region = (G4RegionStore::GetInstance())->GetRegion("crystals");
	region->SetProductionCuts(DetectorCuts);

	G4cout << ">>SetCuts:endof setcut" << G4endl;

	if (verboseLevel > 0) DumpCutValuesTable();
}

void CupPhysicsList::AddPhysicsList(const G4String &name) {
	G4cout << "\n>>>   CupPhysicsList::AddPhysicsList: <<<" << name << ">>> \n" << G4endl;

	if (name == "livermore") {
		G4cout << "Physics : EmLivermore is selected \n";
		emName = name;
	} else {
		G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
			<< " is not defined" << G4endl;
	}
}

void CupPhysicsList::SetCutForGamma(G4double cut) {
	cutForGamma = cut;
	G4cout << "SetCutForGamma=" << cut << G4endl;
	SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

void CupPhysicsList::SetCutForElectron(G4double cut) {
	cutForElectron = cut;
	SetParticleCuts(cutForElectron, G4Electron::Electron());
}

void CupPhysicsList::SetCutForPositron(G4double cut) {
	cutForPositron = cut;
	SetParticleCuts(cutForPositron, G4Positron::Positron());
}

void CupPhysicsList::SetDetectorCut(G4double cut) {
	if (!DetectorCuts) DetectorCuts = new G4ProductionCuts();

	DetectorCuts->SetProductionCut(cut, idxG4GammaCut);
	DetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
	DetectorCuts->SetProductionCut(cut, idxG4PositronCut);

	G4cout << ">>DetectorCuts are set with " << G4BestUnit(cut, "Length") << G4endl;
}

void CupPhysicsList::List() {}
