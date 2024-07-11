#include <iomanip>

#include "AmoreSim/AmorePhysicsList.hh"
#include "AmoreSim/AmorePhysicsOp.hh"
#include "AmoreSim/PhysListEmStandardNR.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4ThermalNeutrons.hh"

#include "G4ParallelWorldProcess.hh"
#include "G4ProcessManager.hh"
#include "G4TransportationManager.hh"

// Constructor /////////////////////////////////////////////////////////////
AmorePhysicsList::AmorePhysicsList() : CupPhysicsList() {
    emIsRegisted = false;
    opIsRegisted = false;
    hadIsRegisted = false;
}

// Destructor //////////////////////////////////////////////////////////////
AmorePhysicsList::~AmorePhysicsList() {
    delete emPhysList;
    delete OpPhysList;
    for(size_t i=0; i<hadronPhys.size(); i++) {delete hadronPhys[i];}
}

void AmorePhysicsList::ConstructProcess() {

    AddTransportation();

    AddParameterisation();

    // Parallel world process
    G4String pwpname = "ParallelWorldProc";
    G4ParallelWorldProcess* theParallelWorldProcess = new G4ParallelWorldProcess(pwpname);
    theParallelWorldProcess->SetParallelWorld("AmoreParallelWorld");
    theParallelWorldProcess->SetLayeredMaterialFlag(true);

    auto particleIterator = GetParticleIterator();
    particleIterator->reset();
    while( (*particleIterator)() ){
        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName          = particle->GetParticleName();

        pmanager->AddProcess(theParallelWorldProcess, -1,0,0);
        /*
        //if(theParallelWorldProcess->IsAtRestRequired(particle)){
            pmanager->SetProcessOrdering(theParallelWorldProcess, idxAtRest);
            pmanager->SetProcessOrderingToSecond(theParallelWorldProcess, idxAlongStep);
            pmanager->SetProcessOrdering(theParallelWorldProcess, idxPostStep);
        }
        */
    }

    // -- ConstructEM
    if (emIsRegisted) emPhysList->ConstructProcess(); 
    else {
        ConstructEM();
        G4cout << "\n EM Physics is ConstructEM(): default! \n" << G4endl;
    }
    // -- ConstructOp
    if (opIsRegisted) OpPhysList->ConstructProcess();
    else {
        ConstructOp();
        G4cout << "\n Op Physics is ConstructOp(): default! \n" << G4endl;
    }

    // -- ConstructHad
    if (hadIsRegisted) {
	// hadronic physics lists
    G4cout << "\n JW: Had Physics is ConstructHad(): " << hadronPhys.size() << " lists \n" << G4endl;
 	for(size_t i=0; i<hadronPhys.size(); i++) {
        G4cout << "JW:  Had Physics " << i  << ": " << hadronPhys[i]->GetPhysicsName() << G4endl;
   	    hadronPhys[i]->ConstructProcess();
        }
    } else {
    	ConstructHad();
        G4cout << "\n Had Physics is ConstructHad(): default! \n" << G4endl;
    }	

    // ConstructHad();

    ConstructGeneral();
}

// G4FastSimulation Processes //////////////////////////////////////////////////////
#include "G4FastSimulationManagerProcess.hh"
#include "G4ProcessManager.hh"

void AmorePhysicsList::AddParameterisation() {
	auto theParticleIterator = GetParticleIterator();

    	G4cout << "\n>>>   AmorePhysicsList::AddParameterisation: <<< FastSimulationModel >>> \n" << G4endl;

	G4FastSimulationManagerProcess *theFastSimulationManagerProcess =
		//new G4FastSimulationManagerProcess( "G4FSMP" );
		new G4FastSimulationManagerProcess( "G4FSMP", "AmoreParallelWorld");
	theParticleIterator->reset();
	while ((*theParticleIterator)()) {
		G4ParticleDefinition *particle = theParticleIterator->value();
		G4ProcessManager *pmanager     = particle->GetProcessManager();
        //pmanager->AddDiscreteProcess( theFastSimulationManagerProcess );   // No parallel geometry
        pmanager->AddProcess(theFastSimulationManagerProcess, -1, 0, 0);  // General
	}
}

// Add Physics List ////////////////////////////////////////////////////////
void AmorePhysicsList::AddPhysicsList(const G4String &name) {
    G4cout << "\n>>>   AmorePhysicsList::AddPhysicsList: <<<" << name << ">>> \n" << G4endl;

    // EM physics
    if (name == "livermore") {
        fEMName = name;
	emPhysList = new G4EmLivermorePhysics();
	emIsRegisted = true;
    } else if (name == "emstandardNR") {
        fEMName = name;
	emPhysList = new PhysListEmStandardNR();
	emIsRegisted = true;
    // Op physics
    } else if (name == "amorephysicsOp") {
        fOpName = name;
	OpPhysList = new AmorePhysicsOp();
	opIsRegisted = true;
    // Had physics
    // } else if (name == "amorephysicsHad" && !hadIsRegisted) {
    } else if (name == "amorephysicsHad") {
        fHadName = name;
	hadronPhys.push_back(new G4HadronElasticPhysicsHP());
	// hadronPhys.push_back(new G4HadronPhysicsQGSP_BERT_HP());
    hadronPhys.push_back(new G4HadronPhysicsFTFP_BERT_HP());
	hadronPhys.push_back(new G4ThermalNeutrons(0)); 
	hadIsRegisted = true;
    } else {
        G4cout << "AmorePhysicsList::AddPhysicsList: <" << name << ">"
               << " is not defined" << G4endl;
    }
}