/**@file amoresim.cc
  Amore Simulation main program (amoresim.cc)
*/
#include "G4Version.hh"
#include "G4MTRunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"
#include "G4VisExecutive.hh"

#include "AmoreSim/AmoreActionInitialization.hh"
#include "AmoreSim/AmorePhysicsList.hh"
#include "AmoreSim/AmoreParallelWorldConstruction.hh"
#include "AmoreSim/AmoreDetectorConstruction.hh"
//#include "AmoreSim/AmorePLManager.hh"
#include "AmoreSim/AmoreRootNtuple.hh"
#include "CupSim/CupDebugMessenger.hh"
#include "CupSim/CupParam.hh"

#include "AmoreSim/AmoreSimGitRevision.hh"
#include "CupSim/CupSimGitRevision.hh"
#include "MCObjs/MCObjsGitRevision.hh"

#ifdef G4VIS_USE
#include "CupSim/CupVisMessenger.hh"
#include "G4VisExecutive.hh"
#endif
#include "G4Run.hh"
#include "Rtypes.h"
#include <cstdlib>


#define STRINGFY(X) #X
#define TOSTRING(X) STRINGFY(X)

using namespace std;

int main(int argc, char **argv) {
    ROOT::EnableThreadSafety();

    cout << "Version information for amoresim:" << endl;
    cout << "AmoreSim executable: "
         << " <Branch: " << TOSTRING(AmoreSim_GIT_BRANCH)
         << ", Revision hash: " << TOSTRING(AmoreSim_GIT_COMMIT_HASH) << ">" << endl;
    AmoreSimGitRevision::PrintGitInfo();
    CupSimGitRevision::PrintGitInfo();
    MCObjsGitRevision::PrintGitInfo();
    cout << endl;
    if (argc == 2 && strcmp(argv[1], "git") == 0) return 0;
	
	// Run manager
    // auto* theRunManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
    auto* theRunManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    // auto theRunManager = new G4MTRunManager;
    theRunManager->SetNumberOfThreads(1);
    // auto theRunManager = new G4RunManager;

    // -- database
    CupParam &db(CupParam::GetDB());
    if (getenv("AmoreDATA") != NULL)
        db.ReadFile((G4String(getenv("AmoreDATA")) + "/settings.dat").c_str());
    else
        db.ReadFile("data/settings.dat");

    // UserInitialization classes
    AmoreDetectorConstruction *theAmoreDetectorConstruction = new AmoreDetectorConstruction;
    theRunManager->SetUserInitialization(theAmoreDetectorConstruction);

    // Parallel world
    G4String parallelWorldName = "AmoreParallelWorld";
    AmoreParallelWorldConstruction *theParallelWorld = new AmoreParallelWorldConstruction(parallelWorldName);
    theAmoreDetectorConstruction->RegisterParallelWorld(theParallelWorld); 	
    //theRunManager->SetUserInitialization(new AmorePhysicsList());

    theRunManager->SetUserInitialization(new AmorePhysicsList());

    // Create the AmoreRecorderBase object
    AmoreRootNtuple *myRecords = new AmoreRootNtuple; // EJ

    // ActionInitialization
    theRunManager->SetUserInitialization(
			new AmoreActionInitialization(myRecords, theAmoreDetectorConstruction));

    // an additional "messenger" class for user diagnostics
    CupDebugMessenger theDebugMessenger(theAmoreDetectorConstruction);

    // Visualization, only if you choose to have it!
#ifdef G4VIS_USE
    G4VisManager *theVisManager      = new G4VisExecutive();
    CupVisMessenger *theVisMessenger = new CupVisMessenger(theVisManager);
    theVisManager->Initialize();
#endif
    // user interface
    G4UImanager *theUI = G4UImanager::GetUIpointer();

    // interactive or batch according to command-line args
    if (argc == 1) {
        // G4UIterminal is a (dumb) terminal.
        // ..but it can be made smart by adding a "shell" to it
        G4UIsession *theSession = new G4UIterminal(new G4UItcsh);
        theSession->SessionStart();
        delete theSession;
    } else { // Batch mode, with optional user interaction
        if (strcmp(argv[1], "gui.mac") == 0 || strcmp(argv[1], "gui") == 0) {
            G4UIExecutive *theSession = new G4UIExecutive(argc, argv);
            theUI->ApplyCommand("/control/execute init_vis.mac");
            if (theSession->IsGUI()) {
                theUI->ApplyCommand("/control/execute gui.mac");
            }

            theSession->SessionStart();
            delete theSession;
        } else {
            G4String command = "/control/execute ";
            for (int iarg = 1; iarg < argc; iarg++)
            {
                G4String fileName = argv[iarg];
                theUI->ApplyCommand(command + fileName);
            }
        }
    }

#ifdef G4VIS_USE
    delete theVisManager;
    delete theVisMessenger;
#endif

    myRecords->CloseFile();

    delete theRunManager;
    delete myRecords; // EJ

    return 0;
}
