//
//  Concrete implementation of G4UserSteppingAction
//

#include "AmoreSim/AmoreSteppingAction.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CupSim/CupPrimaryGeneratorAction.hh"
#include "CupSim/CupRecorderBase.hh" // EJ
#include "CupSim/CupScintillation.hh"
#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VSolid.hh"
#include "G4VisExtent.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4ParallelWorldProcess.hh"

AmoreSteppingAction::AmoreSteppingAction(AmoreRootNtuple *r) : CupSteppingAction(r){}
AmoreSteppingAction::AmoreSteppingAction(AmoreRootNtuple *r, CupPrimaryGeneratorAction *p)
    : CupSteppingAction(r, p){}

void AmoreSteppingAction::UserSteppingAction(const G4Step *aStep) {
    CupSteppingAction::UserSteppingAction(aStep);
    G4Track* trk = aStep->GetTrack();
    const G4Step* pStep = &(*aStep);
    const G4Step* hStep = G4ParallelWorldProcess::GetHyperStep();
    if(hStep != nullptr) pStep = hStep;

    if(trk->GetDefinition()->GetParticleName() == "opticalphoton") {
        if(pStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary){
    	    G4VPhysicalVolume* thePrePV  = pStep->GetPreStepPoint()->GetPhysicalVolume();
    	    G4VPhysicalVolume* thePostPV = pStep->GetPostStepPoint()->GetPhysicalVolume();
    	    G4String prevolName = thePrePV->GetName();
    	    G4String postvolName = thePostPV->GetName();
    	    G4String prematName = thePrePV->GetLogicalVolume()->GetMaterial()->GetName();
    	    G4String postmatName = thePostPV->GetLogicalVolume()->GetMaterial()->GetName();
    	    if (prematName == "PMT_Vac" && postmatName == "PMT_Vac") {
              	trk->SetTrackStatus(fStopAndKill);
        	}
        }
    }
    G4VPhysicalVolume*  trkVolume = aStep->GetTrack()->GetNextVolume();
    if(trkVolume->GetName() == "physWorld") {
        trk->SetTrackStatus(fStopAndKill);
        return;
    }    
    /*
    const G4Track* trk = aStep->GetTrack();
    G4cout << "PreStep " << aStep->GetPreStepPoint()->GetPosition()
         << " -- PostStep " << aStep->GetPostStepPoint()->GetPosition()
         << " >>> Material= " << trk->GetMaterial()->GetName() << G4endl;
         */
}
