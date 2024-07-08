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

AmoreSteppingAction::AmoreSteppingAction(AmoreRootNtuple *r) : CupSteppingAction(r){}
AmoreSteppingAction::AmoreSteppingAction(AmoreRootNtuple *r, CupPrimaryGeneratorAction *p)
    : CupSteppingAction(r, p){}

void AmoreSteppingAction::UserSteppingAction(const G4Step *aStep) {
    CupSteppingAction::UserSteppingAction(aStep);
    /*
    const G4Track* trk = aStep->GetTrack();
    G4cout << "PreStep " << aStep->GetPreStepPoint()->GetPosition()
         << " -- PostStep " << aStep->GetPostStepPoint()->GetPosition()
         << " >>> Material= " << trk->GetMaterial()->GetName() << G4endl;
         */
}
