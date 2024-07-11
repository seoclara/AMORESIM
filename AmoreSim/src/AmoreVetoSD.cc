#include "AmoreSim/AmoreVetoSD.hh"
#include "AmoreSim/AmoreDetectorConstruction.hh"
#include "CupSim/CupVetoSD.hh"

#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"

#include <fstream>
#include <sstream>
using namespace std;

AmoreVetoSD::AmoreVetoSD(G4String name, int arg_max_tgs) : CupVetoSD(name, arg_max_tgs) {}
AmoreVetoSD::~AmoreVetoSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4bool AmoreVetoSD::ProcessHits(G4Step *aStep, G4TouchableHistory * /*ROhist*/) 
{
    G4EmSaturation *emSaturation = G4LossTableManager::Instance()->EmSaturation();
    //G4double edep_quenched = emSaturation->VisibleEnergyDeposition(aStep);
    G4double edep_quenched = emSaturation->VisibleEnergyDepositionAtAStep(aStep);

    G4double edep = aStep->GetTotalEnergyDeposit();
    G4ParticleDefinition *particleType = aStep->GetTrack()->GetDefinition();
    G4String particleName = particleType->GetParticleName();

    if (edep == 0. || particleName == "opticalphoton") {
        return true;
    }
    // G4cout << "JW: AmoreVetoSD::ProcessHits ---> edep_quenched= " << edep_quenched << G4endl;

    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    G4int copyNo = theTouchable->GetCopyNumber();
    G4int motherCopyNo = theTouchable->GetCopyNumber(1); 
    G4int envelopeCopyNo = theTouchable->GetCopyNumber(2);
    G4VPhysicalVolume *thePhysical = theTouchable->GetVolume();
    G4VPhysicalVolume *theMotherPhysical = theTouchable->GetVolume(1);
    G4VPhysicalVolume *theEnvelopePhysical = theTouchable->GetVolume(2);
    G4String VolName = thePhysical->GetName();
    G4String motherVolName = theMotherPhysical->GetName();
    G4String envelopeVolName = theEnvelopePhysical->GetName();

    // G4cout << "JW: MuonVeto copyNo= " << copyNo << G4endl;
    switch(AmoreDetectorConstruction::GetDetGeometryType()){
        case eDetGeometry::kDetector_AMoRE200:{
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

            copyNo = (strstr(VolName, "PlasticScintO")) ? envelopeCopyNo : envelopeCopyNo + maxPSNo;
            break;
        }
        case eDetGeometry::kDetector_AMoRE_I:{
            copyNo = (strstr(motherVolName, "Envelope")) ? motherCopyNo : copyNo;
        }

    }
    // G4cout << "             VolName= " << VolName << G4endl;
    // G4cout << "             motherCopyNo= " << motherCopyNo << G4endl;
    // G4cout << "             motherVolName= " << motherVolName << G4endl;
    // G4cout << "             envelopeCopyNo= " << envelopeCopyNo << G4endl;
    // G4cout << "             envelopeVolName= " << envelopeVolName << G4endl;
    // G4cout << "      changed CopyNo = " << copyNo << G4endl;

    CupVetoHit *aHit = (*hitsCollection)[copyNo];
    if (!(aHit->GetLogV())) {
        aHit->SetLogV(thePhysical->GetLogicalVolume());
        G4AffineTransform aTrans = theTouchable->GetHistory()->GetTopTransform();
        aTrans.Invert();
        aHit->SetRot(aTrans.NetRotation());
        aHit->SetPos(aTrans.NetTranslation());
    }

    aHit->AddEdep(edep);
    aHit->AddEdepQuenched(edep_quenched);

    // return CupVetoSD::ProcessHits(aStep, ROhist);
    return true;
}
