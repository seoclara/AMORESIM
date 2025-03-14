#include "AmoreSim/AmoreVetoSD.hh"
#include "AmoreSim/AmoreDetectorConstruction.hh"
#include "CupSim/CupVetoSD.hh"
#include "AmoreSim/AmoreScintillation.hh"

#include "G4SDManager.hh"
#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"

#include <fstream>
#include <sstream>
using namespace std;

AmoreVetoSD::AmoreVetoSD(G4String name, int arg_max_tgs) : CupVetoSD(name, arg_max_tgs) {
    max_tgs = arg_max_tgs;
    G4String HCname;
    collectionName.insert(HCname = "PSMDColl");
    HCID = -1;
}
AmoreVetoSD::~AmoreVetoSD() {}

void AmoreVetoSD::Initialize(G4HCofThisEvent *HCE){
    // hitsCollection = new CupVetoHitsCollection(SensitiveDetectorName, collectionName[0]);
    hitsCollection = new CupVetoHitsCollection(SensitiveDetectorName, collectionName[1]);
    // G4cout << "JW: AmoreVetoSD::Initialize ---> SensitiveDetectorName= " << SensitiveDetectorName << G4endl;
    // G4cout << "JW: AmoreVetoSD::Initialize ---> collectionName[0]= " << collectionName[0] << G4endl;
    // G4cout << "JW: FullPathName= " << GetFullPathName() << G4endl;
    // G4cout << "JW: collectionName[1]=" << collectionName[1] << G4endl;
    if (HCID < 0) {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection);
    }
    HCE->AddHitsCollection(HCID, hitsCollection);

    for (G4int i = 0; i < max_tgs; i++){
        //CupVetoHit *aHit = new CupVetoHit(i);
        CupVetoHit *aHit = new CupVetoHit();
        hitsCollection->insert(aHit);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool AmoreVetoSD::ProcessHits(G4Step *aStep, G4TouchableHistory * /*ROhist*/) 
{
    // G4cout << "JW:: AmoreVetoSD ProcessHits test" << G4endl;
    G4double edep_quenched             = AmoreScintillation::GetTotEdepQuenched();
    // G4EmSaturation *emSaturation = G4LossTableManager::Instance()->EmSaturation();
    //G4double edep_quenched = emSaturation->VisibleEnergyDeposition(aStep);
    // G4double edep_quenched = emSaturation->VisibleEnergyDepositionAtAStep(aStep);
    G4double edep                      = aStep->GetTotalEnergyDeposit();
    G4ParticleDefinition *particleType = aStep->GetTrack()->GetDefinition();
    G4String particleName              = particleType->GetParticleName();

    // G4cout << "JW: AmoreVetoSD::ProcessHits ---> edep= " << edep << ", edep_quenched= " << edep_quenched << G4endl;
    // G4cout << "JW: AmoreVetoSD::ProcessHits ---> particleName= " << particleName << G4endl;
    if (edep == 0. || particleName == "opticalphoton") return true;
    // if (particleName == "opticalphoton") return true;

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
    // G4cout << "JW: AmoreVetoSD::ProcessHits ---> VolName= " << VolName << G4endl;
    // G4cout << "JW: AmoreVetoSD::ProcessHits ---> motherVolName= " << motherVolName << G4endl;
    // G4cout << "JW: AmoreVetoSD::ProcessHits ---> envelopeVolName= " << envelopeVolName << G4endl;
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
            break;
        }
        default:
            break;

    }
    // G4cout << "JW: MuonVeto copyNo= " << copyNo << G4endl;
    // G4cout << "             VolName= " << VolName << G4endl;
    // G4cout << "             motherCopyNo= " << motherCopyNo << G4endl;
    // G4cout << "             motherVolName= " << motherVolName << G4endl;
    // G4cout << "             envelopeCopyNo= " << envelopeCopyNo << G4endl;
    // G4cout << "             envelopeVolName= " << envelopeVolName << G4endl;
    G4int cellID = -1;
    std::string prefix = "MuonVeto_Envelope";
    if (envelopeVolName.find(prefix) == 0){
        std::string cellIDStr = envelopeVolName.substr(prefix.length());
        cellID = std::stoi(cellIDStr);
        // cellID = (strstr(VolName, "PlasticScintO")) ? std::stoi(cellIDStr) : std::stoi(cellIDStr) + maxPSNo;
        // G4cout << "             cellID= " << cellID << G4endl;
    }

    CupVetoHit *aHit = (*hitsCollection)[copyNo];
    if (!(aHit->GetLogV())) {
        aHit->SetLogV(thePhysical->GetLogicalVolume());
        // G4cout << "             LogV= " << thePhysical->GetLogicalVolume()->GetName() << G4endl;
        G4AffineTransform aTrans = theTouchable->GetHistory()->GetTopTransform();
        aTrans.Invert();
        aHit->SetRot(aTrans.NetRotation());
        aHit->SetPos(aTrans.NetTranslation());
    }

    aHit->AddEdep(edep);
    aHit->AddEdepQuenched(edep_quenched);
    aHit->SetCellID(cellID);

    // return CupVetoSD::ProcessHits(aStep, ROhist);
    return true;
}

void AmoreVetoSD::EndOfEvent(G4HCofThisEvent * /*HCE*/){}
