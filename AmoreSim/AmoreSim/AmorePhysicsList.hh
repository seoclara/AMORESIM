#ifndef AmorePhysicsList_h
#define AmorePhysicsList_h 1

#include "CupSim/CupPhysicsList.hh"

#include "G4VUserPhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class AmorePhysicsList : public CupPhysicsList {
  public:
    AmorePhysicsList();
    ~AmorePhysicsList();
    virtual void AddParameterisation();
    virtual void AddPhysicsList(const G4String &name);
    virtual void ConstructProcess();

  private:
    G4String fEMName;
    G4bool emIsRegisted;
    G4VPhysicsConstructor*               emPhysList;
    G4String fOpName;
    G4bool opIsRegisted;
    G4VPhysicsConstructor*               OpPhysList;
    G4String fHadName;
    G4bool hadIsRegisted;
    std::vector<G4VPhysicsConstructor*>  hadronPhys;
};

#endif
