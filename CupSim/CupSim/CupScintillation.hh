#ifndef CupScintillation_h
#define CupScintillation_h 1

/////////////
// Includes
/////////////

#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleMomentum.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4Poisson.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4VRestDiscreteProcess.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "templates.hh"
#include "G4Version.hh"

#include "G4EmSaturation.hh"
#include "G4UImessenger.hh"

// Class Description:
// RestDiscrete Process - Generation of Scintillation Photons.
// Class inherits publicly from G4VRestDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

// EJ: start
class G4UIcommand;
class G4UIdirectory;

class CupScintillation : public G4VRestDiscreteProcess,
                         public G4UImessenger

{

  public:
    ////////////////////////////////
    // Constructors and Destructor
    ////////////////////////////////

    CupScintillation(const G4String &processName = "Scintillation",
                     G4ProcessType type          = fElectromagnetic);
    ~CupScintillation();

  private:
    CupScintillation(const CupScintillation &right);

    //////////////
    // Operators
    //////////////

    CupScintillation &operator=(const CupScintillation &right);

  public:
    ////////////
    // Methods
    ////////////

    // CupScintillation Process has both PostStepDoIt (for energy
    // deposition of particles in flight) and AtRestDoIt (for energy
    // given to the medium by particles at rest)

    G4bool IsApplicable(const G4ParticleDefinition &aParticleType);
    // Returns true -> 'is applicable', for any particle type except
    // for an 'opticalphoton' and for short-lived particles

    G4double GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *);
    // Returns infinity; i. e. the process does not limit the step,
    // but sets the 'StronglyForced' condition for the DoIt to be
    // invoked at every step.

    G4double GetMeanLifeTime(const G4Track &aTrack, G4ForceCondition *);
    // Returns infinity; i. e. the process does not limit the time,
    // but sets the 'StronglyForced' condition for the DoIt to be
    // invoked at every step.

    G4VParticleChange *PostStepDoIt(const G4Track &aTrack, const G4Step &aStep);
    G4VParticleChange *AtRestDoIt(const G4Track &aTrack, const G4Step &aStep);

    // These are the methods implementing the scintillation process.

    void SetTrackSecondariesFirst(const G4bool state);
    // If set, the primary particle tracking is interrupted and any
    // produced scintillation photons are tracked next. When all
    // have been tracked, the tracking of the primary resumes.

    void SetFiniteRiseTime(const G4bool state);
    // If set, the CupScintillation process expects the user to have
    // set the constant material property FAST/SLOWSCINTILLATIONRISETIME.

    G4bool GetTrackSecondariesFirst() const;
    // Returns the boolean flag for tracking secondaries first.

    G4bool GetFiniteRiseTime() const;
    // Returns the boolean flag for a finite scintillation rise time.

    void SetScintillationYieldFactor(const G4double yieldfactor);
    // Called to set the scintillation photon yield factor, needed when
    // the yield is different for different types of particles. This
    // scales the yield obtained from the G4MaterialPropertiesTable.

    G4double GetScintillationYieldFactor() const;
    // Returns the photon yield factor.

    void SetScintillationExcitationRatio(const G4double excitationratio);
    // Called to set the scintillation exciation ratio, needed when
    // the scintillation level excitation is different for different
    // types of particles. This overwrites the YieldRatio obtained
    // from the G4MaterialPropertiesTable.

    G4double GetScintillationExcitationRatio() const;
    // Returns the scintillation level excitation ratio.

    G4PhysicsTable *GetFastIntegralTable() const;
    // Returns the address of the fast scintillation integral table.

    G4PhysicsTable *GetSlowIntegralTable() const;
    // Returns the address of the slow scintillation integral table.

    void AddSaturation(G4EmSaturation *sat) { emSaturation = sat; }
    // Adds Birks Saturation to the process.

    void RemoveSaturation() { emSaturation = NULL; }
    // Removes the Birks Saturation from the process.

    G4EmSaturation *GetSaturation() const { return emSaturation; }
    // Returns the Birks Saturation.

    void SetScintillationByParticleType(const G4bool);
    // Called by the user to set the scintillation yield as a function
    // of energy deposited by particle type

    G4bool GetScintillationByParticleType() const { return scintillationByParticleType; }
    // Return the boolean that determines the method of scintillation
    // production

    void DumpPhysicsTable() const;
    // Prints the fast and slow scintillation integral tables.

    // EJ: start
    // following two methods are for G4UImessenger
    void SetNewValue(G4UIcommand *command, G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);
    void SetVerboseLevel(G4int level);
    G4int GetVerboseLevel(void) const;
    void PrintVerboseLevel();  // Method to print the current verbose level

    // following are for energy deposition diagnosis
    static void ResetTotEdep() {
        totEdep = totEdep_quenched = 0.0;
        nScintPhotons              = 0;
        // scintCentroidSum *= 0.0;
    }
    static G4double GetTotEdep() { return totEdep; }
    static G4double GetTotEdepQuenched() { return totEdep_quenched; }
    static G4int GetNScintPhotons() { return nScintPhotons; }
    static G4bool GetDoScintillation() { return doScintillation; }
    static G4ThreeVector GetScintCentroid() { return scintCentroidSum * (1.0 / totEdep_quenched); }
    // EJ: end

  protected:
    void BuildThePhysicsTable();
    // It builds either the fast or slow scintillation integral table;
    // or both.

    ///////////////////////
    // Class Data Members
    ///////////////////////

    G4PhysicsTable *theSlowIntegralTable;
    G4PhysicsTable *theFastIntegralTable;

    G4bool fTrackSecondariesFirst;
    G4bool fFiniteRiseTime;

    G4double YieldFactor;

    G4double ExcitationRatio;

    G4bool scintillationByParticleType;

    // private:

    G4double single_exp(G4double t, G4double tau2);
    G4double bi_exp(G4double t, G4double tau1, G4double tau2);

    // emission time distribution when there is a finite rise time
    G4double sample_time(G4double tau1, G4double tau2);

    G4EmSaturation *emSaturation;

    // EJ: start
    static G4int verboseLevel;
    static G4UIdirectory *CupScintDir;
    static G4bool doScintillation;
    static G4ThreadLocal G4double totEdep;
    static G4ThreadLocal G4double totEdep_quenched;
    static G4ThreadLocal G4int nScintPhotons;
    static G4ThreadLocal G4ThreeVector scintCentroidSum;
    // EJ: end
};

////////////////////
// Inline methods
////////////////////

inline G4bool CupScintillation::IsApplicable(const G4ParticleDefinition &aParticleType) {
    if (aParticleType.GetParticleName() == "opticalphoton") return false;
    if (aParticleType.IsShortLived()) return false;

    return true;
}

inline void CupScintillation::SetTrackSecondariesFirst(const G4bool state) {
    fTrackSecondariesFirst = state;
}

inline void CupScintillation::SetFiniteRiseTime(const G4bool state) { fFiniteRiseTime = state; }

inline G4bool CupScintillation::GetTrackSecondariesFirst() const { return fTrackSecondariesFirst; }

inline G4bool CupScintillation::GetFiniteRiseTime() const { return fFiniteRiseTime; }

inline void CupScintillation::SetScintillationYieldFactor(const G4double yieldfactor) {
    YieldFactor = yieldfactor;
}

inline G4double CupScintillation::GetScintillationYieldFactor() const { return YieldFactor; }

inline void CupScintillation::SetScintillationExcitationRatio(const G4double excitationratio) {
    ExcitationRatio = excitationratio;
}

inline G4double CupScintillation::GetScintillationExcitationRatio() const {
    return ExcitationRatio;
}

inline G4PhysicsTable *CupScintillation::GetSlowIntegralTable() const {
    return theSlowIntegralTable;
}

inline G4PhysicsTable *CupScintillation::GetFastIntegralTable() const {
    return theFastIntegralTable;
}

inline void CupScintillation::DumpPhysicsTable() const {
    if (theFastIntegralTable) {
        G4int PhysicsTableSize = theFastIntegralTable->entries();
        G4PhysicsOrderedFreeVector *v;

        for (G4int i = 0; i < PhysicsTableSize; i++) {
            v = (G4PhysicsOrderedFreeVector *)(*theFastIntegralTable)[i];
            v->DumpValues();
        }
    }

    if (theSlowIntegralTable) {
        G4int PhysicsTableSize = theSlowIntegralTable->entries();
        G4PhysicsOrderedFreeVector *v;

        for (G4int i = 0; i < PhysicsTableSize; i++) {
            v = (G4PhysicsOrderedFreeVector *)(*theSlowIntegralTable)[i];
            v->DumpValues();
        }
    }
}

inline G4double CupScintillation::single_exp(G4double t, G4double tau2) {
    return std::exp(-1.0 * t / tau2) / tau2;
}

inline G4double CupScintillation::bi_exp(G4double t, G4double tau1, G4double tau2) {
    return std::exp(-1.0 * t / tau2) * (1 - std::exp(-1.0 * t / tau1)) / tau2 / tau2 *
           (tau1 + tau2);
}

inline void CupScintillation::SetVerboseLevel(int level) { verboseLevel= level; }
inline G4int CupScintillation::GetVerboseLevel(void) const { return verboseLevel; }

#endif /* CupScintillation_h */
