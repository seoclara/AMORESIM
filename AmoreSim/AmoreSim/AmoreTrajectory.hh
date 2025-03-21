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
/// \file runAndEvent/Amore/include/AmoreTrajectory.hh
/// \brief Definition of the AmoreTrajectory class
//
// $Id: $
//
#ifndef AmoreTrajectory_h
#define AmoreTrajectory_h 1

#include "AmoreTrajectoryPoint.hh" // Include from 'tracking'

#include "G4Allocator.hh"
#include "G4ParticleDefinition.hh" // Include from 'particle+matter'
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTrajectory.hh"
#include "G4Version.hh"
#include "G4ios.hh"   // Include from 'system'
#include "globals.hh" // Include from 'global'
#include <stdlib.h>   // Include from 'system'
#include <vector>     // G4RWTValOrderedVector

class G4Polyline; // Forward declaration.
class G4TrackingManager;

typedef std::vector<G4VTrajectoryPoint *> TrajectoryPointContainer;

//
/// User trajectory class
///
/// - new, delete and "==" operators are overwritten
///
/// - get functions
///     G4int GetTrackID() const, G4int GetParentID() const,
///     G4String GetParticleName() const, G4double GetCharge() const,
///     G4int GetPDGEncoding() const, G4double GetInitialKineticEnergy() const
///     and G4ThreeVector GetInitialMomentum() const
///
/// - void ShowTrajectory(std::ostream& os=G4cout) const
///     invokes the default implementation
///
/// - void DrawTrajectory(G4int i_mode = 0) const
///     invokes the default implementation
///
/// - void AppendStep(const G4Step* aStep)
///     adds a user trajectory point object, AmoreTrajectoryPoint
///
/// - int GetPointEntries() const
///     returns the number of point entries
///
/// - G4VTrajectoryPoint* GetPoint(G4int i) const
///     gets the i-th trajectory point
///
/// - void MergeTrajectory(G4VTrajectory* secondTrajectory)
///     adds a trajectory to a TrajectoryPointContainer, fPositionRecord
///
/// - G4ParticleDefinition* GetParticleDefinition()
///     get a particle definition from G4ParticleTable
///
/// - const std::map<G4String,G4AttDef>* GetAttDefs() const
///    defines the track ID, the parent ID, the particle name, the charge,
///    the PDG encoding, the initial kinetic energy, the initial momentum,
///    the initial momentum magnitude and the number of points as attiributes
///
/// - std::vector<G4AttValue>* CreateAttValues() const
///    sets and returns the attributes
//
///////////////////
class AmoreTrajectory : public G4VTrajectory
///////////////////
{

    //--------
  public:
    //--------

    // Constructor/Destrcutor

    AmoreTrajectory();

    AmoreTrajectory(const G4Track *aTrack, G4TrackingManager *aMan);
    AmoreTrajectory(AmoreTrajectory &);
    virtual ~AmoreTrajectory();

    // Operators
    inline void *operator new(size_t);
    inline void operator delete(void *);
    inline int operator==(const AmoreTrajectory &right) const { return (this == &right); }

    // Get/Set functions
    inline virtual G4int GetTrackID() const { return fTrackID; }
    inline virtual G4int GetParentID() const { return fParentID; }
    inline virtual G4String GetParticleName() const { return fParticleName; }
    inline virtual G4double GetCharge() const { return fPDGCharge; }
    inline virtual G4int GetPDGEncoding() const { return fPDGEncoding; }
    inline virtual G4double GetInitialKineticEnergy() const { return fInitialKineticEnergy; }
    inline virtual G4ThreeVector GetInitialMomentum() const { return fInitialMomentum; }

    // Other member functions
    virtual void ShowTrajectory(std::ostream &os = G4cout) const;
    virtual void DrawTrajectory() const;
    virtual void AppendStep(const G4Step *aStep);
    virtual int GetPointEntries() const { return fPositionRecord->size(); }
    virtual G4VTrajectoryPoint *GetPoint(G4int i) const { return (*fPositionRecord)[i]; }
    virtual void MergeTrajectory(G4VTrajectory *secondTrajectory);

    G4ParticleDefinition *GetParticleDefinition();

    virtual const std::map<G4String, G4AttDef> *GetAttDefs() const;
    virtual std::vector<G4AttValue> *CreateAttValues() const;

    //---------
  private:
    //---------

    G4TrackingManager *fpMan;
    TrajectoryPointContainer *fPositionRecord;
    G4int fTrackID;
    G4int fParentID;
    G4int fPDGEncoding;
    G4double fPDGCharge;
    G4String fParticleName;
    G4double fInitialKineticEnergy;
    G4ThreeVector fInitialMomentum;

    G4int prevSecNum;
};

#if G4VERSION_NUMBER <= 999
extern G4Allocator<AmoreTrajectory> *faTrajAllocator;
#else
extern G4ThreadLocal G4Allocator<AmoreTrajectory> *faTrajAllocator;
#endif

inline void *AmoreTrajectory::operator new(size_t) {
    void *aTrajectory;
    if (faTrajAllocator == nullptr) faTrajAllocator = new G4Allocator<AmoreTrajectory>;
    aTrajectory = (void *)faTrajAllocator->MallocSingle();
    return aTrajectory;
}

inline void AmoreTrajectory::operator delete(void *aTrajectory) {
    faTrajAllocator->FreeSingle((AmoreTrajectory *)aTrajectory);
}

#endif
