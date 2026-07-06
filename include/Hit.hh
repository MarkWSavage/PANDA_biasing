#ifndef Hit_h
#define Hit_h 1

#include "G4ThreeVector.hh"
#include "globals.hh"

struct Hit
{
    G4ThreeVector pos;   // deposition position
    G4double edep;       // raw deposited energy
    G4double time;       // global time

    G4String particle;   // particle species
    G4String process;   // NEW
    G4int trackID;       // unique track
    G4int parentID;      // ancestry
    G4int stepNumber;    // step index

    // Recoil-hit export fields (see EventAction::MergeRecoilHitsOutputs()
    // and /sim/logRecoilHits). Populated for every hit regardless of
    // species -- filtering to recoils only (excluding proton/e-) happens
    // at export time in EventAction, not here.
    G4double stepLength; // this step's path length (native G4 units)
    G4int z;              // atomic number (0 for non-ion species)
    G4int a;              // mass number (0 for non-ion species)
    G4double let;         // LET in MeV*cm2/mg, precomputed using the
                          // hit volume's material density
    G4double weight;      // track weight at time of hit (bias correction,
                          // same value as passed to UpdateEventWeight)
};

#endif
