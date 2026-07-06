#ifndef EventAction_h
#define EventAction_h 1

#include "G4SystemOfUnits.hh"
#include "G4UserEventAction.hh"
#include "G4GenericMessenger.hh"
#include "globals.hh"
#include "Hit.hh"

#include <vector>
#include <fstream>

class EventAction : public G4UserEventAction
{
public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    // Energy accumulation. These now accumulate RAW, UNWEIGHTED
    // physical energy -- the true value for this specific simulated
    // history. Do NOT multiply by track weight here.
    //
    // IMPORTANT (bug history): an earlier version multiplied edep by
    // weight before summing into fTotalEdep/etc. That is WRONG under
    // cross-section biasing: PANDA_Analyze.py separately applies
    // EventWeight again when building the histogram (weights=w). If
    // the charge value itself is already scaled by weight, a rare
    // biased interaction with weight ~1/1000 gets its charge shrunk
    // by ~1000x AND its histogram count shrunk by ~1000x -- it lands
    // in the wrong (far too low) charge bin instead of contributing
    // rarely-but-correctly to the true high-charge bin. This silently
    // erases the nuclear-recoil tail once biasing is active, while
    // being completely invisible in a weight=1.0 sanity check (since
    // edep*1.0 == edep). The fix: accumulate raw edep here; track the
    // event's representative weight separately via UpdateEventWeight,
    // and let Python apply that weight exactly once, to the histogram
    // count, using the raw charge value for bin placement.
    void AddProtonEdep(G4double edep);
    void AddElectronEdep(G4double edep);
    void AddRecoilEdep(G4double edep);
    void AddCollectedCharge(G4double q);

    // Tracks the representative weight for this event, used to build
    // a correctly weighted cross-section histogram in
    // PANDA_Analyze.py (see events.csv EventWeight column). Call once
    // per step from SteppingAction, passing that step's raw edep and
    // the depositing track's current weight. The event's weight is
    // taken from whichever single step contributed the most edep so
    // far -- exact when at most one biased hadronic interaction
    // drives the event's dominant energy deposit (the overwhelmingly
    // common case here, since even after 1000x-scale biasing these
    // reactions remain rare per event). If multiple independently-
    // biased interactions ever contribute comparably within a single
    // event, this simplification should be revisited.
    void UpdateEventWeight(G4double edep, G4double weight);

    // Event hit storage
    void AddHit(const Hit& hit);

    // Controls
    void SetCriticalCharge(G4double qc);
    void SetVerbose(G4bool val);

    // Geant4 MT creates one EventAction per worker thread, each with
    // its own events.csv writer (see constructor) -- required, since
    // multiple threads truncating/writing the same file concurrently
    // would corrupt it. Call this exactly once, after the run
    // finishes and all worker threads have been joined (e.g.
    // alongside SEEBiasingOperator::PrintTotals() in PANDA.cc), to
    // merge every thread's events_t<N>.csv into the one canonical
    // events.csv PANDA_Analyze.py/PANDAEX_Analyze.py expect, then
    // delete the per-thread files.
    static void MergeThreadOutputs();

    // Same per-thread-file-then-merge pattern as MergeThreadOutputs(),
    // for the optional recoil-hits export (see /sim/logRecoilHits).
    // Call alongside MergeThreadOutputs() in PANDA.cc.
    static void MergeRecoilHitsOutputs();

private:
    std::vector<Hit> fHits;

// Collected charge accumulator
    G4double fCollectedCharge = 0.0;

// Energy totals
    G4double fTotalEdep = 0.0;
    G4double fProtonEdep = 0.0;
    G4double fElectronEdep = 0.0;
    G4double fRecoilEdep = 0.0;

    // See UpdateEventWeight() above.
    G4double fEventWeight = 1.0;
    G4double fMaxSingleEdep = 0.0;

    // Upset logic (G4double, not G4int: accumulates event weights,
    // which are fractional under cross-section biasing)
    G4double fUpsetCount = 0.0;
    G4double fCriticalCharge = 150.0 * 1.0e-15 * CLHEP::coulomb;

    // Output + controls
    G4GenericMessenger* fMessenger = nullptr;
    std::ofstream fCSV;
    G4bool fVerbose = false;

    // Opt-in export of per-hit recoil data (species, Z/A, LET, position)
    // to recoil_hits.csv, for studying the recoil-species/LET spectrum
    // rather than just the aggregate Proton/Electron/Recoil_keV sums --
    // see /sim/logRecoilHits. Off by default: writing every qualifying
    // step to a second file adds real overhead most runs don't need.
    G4bool fLogRecoilHits = false;
    std::ofstream fRecoilHitsCSV;
};

#endif
