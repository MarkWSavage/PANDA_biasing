#ifndef SEEBiasingOperator_h
#define SEEBiasingOperator_h 1

#include "G4VBiasingOperator.hh"
#include "G4BOptnChangeCrossSection.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include <map>
#include <vector>
#include <atomic>

// Biases the hadronic INELASTIC process for one or more chosen particle
// species (e.g. protons as PANDA's primary, plus secondary neutrons
// produced by the primary's own reactions) by artificially multiplying
// each species' interaction cross section by its own fixed factor. This
// follows the standard Geant4 "GB01"-style non-splitting cross-section
// biasing recipe: G4VBiasingOperator selects a G4BOptnChangeCrossSection
// operation for the process to be biased; Geant4's biasing framework
// automatically adjusts the interacting track's weight to compensate, so
// the biased simulation remains statistically unbiased in expectation.
//
// IMPORTANT: Geant4 allows only ONE G4VBiasingOperator to be attached to
// a given logical volume at a time (AttachTo() on a second operator for
// an already-attached volume silently fails with a console warning and
// never actually engages). Biasing multiple species in the same volume
// -- e.g. a proton primary plus its secondary neutrons in the sensitive/
// dead/surrounding volumes -- therefore requires ONE operator instance
// covering all of them, not multiple separate instances. Hence this
// class takes a list of (species, factor) pairs rather than just one.
//
// IMPORTANT (bug history -- read before touching EventAction/
// SteppingAction scoring): an earlier version of this codebase
// multiplied every accumulated physical quantity (edep, collected
// charge) by the track's weight directly in C++, on the theory that
// "every scored quantity must be weighted". That is WRONG: the charge
// value itself must stay the RAW, unweighted physical quantity for
// that specific simulated history. PANDA_Analyze.py separately
// applies EventWeight when building the histogram (weights=w in
// np.histogram) -- applying weight in BOTH places double-weights every
// biased event, shrinking its charge value by ~1/factor AND its
// histogram count by another ~1/factor, which silently erases the
// entire nuclear-recoil tail once biasing is active. This is invisible
// in a weight=1.0 sanity check (edep*1.0 == edep) and only manifests
// once real fractional weights exist. The correct design: accumulate
// raw edep/charge in EventAction, track the event's representative
// weight separately via EventAction::UpdateEventWeight(), and let
// weight enter the calculation exactly once, at histogram time in
// Python, using the raw charge value for bin placement.
//
// Only the hadronic inelastic process is biased here (see
// ProposeOccurenceBiasingOperation). Elastic scattering and every
// other process for each biased species is left completely unbiased,
// since elastic proton-nucleus recoils in silicon are far too low in
// energy to meaningfully populate the high-charge SEE tail.
class SEEBiasingOperator : public G4VBiasingOperator
{
public:
    // One species to bias: e.g. { G4Proton::ProtonDefinition(), 1041,
    // false } for PANDA's primary, or { G4Neutron::NeutronDefinition(),
    // 1041, true } for secondary neutrons produced by that primary's
    // own reactions.
    //
    // crossSectionFactor: multiplier applied to the analog (physical)
    //   cross section of this species' hadronic inelastic process,
    //   e.g. 1041 to match CREME-MC's "Hadronic Cross Section
    //   Multiplier". A factor of 1.0 is a no-op and should reproduce
    //   the unbiased baseline -- always verify with factor=1.0 before
    //   trusting a boosted run.
    // isSecondaryRole: true if this species is a SECONDARY particle
    //   (e.g. neutrons produced by a biased proton primary's own
    //   reactions), not the run's primary species. Ground-truth
    //   counters are tracked separately per role (see
    //   fsNumProposedSecondary etc.) so PrintTotals() can report both
    //   without conflating two different species/factors into one
    //   number.
    struct SpeciesBias
    {
        const G4ParticleDefinition* particle;
        G4double crossSectionFactor;
        G4bool isSecondaryRole;
    };

    SEEBiasingOperator(const std::vector<SpeciesBias>& speciesBiases,
                        const G4String& name = "SEEBiasingOperator");
    virtual ~SEEBiasingOperator();

    virtual void StartRun() override;

    // Prints the GLOBAL (all-threads-combined) ground-truth interaction
    // counts to G4cout. Call this explicitly from PANDA.cc right after
    // the run finishes (e.g. after UImanager->ApplyCommand(...) in
    // batch mode), NOT relying on the destructor.
    //
    // BUG HISTORY: an earlier version printed these counts from
    // ~SEEBiasingOperator() and per-instance (non-static) members.
    // That never fired: DetectorConstruction::ConstructSDandField()
    // creates each thread's operator with `new` and never stores or
    // deletes the pointer (a deliberate/common leak pattern in Geant4
    // biasing examples, since operators are meant to live for the
    // run's lifetime) -- so the destructor is simply never called,
    // and the diagnostic silently never printed. Moving to static
    // counters + an explicit call site sidesteps object-lifetime
    // questions entirely.
    static void PrintTotals();

private:
    virtual G4VBiasingOperation*
    ProposeOccurenceBiasingOperation(
        const G4Track* track,
        const G4BiasingProcessInterface* callingProcess) override;

    virtual G4VBiasingOperation*
    ProposeFinalStateBiasingOperation(
        const G4Track*,
        const G4BiasingProcessInterface*) override
    { return nullptr; }

    virtual G4VBiasingOperation*
    ProposeNonPhysicsBiasingOperation(
        const G4Track*,
        const G4BiasingProcessInterface*) override
    { return nullptr; }

    virtual void
    OperationApplied(
        const G4BiasingProcessInterface* callingProcess,
        G4BiasingAppliedCase biasingCase,
        G4VBiasingOperation* occurenceOperationApplied,
        G4double weightForOccurenceInteraction,
        G4VBiasingOperation* finalStateOperationApplied,
        const G4VParticleChange* particleChangeProduced) override;

    std::vector<SpeciesBias> fSpeciesBiases;

    // Hard counters, independent of any downstream event/CSV scoring.
    // STATIC and ATOMIC: shared across every worker thread's operator
    // instance (Geant4 MT creates one SEEBiasingOperator per thread
    // via ConstructSDandField(), all biasing the same species/factors
    // -- these counters give the true combined total across the whole
    // run, read via PrintTotals()). Tracked separately per role
    // (primary vs. secondary) so a secondary species' stats never mix
    // into the primary's totals.
    //
    // fsNumProposed: how many times ProposeOccurenceBiasingOperation
    //   actually returned a real (non-null) biasing operation for a
    //   primary-role species' inelastic process, i.e. how many times
    //   the biasing framework was engaged at all for a step.
    // fsNumConfirmedInteractions: how many of those proposed
    //   operations were subsequently confirmed as an ACTUAL
    //   interaction by OperationApplied() (via
    //   SetInteractionOccured()) -- this is the true count of biased
    //   hadronic-inelastic reactions that occurred in this run, ground
    //   truth independent of Recoil_keV/Proton_keV classification
    //   ambiguity in SteppingAction.cc.
    static std::atomic<G4long> fsNumProposed;
    static std::atomic<G4long> fsNumConfirmedInteractions;
    static std::atomic<G4double> fsCrossSectionFactor;

    // Same three counters, but for secondary-role species.
    // fsSecondaryRoleUsed records whether any secondary-role species
    // was ever configured, so PrintTotals() can omit this section
    // entirely for runs that never used one (the default).
    static std::atomic<G4long> fsNumProposedSecondary;
    static std::atomic<G4long> fsNumConfirmedInteractionsSecondary;
    static std::atomic<G4double> fsCrossSectionFactorSecondary;
    static std::atomic<bool> fsSecondaryRoleUsed;

    // Per-wrapped-process-interface bookkeeping. Each species' process
    // manager has its own distinct G4BiasingProcessInterface object(s)
    // for its wrapped process(es) (populated by
    // G4GenericBiasingPhysics::PhysicsBias(...) in PANDA.cc), so keying
    // by that pointer alone already naturally distinguishes species --
    // a proton track will never present a neutron's callingProcess.
    // Populated in StartRun(), one entry per species per wrapped
    // process; ProposeOccurenceBiasingOperation()/OperationApplied()
    // look up the right operation/factor/role via this map at runtime.
    struct BiasEntry
    {
        G4BOptnChangeCrossSection* operation;
        G4double crossSectionFactor;
        G4bool isSecondaryRole;
    };

    std::map<const G4BiasingProcessInterface*, BiasEntry>
        fChangeCrossSectionOperations;
};

#endif
