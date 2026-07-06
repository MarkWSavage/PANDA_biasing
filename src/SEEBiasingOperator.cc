#include "SEEBiasingOperator.hh"

#include "G4BiasingProcessInterface.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessType.hh"
#include "G4HadronicProcessType.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include "G4Threading.hh"

#include <cfloat>

// Static member definitions (shared across every worker thread's
// operator instance -- see SEEBiasingOperator.hh for why).
std::atomic<G4long> SEEBiasingOperator::fsNumProposed{0};
std::atomic<G4long> SEEBiasingOperator::fsNumConfirmedInteractions{0};
std::atomic<G4double> SEEBiasingOperator::fsCrossSectionFactor{1.0};

std::atomic<G4long> SEEBiasingOperator::fsNumProposedSecondary{0};
std::atomic<G4long> SEEBiasingOperator::fsNumConfirmedInteractionsSecondary{0};
std::atomic<G4double> SEEBiasingOperator::fsCrossSectionFactorSecondary{1.0};
std::atomic<bool> SEEBiasingOperator::fsSecondaryRoleUsed{false};

SEEBiasingOperator::SEEBiasingOperator(
    const std::vector<SpeciesBias>& speciesBiases,
    const G4String& name)
: G4VBiasingOperator(name),
  fSpeciesBiases(speciesBiases)
{
    // All threads construct their operator with the same factors (they
    // come from the same run.mac settings), so these writes are racy
    // in principle but harmless in practice -- every thread writes the
    // same values.
    for (const auto& sb : fSpeciesBiases)
    {
        if (sb.isSecondaryRole)
        {
            fsCrossSectionFactorSecondary.store(sb.crossSectionFactor);
            fsSecondaryRoleUsed.store(true);
        }
        else
        {
            fsCrossSectionFactor.store(sb.crossSectionFactor);
        }
    }
}

SEEBiasingOperator::~SEEBiasingOperator()
{
    // NOTE: do not rely on this destructor for diagnostic output --
    // see PrintTotals() and the bug-history comment in the header for
    // why. This destructor may never run in practice (the operator is
    // created with `new` in DetectorConstruction::ConstructSDandField()
    // and never explicitly deleted, matching common Geant4 biasing
    // example patterns).
    for (auto& kv : fChangeCrossSectionOperations)
        delete kv.second.operation;
}

void SEEBiasingOperator::PrintTotals()
{
    G4long proposed = fsNumProposed.load();
    G4long confirmed = fsNumConfirmedInteractions.load();
    G4double factor = fsCrossSectionFactor.load();

    G4cout << G4endl;
    G4cout << "=======================================================" << G4endl;
    G4cout << "SEEBiasingOperator -- GLOBAL FINAL COUNTS (primary)" << G4endl;
    G4cout << "(all worker threads combined; ground truth, independent"
           << " of CSV/Recoil_keV/Proton_keV scoring)" << G4endl;
    G4cout << "    Cross section factor used   : " << factor << G4endl;
    G4cout << "    Proposed biasing operations : " << proposed << G4endl;
    G4cout << "    Confirmed real interactions : " << confirmed << G4endl;
    if (proposed > 0)
    {
        G4cout << "    Confirmed/Proposed ratio    : "
               << (G4double)confirmed / proposed
               << " (expect ~= 1/crossSectionFactor = "
               << 1.0 / factor << ")" << G4endl;
    }
    else
    {
        G4cout << "    WARNING: zero proposed operations -- biasing "
               << "never engaged at all. Check StartRun() output above "
               << "for 'wrapped process(es)' confirmation." << G4endl;
    }
    G4cout << "=======================================================" << G4endl;
    G4cout << G4endl;

    if (fsSecondaryRoleUsed.load())
    {
        G4long proposedSec = fsNumProposedSecondary.load();
        G4long confirmedSec = fsNumConfirmedInteractionsSecondary.load();
        G4double factorSec = fsCrossSectionFactorSecondary.load();

        G4cout << "=======================================================" << G4endl;
        G4cout << "SEEBiasingOperator -- GLOBAL FINAL COUNTS (secondary neutron)" << G4endl;
        G4cout << "    Cross section factor used   : " << factorSec << G4endl;
        G4cout << "    Proposed biasing operations : " << proposedSec << G4endl;
        G4cout << "    Confirmed real interactions : " << confirmedSec << G4endl;
        if (proposedSec > 0)
        {
            G4cout << "    Confirmed/Proposed ratio    : "
                   << (G4double)confirmedSec / proposedSec
                   << " (expect ~= 1/crossSectionFactor = "
                   << 1.0 / factorSec << ")" << G4endl;
        }
        else
        {
            G4cout << "    WARNING: zero proposed operations -- secondary-"
                   << "neutron biasing never engaged at all." << G4endl;
        }
        G4cout << "=======================================================" << G4endl;
        G4cout << G4endl;
    }
}

void SEEBiasingOperator::StartRun()
{
    // Create one change-cross-section operation for every physics
    // process wrapped for biasing, for every species this operator
    // covers (this is populated by
    // G4GenericBiasingPhysics::PhysicsBias(...) called in PANDA.cc,
    // which wraps each species' own inelastic process name -- e.g.
    // "protonInelastic", "neutronInelastic"). Each is a candidate;
    // ProposeOccurenceBiasingOperation below decides at runtime
    // whether a given wrapped process is actually the hadronic
    // inelastic process we want to bias -- all others are left
    // unbiased (nullptr returned for them).
    for (const auto& sb : fSpeciesBiases)
    {
        const G4ProcessManager* processManager =
            sb.particle->GetProcessManager();

        const G4BiasingProcessSharedData* sharedData =
            G4BiasingProcessInterface::GetSharedData(processManager);

        if (!sharedData)
        {
            G4cout << "SEEBiasingOperator::StartRun() -- WARNING: no biasing "
                   << "shared data found for " << sb.particle->GetParticleName()
                   << ". Biasing will have NO EFFECT for this species. Check "
                   << "that PhysicsBias(...) was called for it in PANDA.cc."
                   << G4endl;
            continue;
        }

        G4cout << "SEEBiasingOperator::StartRun() -- wrapped process(es) for "
               << sb.particle->GetParticleName() << " biasing:" << G4endl;

        for (std::size_t i = 0;
             i < sharedData->GetPhysicsBiasingProcessInterfaces().size();
             ++i)
        {
            const G4BiasingProcessInterface* wrapperProcess =
                sharedData->GetPhysicsBiasingProcessInterfaces()[i];

            G4String wrappedName =
                wrapperProcess->GetWrappedProcess()->GetProcessName();

            G4cout << "    - " << wrappedName << G4endl;

            G4String operationName = "XSboost-" + wrappedName;

            fChangeCrossSectionOperations[wrapperProcess] = BiasEntry{
                new G4BOptnChangeCrossSection(operationName),
                sb.crossSectionFactor,
                sb.isSecondaryRole
            };
        }
    }
}

G4VBiasingOperation* SEEBiasingOperator::ProposeOccurenceBiasingOperation(
    const G4Track* /*track*/,
    const G4BiasingProcessInterface* callingProcess)
{
    // Keying by callingProcess alone already selects the right species:
    // each species has its own distinct wrapped-process interface
    // object(s) (populated per-species in StartRun()), so a track of a
    // species this operator doesn't cover simply won't be found here.
    auto it = fChangeCrossSectionOperations.find(callingProcess);
    if (it == fChangeCrossSectionOperations.end())
        return nullptr;

    // Only bias the hadronic INELASTIC process. Elastic scattering
    // and every other process for this particle is left completely
    // unbiased -- returning nullptr restores normal, analog physics
    // for them.
    const G4VProcess* wrapped = callingProcess->GetWrappedProcess();

    if (wrapped->GetProcessType() != fHadronic)
        return nullptr;

    if (wrapped->GetProcessSubType() != fHadronInelastic)
        return nullptr;

    G4double analogInteractionLength =
        wrapped->GetCurrentInteractionLength();

    // Process not active for this step/track (e.g. below threshold,
    // or particle about to leave the biased volume): leave it alone.
    if (analogInteractionLength > DBL_MAX / 10.0)
        return nullptr;

    G4double analogXS = 1.0 / analogInteractionLength;

    BiasEntry& entry = it->second;

    entry.operation->SetBiasedCrossSection(entry.crossSectionFactor * analogXS);
    entry.operation->Sample();

    if (entry.isSecondaryRole)
        ++fsNumProposedSecondary;
    else
        ++fsNumProposed;

    return entry.operation;
}

void SEEBiasingOperator::OperationApplied(
    const G4BiasingProcessInterface* callingProcess,
    G4BiasingAppliedCase /*biasingCase*/,
    G4VBiasingOperation* occurenceOperationApplied,
    G4double /*weightForOccurenceInteraction*/,
    G4VBiasingOperation* /*finalStateOperationApplied*/,
    const G4VParticleChange* /*particleChangeProduced*/)
{
    auto it = fChangeCrossSectionOperations.find(callingProcess);
    if (it == fChangeCrossSectionOperations.end())
        return;

    BiasEntry& entry = it->second;

    if (entry.operation == occurenceOperationApplied)
    {
        entry.operation->SetInteractionOccured();

        if (entry.isSecondaryRole)
            ++fsNumConfirmedInteractionsSecondary;
        else
            ++fsNumConfirmedInteractions;
    }
}
