#include "G4RunManagerFactory.hh"
#include "G4Threading.hh"
#include "G4UImanager.hh"
 #include "QGSP_BIC_HP.hh"
//#include "QGSP_INCLXX_HP.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4GenericBiasingPhysics.hh"
#include "SEEBiasingOperator.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "EventAction.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include <algorithm>
#include <vector>

int main(int argc, char** argv)
{
    // Multithreaded run manager: each worker thread gets its own
    // DetectorConstruction::ConstructSDandField() call (already
    // written to expect this -- see its header comment) and its own
    // EventAction with a per-thread events_t<N>.csv (see
    // EventAction::MergeThreadOutputs(), called below once the run
    // finishes, to merge them into the one canonical events.csv).
    // Leave a couple of cores free for the OS/other work rather than
    // claiming every hardware thread.
    G4int nThreads = std::max(
        1, G4Threading::G4GetNumberOfCores() - 2
    );

    auto* runManager =
        G4RunManagerFactory::CreateRunManager(
            G4RunManagerType::MT
        );
    runManager->SetNumberOfThreads(nThreads);

    // Detector geometry
    auto* detector = new DetectorConstruction();
    runManager->SetUserInitialization(detector);

    // Physics list
    auto* physicsList = new QGSP_BIC_HP();
        //new QGSP_INCLXX_HP()

    // QGSP_BIC_HP's own EM constructor (G4EmStandardPhysics_option4)
    // registers GenericIon's ionisation with G4LindhardSorensenIonModel
    // below ~2 MeV/nucleon and G4BetheBlochModel above it -- neither is
    // a per-species heavy-ion stopping-power table. Cross-checked
    // against SRIM (Au197 in Si, 200 MeV): that combination came out
    // 46% low (38.1 vs SRIM's 70.39 MeV*cm2/mg). Passing a non-empty
    // name here flips G4EmStandardPhysics_option4's fUseExternalDEDX
    // flag, swapping in G4IonParametrisedLossModel instead, which uses
    // ICRU73 tables bundled with Geant4 -- confirmed within ~6% of the
    // same SRIM point. Note: ICRU73's own measured data only covers
    // Z=3-18 directly (e.g. z6_14.dat, used as-is for C12); for Z>=19
    // (Au197 included), G4IonDEDXScalingICRU73's default range
    // (minAtomicNumberIon=19) redirects the lookup to iron's (Z=26)
    // table instead, scaled by the ratio of effective (equilibrium)
    // charge^2 between the real ion and Fe at matched velocity -- so
    // z79_14.dat on disk is never actually read for Au in a
    // single-element sensitive material (confirmed via
    // G4IonDEDXHandler::BuildDEDXTable calling AtomicNumberBaseIon()
    // first). The ~6% Au197 agreement above was against that Fe-scaled
    // value, not a native Au table -- and the scaling grows less
    // reliable the further the real ion sits from Fe in Z (e.g. at
    // 16 MeV/u, Au197 vs SRIM opens up to ~10% by this mechanism
    // alone, well above C12's ~8%, which uses native ICRU73 data
    // throughout since Z=6 is below the scaling threshold). Falls back to
    // G4BraggIonModel/G4BetheBlochModel automatically for any
    // ion/material/energy combination ICRU73 doesn't cover, so this is
    // still a strict improvement over the un-swapped default, not a
    // swap with new gaps. Same registered physics name
    // ("G4EmStandard_opt4"), so ReplacePhysics finds and replaces the
    // one QGSP_BIC_HP already added.
    physicsList->ReplacePhysics(new G4EmStandardPhysics_option4(1, "ICRU73"));

    // Note: heavy-ion primaries (e.g. C12, Au197) are NOT pre-created
    // here. G4IonTable::GetIon() needs GenericIon's own process manager
    // already built (it clones that onto the new ion), which physicsList
    // hasn't done yet at this point in main() -- ConstructParticle()/
    // ConstructProcess() only run later, when /run/initialize actually
    // triggers physics construction. See
    // SEEBiasingOperator::StartRun() for where ion creation actually
    // happens instead (runs once per worker thread, confirmed to be the
    // earliest point physics is guaranteed ready -- see the comment in
    // DetectorConstruction::ConstructSDandField() for the full trace).

    // Wrap physics for biasing. Following the official Geant4 biasing
    // examples (GB01/GB07): explicitly name the ONE process to wrap per
    // particle -- its hadronic inelastic process -- rather than
    // wrapping every physics process. This guarantees the ionization
    // process (responsible for essentially all "normal", non-biased
    // energy deposition) is never touched by the biasing framework at
    // all.
    //
    // This registration must happen here, BEFORE run.mac is executed
    // (and therefore before /sim/particle is known), so every particle
    // PANDA might be asked to simulate is wrapped up front. Wrapping a
    // particle that never ends up as the primary is harmless -- its
    // biasing operator (see DetectorConstruction::ConstructSDandField)
    // simply never gets attached/exercised for it.
    //
    // Process names are NOT a consistent "<particle>Inelastic" pattern
    // in Geant4 -- e.g. deuteron is "dInelastic" and triton is
    // "tInelastic", not "deuteronInelastic"/"tritonInelastic". These
    // were taken directly from this physics list's own startup dump
    // ("Hadronic Processes for <particle>"); if QGSP_BIC_HP's process
    // naming ever changes, or you add a species not listed here, check
    // that dump before guessing a name.
    //
    // This does NOT bias anything by itself -- whether biasing
    // actually does anything, and by how much, is entirely controlled
    // by SEEBiasingOperator, created and attached to the relevant
    // logical volumes in DetectorConstruction::ConstructSDandField(),
    // for whichever particle /sim/particle selected, using the cross-
    // section factor set via /sim/biasCrossSectionFactor (default 1.0
    // = no bias).
    auto* biasingPhysics = new G4GenericBiasingPhysics();
    biasingPhysics->PhysicsBias("proton",     {"protonInelastic"});
    biasingPhysics->PhysicsBias("neutron",    {"neutronInelastic"});
    biasingPhysics->PhysicsBias("alpha",      {"alphaInelastic"});
    biasingPhysics->PhysicsBias("deuteron",   {"dInelastic"});
    biasingPhysics->PhysicsBias("triton",     {"tInelastic"});
    biasingPhysics->PhysicsBias("He3",        {"He3Inelastic"});
    biasingPhysics->PhysicsBias("GenericIon", {"ionInelastic"});
    physicsList->RegisterPhysics(biasingPhysics);

    runManager->SetUserInitialization(physicsList);

    // User actions
    runManager->SetUserInitialization(
        new ActionInitialization(detector)
    );

    // Visualization
    auto* visManager = new G4VisExecutive();
    visManager->Initialize();

    // UI manager
    auto* UImanager =
        G4UImanager::GetUIpointer();

    if(argc == 1)
    {
        // Interactive mode
        auto* ui = new G4UIExecutive(argc, argv);
        UImanager->ApplyCommand("/control/execute run.mac");
        ui->SessionStart();
        delete ui;
    }
    else
    {
        // Batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }

    // Print ground-truth biasing interaction counts now, immediately
    // after the run finishes. Deliberately NOT relying on
    // ~SEEBiasingOperator() for this: that destructor may never run
    // (DetectorConstruction::ConstructSDandField() creates each
    // thread's operator with `new` and never deletes it, a common
    // Geant4 biasing example pattern) -- see SEEBiasingOperator.hh for
    // the full explanation. This call site is guaranteed to execute.
    SEEBiasingOperator::PrintTotals();

    // Ground-truth upset tally (weighted count/probability/cross-
    // section at the configured /sim/criticalCharge), accumulated live
    // via static atomics across every worker thread -- an independent
    // C++-side sanity check against PANDA_Analyze.py's post-hoc
    // P(Q>=Qc)/cross-section-at-threshold numbers. Same guaranteed-to-
    // execute call site as PrintTotals() above, for the same reason.
    EventAction::PrintUpsetSummary();

    // Merge each worker thread's events_t<N>.csv into the one
    // canonical events.csv PANDA_Analyze.py/PANDAEX_Analyze.py expect.
    // Safe here: all worker threads have already been joined by this
    // point (G4RunManager cleans them up before ApplyCommand returns),
    // so every per-thread file is fully written and closed.
    EventAction::MergeThreadOutputs();

    // Same merge, for the optional recoil-hits export (/sim/
    // logRecoilHits, default off) -- a no-op if it was never enabled.
    EventAction::MergeRecoilHitsOutputs();

    // Cleanup
    delete visManager;
    delete runManager;

    return 0;
}
