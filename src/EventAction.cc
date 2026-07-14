#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "EventAction.hh"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4ios.hh"

namespace {
    const char* kResultsDir = "Results/Current";

    // Geant4 MT's master thread has G4GetThreadId() == -1; a plain
    // (non-MT) G4RunManager never spawns worker threads at all, so
    // this constructor only ever runs with threadID == -1 in that
    // case -- either way, "events_t-1.csv" is a valid, unique
    // filename, so no special-casing is needed for serial vs MT.
    std::string ThreadEventsPath()
    {
        return std::string(kResultsDir) + "/events_t"
             + std::to_string(G4Threading::G4GetThreadId()) + ".csv";
    }

    std::string ThreadRecoilHitsPath()
    {
        return std::string(kResultsDir) + "/recoil_hits_t"
             + std::to_string(G4Threading::G4GetThreadId()) + ".csv";
    }
}

EventAction::EventAction()
: G4UserEventAction()
{
    fMessenger = new G4GenericMessenger(this, "/sim/", "Simulation control");

    fMessenger->DeclareMethod(
        "criticalCharge",
        &EventAction::SetCriticalCharge,
        "Set critical charge in fC");

    fMessenger->DeclareProperty(
        "verbose",
        fVerbose,
        "Enable verbose output");

    fMessenger->DeclareProperty(
        "logRecoilHits",
        fLogRecoilHits,
        "Export per-hit recoil-species/LET data to recoil_hits.csv "
        "(Particle, Z, A, LET_MeV_cm2_mg, position, EventWeight) for "
        "hits in the sensitive volume, excluding proton/e- steps. Off "
        "by default -- adds real per-step file-write overhead most "
        "runs don't need. See EventAction::MergeRecoilHitsOutputs()."
    );

    std::filesystem::create_directories(kResultsDir);

    fCSV.open(ThreadEventsPath());

    fCSV << "EventID,"
         << "DepositedCharge_fC,"
         << "CollectedCharge_fC,"
         << "Total_keV,"
         << "Proton_keV,"
         << "Electron_keV,"
         << "Recoil_keV,"
         << "UpsetCharge_fC,"
         << "EventWeight"
         << G4endl;
}

EventAction::~EventAction()
{
    if (fCSV.is_open())
        fCSV.close();

    if (fRecoilHitsCSV.is_open())
        fRecoilHitsCSV.close();

    delete fMessenger;
}

void EventAction::BeginOfEventAction(const G4Event*)
{
    fHits.clear();

    fTotalEdep    = 0.0;
    fProtonEdep   = 0.0;
    fElectronEdep = 0.0;
    fRecoilEdep   = 0.0;
    fCollectedCharge = 0.0;

    fEventWeight = 1.0;
    fMaxSingleEdep = 0.0;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    auto detector =
        static_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()
                ->GetUserDetectorConstruction()
        );

    // Charge equivalent of the RAW deposited energy (ideal case: every
    // e-h pair fully collected, no trapping/recombination loss).
    // Pair-creation energy depends on the sensitive volume's material --
    // see DetectorConstruction::GetSensitivePairCreationEnergy().
    G4double depositedCharge =
        (fTotalEdep / detector->GetSensitivePairCreationEnergy()) * CLHEP::eplus;

    G4double depositedCharge_fC =
        depositedCharge / (1.0e-15 * CLHEP::coulomb);

    // Charge equivalent ACTUALLY collected at the electrode, after the
    // drift/trapping model in SteppingAction (always computed there,
    // regardless of which one is used below for the upset criterion).
    G4double collectedCharge_fC =
        fCollectedCharge / (1.0e-15 * CLHEP::coulomb);

    // Which of the two charges is used to decide whether this event
    // causes an upset is controlled by /sim/useCollectionModel.
    G4double upsetCharge =
        detector->GetUseCollectionModel()
            ? fCollectedCharge
            : depositedCharge;

    if (upsetCharge >= fCriticalCharge)
        fUpsetCount += fEventWeight;

    if (fCSV.is_open())
    {
        fCSV
            << event->GetEventID() << ","
            << depositedCharge_fC << ","
            << collectedCharge_fC << ","
            << fTotalEdep / keV << ","
            << fProtonEdep / keV << ","
            << fElectronEdep / keV << ","
            << fRecoilEdep / keV << ","
            << upsetCharge / (1.0e-15 * CLHEP::coulomb) << ","
            << fEventWeight
            << G4endl;
    }

    if (fVerbose)
    {
        G4cout
            << "Event " << event->GetEventID()
            << " | Hits: " << fHits.size()
            << " | Total: " << fTotalEdep / keV << " keV"
            << " | Proton: " << fProtonEdep / keV << " keV"
            << " | Electron: " << fElectronEdep / keV << " keV"
            << " | Recoil: " << fRecoilEdep / keV << " keV"
            << " | Deposited Q: " << depositedCharge_fC << " fC"
            << " | Collected Q: " << collectedCharge_fC << " fC"
            << " | Upset: "
            << (upsetCharge >= fCriticalCharge ? 1 : 0)
            << G4endl;
    }

    if (fLogRecoilHits)
    {
        // Lazy-open: /sim/logRecoilHits is applied via messenger AFTER
        // this object's constructor runs (same MT command-ordering
        // quirk documented for /sim/verbose etc.), so the flag isn't
        // known yet at construction time. Opening here, on first actual
        // use, means recoil_hits_t<N>.csv is only ever created for
        // threads/runs that actually enabled it.
        if (!fRecoilHitsCSV.is_open())
        {
            fRecoilHitsCSV.open(ThreadRecoilHitsPath());

            fRecoilHitsCSV
                << "EventID,"
                << "Particle,"
                << "Z,"
                << "A,"
                << "Edep_keV,"
                << "StepLength_um,"
                << "LET_MeV_cm2_mg,"
                << "TrackID,"
                << "ParentID,"
                << "Process,"
                << "Position_X_um,"
                << "Position_Y_um,"
                << "Position_Z_um,"
                << "EventWeight"
                << G4endl;
        }

        for (const auto& hit : fHits)
        {
            // Recoils only, by default -- proton/e- steps are frequent
            // and low-value for LET/recoil-species analysis, and would
            // dominate this file's size for no benefit. See the design
            // discussion this feature came from.
            if (hit.particle == "proton" || hit.particle == "e-")
                continue;

            fRecoilHitsCSV
                << event->GetEventID() << ","
                << hit.particle << ","
                << hit.z << ","
                << hit.a << ","
                << hit.edep / keV << ","
                << hit.stepLength / um << ","
                << hit.let << ","
                << hit.trackID << ","
                << hit.parentID << ","
                << hit.process << ","
                << hit.pos.x() / um << ","
                << hit.pos.y() / um << ","
                << hit.pos.z() / um << ","
                << hit.weight
                << G4endl;
        }
    }
}

void EventAction::AddHit(const Hit& hit)
{
    fHits.push_back(hit);
}

void EventAction::AddProtonEdep(G4double edep)
{
    fProtonEdep += edep;
    fTotalEdep += edep;
}

void EventAction::AddElectronEdep(G4double edep)
{
    fElectronEdep += edep;
    fTotalEdep += edep;
}

void EventAction::AddRecoilEdep(G4double edep)
{
    fRecoilEdep += edep;
    fTotalEdep += edep;
}

void EventAction::UpdateEventWeight(G4double edep, G4double weight)
{
    if (edep > fMaxSingleEdep)
    {
        fMaxSingleEdep = edep;
        fEventWeight = weight;
    }
}

void EventAction::SetCriticalCharge(G4double qc)
{
    fCriticalCharge = qc * 1.0e-15 * CLHEP::coulomb;
}

void EventAction::SetVerbose(G4bool val)
{
    fVerbose = val;
}

void EventAction::AddCollectedCharge(G4double q)
{
    fCollectedCharge += q;
}

void EventAction::MergeThreadOutputs()
{
    namespace fs = std::filesystem;

    const fs::path dir(kResultsDir);
    const std::string prefix = "events_t";
    const std::string suffix = ".csv";

    std::vector<fs::path> threadFiles;

    for (const auto& entry : fs::directory_iterator(dir))
    {
        const std::string name = entry.path().filename().string();

        if (name.size() > prefix.size() + suffix.size()
            && name.compare(0, prefix.size(), prefix) == 0
            && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
        {
            threadFiles.push_back(entry.path());
        }
    }

    // events.csv has 9 columns (see the header written in the
    // constructor). A per-thread file can end up with a torn/glued
    // line -- e.g. two overlapping PANDA runs both pointed at
    // kResultsDir at once, racing unsynchronized writes onto the same
    // events_t<N>.csv path -- so check field count here rather than
    // trust it's well-formed; a malformed line otherwise silently
    // breaks every downstream pandas read_csv() with an opaque
    // "Expected 9 fields" C-parser error, far from this file.
    constexpr int kExpectedFields = 9;

    std::ofstream merged(dir / "events.csv");
    bool headerWritten = false;
    int droppedLines = 0;

    for (const auto& file : threadFiles)
    {
        std::ifstream in(file);
        std::string line;
        bool firstLine = true;

        while (std::getline(in, line))
        {
            if (firstLine)
            {
                firstLine = false;

                if (headerWritten)
                    continue; // skip repeated header from later files

                headerWritten = true;
                merged << line << "\n";
                continue;
            }

            if (std::count(line.begin(), line.end(), ',') + 1 != kExpectedFields)
            {
                ++droppedLines;
                continue;
            }

            merged << line << "\n";
        }
    }

    merged.close();

    if (droppedLines > 0)
    {
        G4cerr << "EventAction::MergeThreadOutputs() -- WARNING: dropped "
               << droppedLines << " malformed line(s) while merging "
               << "(wrong field count -- check for another PANDA run "
               << "having written to " << dir.string()
               << " concurrently)" << G4endl;
    }

    for (const auto& file : threadFiles)
        fs::remove(file);

    G4cout << "EventAction::MergeThreadOutputs() -- merged "
           << threadFiles.size() << " per-thread file(s) into "
           << (dir / "events.csv").string() << G4endl;
}

void EventAction::MergeRecoilHitsOutputs()
{
    namespace fs = std::filesystem;

    const fs::path dir(kResultsDir);
    const std::string prefix = "recoil_hits_t";
    const std::string suffix = ".csv";

    std::vector<fs::path> threadFiles;

    for (const auto& entry : fs::directory_iterator(dir))
    {
        const std::string name = entry.path().filename().string();

        if (name.size() > prefix.size() + suffix.size()
            && name.compare(0, prefix.size(), prefix) == 0
            && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
        {
            threadFiles.push_back(entry.path());
        }
    }

    // Unlike events.csv, this file is opt-in (/sim/logRecoilHits,
    // default off) -- no per-thread files existing at all is the
    // normal case for a run that never enabled it, not an error.
    if (threadFiles.empty())
        return;

    // See the matching field-count check in MergeThreadOutputs() --
    // same rationale, this file's header has 14 columns.
    constexpr int kExpectedFields = 14;

    std::ofstream merged(dir / "recoil_hits.csv");
    bool headerWritten = false;
    int droppedLines = 0;

    for (const auto& file : threadFiles)
    {
        std::ifstream in(file);
        std::string line;
        bool firstLine = true;

        while (std::getline(in, line))
        {
            if (firstLine)
            {
                firstLine = false;

                if (headerWritten)
                    continue;

                headerWritten = true;
                merged << line << "\n";
                continue;
            }

            if (std::count(line.begin(), line.end(), ',') + 1 != kExpectedFields)
            {
                ++droppedLines;
                continue;
            }

            merged << line << "\n";
        }
    }

    merged.close();

    if (droppedLines > 0)
    {
        G4cerr << "EventAction::MergeRecoilHitsOutputs() -- WARNING: dropped "
               << droppedLines << " malformed line(s) while merging "
               << "(wrong field count -- check for another PANDA run "
               << "having written to " << dir.string()
               << " concurrently)" << G4endl;
    }

    for (const auto& file : threadFiles)
        fs::remove(file);

    G4cout << "EventAction::MergeRecoilHitsOutputs() -- merged "
           << threadFiles.size() << " per-thread file(s) into "
           << (dir / "recoil_hits.csv").string() << G4endl;
}
