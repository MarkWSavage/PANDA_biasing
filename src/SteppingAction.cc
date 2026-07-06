#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"

SteppingAction::SteppingAction(
    DetectorConstruction* detector,
    EventAction* eventAction)
: G4UserSteppingAction(),
  fDetector(detector),
  fEventAction(eventAction)
{
}

SteppingAction::~SteppingAction()
{
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    // Energy deposited this step
    G4double edep = step->GetTotalEnergyDeposit();

    if (edep <= 0.0)
        return;

    // Current volume
    auto volume =
        step->GetPreStepPoint()
            ->GetTouchableHandle()
            ->GetVolume();

    if (!volume)
        return;

    // Only score inside sensitive volume
    if (volume->GetLogicalVolume() != fDetector->GetSensitiveLogical())
        return;

    // Particle name
    G4String pname =
        step->GetTrack()
            ->GetParticleDefinition()
            ->GetParticleName();

    // Current track weight. Always 1.0 unless cross-section biasing
    // (SEEBiasingOperator) is active for this track, in which case
    // Geant4's biasing framework has already adjusted it to keep the
    // simulation statistically unbiased in expectation.
    //
    // IMPORTANT: this weight is tracked separately (via
    // UpdateEventWeight) and applied exactly once, at histogram time
    // in PANDA_Analyze.py -- it must NOT be multiplied into edep here.
    // Doing so double-applies the weight (once to the physical charge
    // value, again to the histogram count in Python), which silently
    // erases the nuclear-recoil tail once biasing is active. See the
    // warning in EventAction.hh for the full explanation.
    G4double weight = step->GetTrack()->GetWeight();

    // Record this step's contribution to the event's representative
    // weight (used only for histogram weighting in Python, never for
    // scaling the physical charge/energy values below).
    fEventAction->UpdateEventWeight(edep, weight);

    // Build hit record
    Hit hit;
    hit.edep = edep;
    hit.particle = pname;
    hit.trackID =
        step->GetTrack()->GetTrackID();
    hit.parentID =
        step->GetTrack()->GetParentID();
    hit.stepNumber =
        step->GetTrack()->GetCurrentStepNumber();

    // Recoil-hit export fields -- see Hit.hh. GetAtomicNumber()/
    // GetAtomicMass() are G4ParticleDefinition members (not G4Ions-
    // specific): they return 0 for proton/e-/etc, which is fine since
    // those species are filtered out at export time anyway.
    hit.stepLength = step->GetStepLength();
    hit.z = step->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
    hit.a = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();
    hit.weight = weight;

    if (hit.stepLength > 0.0)
    {
        G4double dEdx = edep / hit.stepLength;
        G4double density =
            fDetector->GetSensitiveLogical()->GetMaterial()->GetDensity();
        hit.let = (dEdx / density) / (MeV * cm2 / mg);
    }
    else
    {
        hit.let = 0.0;
    }

    auto prePoint  = step->GetPreStepPoint();
    auto postPoint = step->GetPostStepPoint();

    // Energy classification (raw, unweighted -- see comment above)
    if (pname == "proton")
    {   
        //G4cout
        //    << "PROTON STEP | "
        //    << "edep = " << edep/keV << " keV | "
        //    << "pre z = " << prePoint->GetPosition().z()/um << " um | "
        //    << "post z = " << postPoint->GetPosition().z()/um << " um | "
        //    << "volume = " << volume->GetName()
        //    << G4endl;

        fEventAction->AddProtonEdep(edep);
    }
    else if (pname == "e-")
    {
        fEventAction->AddElectronEdep(edep);
    }
    else
    {
     // recoil nuclei, ions, fragments
        fEventAction->AddRecoilEdep(edep);
    }

    // Position/time for the hit record
    G4ThreeVector pos =
        0.5 * (prePoint->GetPosition() + postPoint->GetPosition());

    hit.pos  = pos;
    hit.time = prePoint->GetGlobalTime();
    hit.process =
        postPoint->GetProcessDefinedStep()
            ? postPoint->GetProcessDefinedStep()->GetProcessName()
            : G4String("unknown");

    // Store hit
    fEventAction->AddHit(hit);

    // Collected-charge model.
    // This is ALWAYS evaluated so that both the raw deposited charge
    // (scored via edep in EventAction) and the drift/trapping-corrected
    // collected charge are available for every event. The detector's
    // "useCollectionModel" flag only decides which of the two is used
    // as the criterion for an upset in EventAction — it no longer
    // gates whether the collected charge gets computed at all.

    // Generated charge (Coulombs). Pair-creation energy depends on the
    // sensitive volume's material -- see
    // DetectorConstruction::GetSensitivePairCreationEnergy().
    double Qgen = (edep / fDetector->GetSensitivePairCreationEnergy()) * CLHEP::eplus;

    // Get sensitive layer thickness from detector
    double d_sens = fDetector->GetSensitiveThickness();

    // Sensitive volume centered at z=0, electrode assumed at +d_sens/2
    G4ThreeVector localPos =
        prePoint->GetTouchableHandle()
            ->GetHistory()
            ->GetTopTransform()
             .TransformPoint(pos);

    double z = localPos.z();

    double d = (d_sens / 2.0) - z;

    if (d < 0.0)
        d = 0.0;

    double tau_r = fDetector->GetCarrierLifetime();
    double E     = fDetector->GetElectricField();

    // Drift time: tau_d = d / (mu * E)
    // Mobility/saturation velocity depend on the sensitive volume's
    // material -- see DetectorConstruction::GetSensitiveElectronMobility()
    // / GetSensitiveSaturationVelocity().
    double mu = fDetector->GetSensitiveElectronMobility();

    double v = mu * E;

    double vsat = fDetector->GetSensitiveSaturationVelocity();

    if (v > vsat)
        v = vsat;

    double tau_d = d / v;

    // Collection efficiency: T = exp(-tau_d / tau_r)
    double T = std::exp(-tau_d / tau_r);

    // Collected charge
    double Qeff = T * Qgen;

    // Accumulate (raw, unweighted -- see comment above)
    fEventAction->AddCollectedCharge(Qeff);
}
