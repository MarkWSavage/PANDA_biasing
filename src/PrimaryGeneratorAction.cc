#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "Randomize.hh"

std::atomic<G4double> PrimaryGeneratorAction::fsBeamXY{10*um};

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fParticleGun = new G4ParticleGun(1);

    // Defaults
    fEnergy = 50*MeV;
    fBeamXY = 10*um;

    fMessenger =
        new G4GenericMessenger(this, "/sim/", "Simulation control");

    fMessenger->DeclarePropertyWithUnit("energy", "MeV", fEnergy);
    fMessenger->DeclareMethodWithUnit(
        "beamXY", "um",
        &PrimaryGeneratorAction::SetBeamXY,
        "Beam spot size (square, uniform transverse profile)");
    // /sim/particle is declared by DetectorConstruction, not here --
    // see DetectorConstruction::GetParticleName() for why.

    // Beam direction toward +Z
    fParticleGun->SetParticleMomentumDirection(
        G4ThreeVector(0,0,1)
    );
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
    delete fMessenger;
}

void PrimaryGeneratorAction::SetBeamXY(G4double val)
{
    fBeamXY = val;
    fsBeamXY.store(val);
}

G4double PrimaryGeneratorAction::GetBeamXY()
{
    return fsBeamXY.load();
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // Particle name lives on DetectorConstruction (see
    // DetectorConstruction::GetParticleName() for why) since the
    // biasing operator also needs it.
    auto detector =
        static_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction()
        );

    auto particle =
        G4ParticleTable::GetParticleTable()
            ->FindParticle(detector->GetParticleName());

    fParticleGun->SetParticleDefinition(particle);

    fParticleGun->SetParticleEnergy(fEnergy);

    // Uniform beam over square area
    G4double x =
        (G4UniformRand() - 0.5)*fBeamXY;

    G4double y =
        (G4UniformRand() - 0.5)*fBeamXY;

    // Always start well upstream in vacuum. This is only actually true
    // because DetectorConstruction::Construct() truncates
    // SurroundingVolume's front face flush with the dead layer's own
    // front face (see fActualSurroundingFrontZ there) -- before that
    // fix, SurroundingVolume (built from the sensitive material, not
    // vacuum) extended past the dead layer's front face for any
    // dead+sensitive stack thinner than its ~10-20% auto-grow margin,
    // silently costing every primary a fixed chunk of energy crossing
    // solid material this comment claimed was vacuum (e.g. ~140 keV for
    // a 5.486 MeV alpha, confirmed against the ORTEC BU-014-050-100
    // datasheet's independently-quoted dead-layer loss).
    G4double sourceZ =
        -(detector->GetSensitiveThickness()/2.0
           + detector->GetDeadThickness()
           + 1.0*um);

    fParticleGun->SetParticlePosition(
        G4ThreeVector(x,y,sourceZ)
    );

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
