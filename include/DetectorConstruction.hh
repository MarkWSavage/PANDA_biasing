#ifndef DetectorConstruction_h
#define DetectorConstruction_h
#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4GenericMessenger;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    // Constructs and attaches the SEE biasing operator to the
    // sensitive and dead-layer logical volumes. Called once (serial)
    // or once per worker thread (MT) by the run manager, after
    // Construct(). This is the Geant4-recommended, thread-safe place
    // to attach biasing operators -- see SEEBiasingOperator.hh.
    virtual void ConstructSDandField();

    G4LogicalVolume* GetSensitiveLogical() const;
    G4LogicalVolume* GetDeadLogical() const { return fDeadLogical; }

    // Effective thickness (fSensitiveThickness/cos(fIncidentAngle)),
    // recomputed by Construct() -- see fEffectiveSensitiveThickness.
    // Everything downstream (PrimaryGeneratorAction's beam start Z,
    // SteppingAction's drift-distance model) reads thickness only
    // through this getter, so both stay consistent with the tilted
    // geometry automatically.
    G4double GetSensitiveThickness() const { return fEffectiveSensitiveThickness; }
    G4double GetDeadThickness() const { return fDeadThickness; }
    G4double GetSensitiveXY() const { return fSensitiveXY; }
    G4double GetDeadXY() const { return fDeadXY; }
    G4bool GetUseCollectionModel() const { return fUseCollectionModel; }
    G4double GetCarrierLifetime() const { return fCarrierLifetime; }
    G4double GetElectricField() const { return fElectricField; }
    G4double GetBiasCrossSectionFactor() const { return fBiasCrossSectionFactor; }

    // Owns /sim/particle (moved here from PrimaryGeneratorAction): the
    // biasing operator attached in ConstructSDandField() needs to know
    // which particle to bias, and a G4UI command path can only be
    // declared by one messenger, so this is the single source of truth
    // -- PrimaryGeneratorAction reads it back via this getter.
    const G4String& GetParticleName() const { return fParticleName; }

    const G4String& GetSensitiveMaterialName() const { return fSensitiveMaterialName; }
    const G4String& GetDeadMaterialName() const { return fDeadMaterialName; }

    // Material-dependent constants used by SteppingAction's MeV->charge
    // conversion (pair-creation energy) and collection-efficiency drift
    // model (mobility/saturation velocity). Keyed off whatever
    // /sim/sensitiveMaterial selected -- see the .cc for values and
    // sources. Only the sensitive volume's material matters here:
    // SteppingAction only computes Qgen/Qeff for steps inside the
    // sensitive volume (see its early-return on GetSensitiveLogical()).
    G4double GetSensitivePairCreationEnergy() const;
    G4double GetSensitiveElectronMobility() const;
    G4double GetSensitiveSaturationVelocity() const;

private:
    // Resolves a material name (e.g. "Si", "GaAs", "SiO2") to a
    // G4Material, either via Geant4's NIST database or, for compounds
    // NIST doesn't carry (SiC, GaN), by building it from elements. Used
    // by Construct() for both the sensitive and dead volumes -- either
    // may be set to any name this recognizes.
    G4Material* ResolveMaterial(const G4String& name);
    G4LogicalVolume* fSensitiveLogical;
    G4LogicalVolume* fDeadLogical = nullptr;
    G4LogicalVolume* fSurroundingLogical = nullptr;
    G4GenericMessenger* fMessenger;

    G4double fStepSize = 0.018*um;
    G4double fSensitiveXY = 10*um;
    G4double fDeadXY = 10*um;
    G4double fSensitiveThickness = 10*um;
    G4double fDeadThickness = 5*um;
    G4String fParticleName = "proton";

    // Angle of the primary beam from the sensitive volume's normal, set
    // via /sim/incidentAngle (default 0 = normal incidence, preserving
    // existing behavior for every macro that doesn't set it). Rather
    // than actually rotating the beam or geometry -- which would also
    // need to reshape the dead layer, surrounding volume, and
    // secondary-neutron-biasing region to stay consistent, a much
    // larger lift -- this approximates the tilt by lengthening the
    // sensitive volume's effective thickness by 1/cos(theta): the same
    // chord-length-elongation model used for tilt-angle corrections in
    // heavy-ion SEE testing (LET_eff = LET(0 deg)/cos(theta)). Only
    // valid for beam-footprint-covers-device test conditions (see
    // Documentation/PANDA_MASTER_DESIGN) -- the lateral (XY) footprint
    // is deliberately left unchanged.
    G4double fIncidentAngle = 0*deg;

    // Recomputed by Construct() every call as
    // fSensitiveThickness/cos(fIncidentAngle) -- NOT written back into
    // fSensitiveThickness itself, since Construct() runs twice per
    // invocation (every macro calls /run/initialize then
    // /run/reinitializeGeometry) and overwriting the raw configured
    // value would double-apply the 1/cos(theta) scaling on the second
    // call. This is the value actually used to build the sensitive
    // volume's G4Box and everything derived from it (dead-layer
    // placement, surrounding-volume auto-grow sizing, production-cut
    // length), and the value GetSensitiveThickness() returns.
    G4double fEffectiveSensitiveThickness = 10*um;

    // Bulk material surrounding the dead+sensitive stack, matching the
    // sensitive volume's material (e.g. bulk Si around a Si junction).
    // Needed so nearby nuclear reactions in that surrounding material
    // can contribute recoils into the sensitive volume -- without it,
    // the stack sits in vacuum and PANDA under-predicts high-deposited-
    // energy events (see the McNulty et al. comparison discrepancy in
    // Documentation/PANDA_MASTER_DESIGN's Known Limitations). 100 um is
    // a practical default: large enough to capture reactions close
    // enough to matter, without the prohibitive overhead of containing
    // a primary's full multi-cm range (a 1 mm surrounding volume made
    // runs impractically slow -- most of that extra volume is too far
    // from the sensitive volume to ever contribute a reaching recoil
    // anyway). Construct() grows these automatically if the sensitive/
    // dead stack itself is larger (e.g. the 5000 um sensitiveXY preset).
    G4double fSurroundingXY = 100*um;
    G4double fSurroundingThickness = 100*um;

    // Defaults preserve pre-existing behavior (both volumes were
    // hardcoded to silicon) for any macro that doesn't set these.
    G4String fSensitiveMaterialName = "Si";
    G4String fDeadMaterialName = "Si";

    G4bool   fUseCollectionModel = true;
    G4double fCarrierLifetime    = 10.0 * ns;
    G4double fElectricField      = 1.0 * kilovolt/cm;

    // Multiplier applied to the hadronic inelastic cross section for
    // whichever particle fParticleName selects, in the sensitive + dead
    // volumes (see SEEBiasingOperator). 1.0 = no bias (default); set
    // via /sim/biasCrossSectionFactor. ALWAYS verify a factor=1.0 run
    // reproduces the unbiased baseline before trusting results from a
    // boosted factor.
    G4double fBiasCrossSectionFactor = 1.0;

    // Multiplier applied to neutronInelastic for SECONDARY neutrons
    // (i.e. not the run's primary particle) in the sensitive + dead +
    // surrounding volumes. 1.0 = no secondary-neutron bias (default,
    // preserves existing behavior for every macro that doesn't set
    // this). The primary-species operator above only ever biases
    // tracks matching fParticleName -- a proton primary's secondary
    // neutrons pass through completely unbiased otherwise, and with
    // silicon's cm-scale neutron mean free path, a second reaction
    // near the sensitive volume from one is exceedingly rare to sample
    // without this. Set via /sim/secondaryNeutronBiasFactor. Ignored
    // (no second operator created) when fParticleName is itself
    // "neutron", since the primary-species operator already biases
    // every neutron track, primary and secondary alike, in that case.
    G4double fSecondaryNeutronBiasFactor = 1.0;

};

#endif
