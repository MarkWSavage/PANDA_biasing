#include "DetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4UserLimits.hh"
#include "G4RunManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "SEEBiasingOperator.hh"

#include <algorithm>
#include <cmath>
#include <vector>

DetectorConstruction::DetectorConstruction()
{
    fSensitiveLogical = nullptr;

    fSensitiveXY = 10*um;
    fSensitiveThickness = 2*um;
    fDeadThickness = 5*um;
    fStepSize = 0.01*um;

    fMessenger =
        new G4GenericMessenger(this, "/sim/", "Simulation control");

auto& xyCmd =
    fMessenger->DeclarePropertyWithUnit(
        "sensitiveXY",
        "um",
        fSensitiveXY);
xyCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& deadXYCmd =
    fMessenger->DeclarePropertyWithUnit(
        "deadXY",
        "um",
        fDeadXY,
        "Lateral (X/Y) size of the dead layer. Independent of "
        "sensitiveXY -- defaults to the same value (10 um) but does "
        "NOT track it, so set both explicitly if they should match.");
deadXYCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& particleCmd =
    fMessenger->DeclareProperty(
        "particle",
        fParticleName,
        "Primary particle name (e.g. proton, neutron, alpha, deuteron, "
        "triton, He3, GenericIon, or one of the heavy-ion primaries: "
        "C12, F19, Cl35, Ni58, I127, Au197). Also selects which species' "
        "hadronic inelastic process the SEEBiasingOperator biases -- see "
        "ConstructSDandField() and PANDA.cc.");
particleCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& sensMatCmd =
    fMessenger->DeclareProperty(
        "sensitiveMaterial",
        fSensitiveMaterialName,
        "Sensitive volume material: Si, Ge, GaAs, SiC, or GaN. Also "
        "selects the pair-creation energy and mobility/saturation-"
        "velocity constants SteppingAction uses for the MeV->charge "
        "conversion -- see DetectorConstruction::GetSensitivePairCreation"
        "Energy() and friends for values/sources.");
sensMatCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& deadMatCmd =
    fMessenger->DeclareProperty(
        "deadMaterial",
        fDeadMaterialName,
        "Dead-layer/electrode material: Si, Ge, GaAs, SiC, GaN, SiO2, "
        "Al2O3, TiO2, Au, W, Pb, or Ta. Independent of sensitiveMaterial. "
        "Au/W/Pb/Ta are for modeling metallization/via-plug materials "
        "directly against the sensitive volume (e.g. gold contact "
        "layers, tungsten plugs -- tungsten plugs specifically have "
        "been reported to undergo neutron-induced fission, a distinct "
        "mechanism from proton-induced reactions in Si).");
deadMatCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& thickCmd =
    fMessenger->DeclarePropertyWithUnit(
        "sensitiveThickness",
        "um",
        fSensitiveThickness);
thickCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& angleCmd =
    fMessenger->DeclarePropertyWithUnit(
        "incidentAngle",
        "deg",
        fIncidentAngle,
        "Angle of the primary beam from the sensitive volume's normal "
        "(0 = normal incidence, default). Approximates a tilted beam by "
        "dividing the sensitive volume's effective thickness by "
        "cos(incidentAngle) -- the standard chord-length-elongation "
        "model used for tilt-angle corrections in heavy-ion SEE testing "
        "-- instead of actually rotating the beam or geometry. Valid "
        "range [0, 90) deg. The lateral (XY) footprint is NOT adjusted "
        "-- see Documentation/PANDA_MASTER_DESIGN.");
angleCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& collectCmd =
    fMessenger->DeclareProperty(
        "useCollectionModel",
        fUseCollectionModel
    );
collectCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& lifetimeCmd =
    fMessenger->DeclarePropertyWithUnit(
        "carrierLifetime",
        "ns",
        fCarrierLifetime
    );
lifetimeCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& fieldCmd =
    fMessenger->DeclarePropertyWithUnit(
        "electricField",
        "volt/m",
        fElectricField
    );

fieldCmd.SetStates(G4State_PreInit, G4State_Idle);


auto& deadCmd =
    fMessenger->DeclarePropertyWithUnit(
        "deadThickness",
        "um",
        fDeadThickness);
deadCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& stepCmd =
    fMessenger->DeclarePropertyWithUnit(
        "stepSize",
        "um",
        fStepSize
    );
stepCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& surrXYCmd =
    fMessenger->DeclarePropertyWithUnit(
        "surroundingXY",
        "um",
        fSurroundingXY,
        "Lateral (X/Y) size of the bulk material surrounding the dead+"
        "sensitive stack, matching the sensitive volume's material. "
        "Auto-grown at construction time if the dead/sensitive stack "
        "itself is larger. Default 100 um -- large enough to capture "
        "nearby nuclear reactions contributing recoils into the "
        "sensitive volume, without the overhead of containing a "
        "primary particle's full range (typically cm-scale, but "
        "reactions that far away essentially never reach the "
        "sensitive volume anyway).");
surrXYCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& surrThickCmd =
    fMessenger->DeclarePropertyWithUnit(
        "surroundingThickness",
        "um",
        fSurroundingThickness,
        "Thickness of the surrounding bulk material -- see "
        "surroundingXY.");
surrThickCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& biasCmd =
    fMessenger->DeclareProperty(
        "biasCrossSectionFactor",
        fBiasCrossSectionFactor,
        "Multiplier on the hadronic inelastic cross section for "
        "whichever particle /sim/particle selects, in the sensitive+"
        "dead volumes. 1.0 = no bias (default). Set to match CREME-MC's "
        "Hadronic Cross Section Multiplier for direct comparability. "
        "ALWAYS verify a run with this left at 1.0 reproduces the "
        "unbiased baseline before trusting a boosted run."
    );
biasCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& secondaryNeutronBiasCmd =
    fMessenger->DeclareProperty(
        "secondaryNeutronBiasFactor",
        fSecondaryNeutronBiasFactor,
        "Multiplier on neutronInelastic for SECONDARY neutrons (not "
        "the run's primary particle) in the sensitive+dead+surrounding "
        "volumes. 1.0 = no secondary-neutron bias (default). The "
        "primary-species biasing operator only ever biases tracks "
        "matching /sim/particle -- a proton primary's secondary "
        "neutrons pass through completely unbiased otherwise, and "
        "with silicon's cm-scale neutron mean free path, a second "
        "reaction near the sensitive volume from one is exceedingly "
        "rare to sample without this. Ignored if /sim/particle is "
        "itself neutron (already covered by the primary operator)."
    );
secondaryNeutronBiasCmd.SetStates(G4State_PreInit, G4State_Idle);

}

DetectorConstruction::~DetectorConstruction()
{
    delete fMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    auto nist = G4NistManager::Instance();

    auto vacuum =
        nist->FindOrBuildMaterial("G4_Galactic");

    auto sensitiveMaterial = ResolveMaterial(fSensitiveMaterialName);
    auto deadMaterial      = ResolveMaterial(fDeadMaterialName);

    // Effective (tilt-adjusted) sensitive thickness -- see
    // fEffectiveSensitiveThickness in the header for why this is a
    // fresh recompute from fSensitiveThickness every call, not an
    // in-place update. 0 deg (the default) gives cos(0)=1, exactly
    // reproducing the untilted baseline for every macro that doesn't
    // set /sim/incidentAngle.
    if (fIncidentAngle < 0.0*deg || fIncidentAngle >= 90.0*deg)
    {
        G4Exception(
            "DetectorConstruction::Construct()",
            "InvalidIncidentAngle",
            FatalException,
            "/sim/incidentAngle must be in [0, 90) deg -- cos(theta) "
            "must stay positive and finite for the thickness/cos(theta) "
            "approximation to make sense."
        );
    }

    fEffectiveSensitiveThickness =
        fSensitiveThickness / std::cos(fIncidentAngle);

    // Same elongation applied to the dead layer -- see
    // fEffectiveDeadThickness in the header for why this is needed too
    // (a tilted beam traverses more physical dead-layer material before
    // reaching the sensitive volume, not just more sensitive-volume
    // material).
    fEffectiveDeadThickness =
        fDeadThickness / std::cos(fIncidentAngle);

    // The chord-length-elongation model only stretches thickness, never
    // the lateral (XY) footprint -- see fIncidentAngle in the header.
    // At extreme angles this can make the sensitive volume's elongated
    // thickness far exceed its own (un-widened) XY footprint, turning
    // it into an artificially tall, narrow column. Real multiple-
    // Coulomb/nuclear scattering accumulated over the whole dead+
    // sensitive path can then push a primary out the side of that
    // column before it finishes traversing the sensitive volume, so it
    // never deposits the chord-length-elongated energy the model
    // intends there -- biasing deposited-energy results low. Confirmed
    // empirically (2026-07-19, 1um sensitive/10um dead/10um XY stack):
    // fraction of events with zero sensitive-volume deposit stayed
    // under 1% through 45 deg, jumped to 5.7% at 89 deg (sensitive
    // elongation ~57x = 5.7x its own 10um XY), and 91% at 89.9 deg
    // (~573x). Checked against the DEAD layer's own ratio too -- it
    // does NOT correlate as cleanly (e.g. 80 deg already exceeds 5x the
    // dead layer's XY but showed no elevated zero-fraction at all,
    // 0.26%, actually below the 0-deg baseline) -- what apparently
    // matters is staying in-bounds specifically during the SCORED
    // (sensitive-volume) traversal, not the dead layer's own individual
    // ratio, so only the sensitive volume is checked here. Not
    // something this model fixes (would require also widening the
    // lateral footprint at high angle, a separate, larger lift) -- a
    // real device at a grazing angle extreme enough to trigger this
    // would have to traverse neighboring sensitive volumes, isolation
    // oxide, and packaging first anyway, so it isn't a physically
    // realizable single-junction test condition. Just a heads-up: a
    // simple, empirically-calibrated aspect-ratio heuristic (elongated
    // sensitive thickness > 5x its own lateral footprint), not a
    // precise scattering-length calculation.
    if (fEffectiveSensitiveThickness > 5.0 * fSensitiveXY)
    {
        G4Exception(
            "DetectorConstruction::Construct()",
            "IncidentAngleLateralEscape",
            JustWarning,
            "incidentAngle has elongated the sensitive volume's "
            "thickness to more than 5x its own lateral (XY) footprint. "
            "Real multiple-Coulomb/nuclear scattering over that "
            "artificially long, narrow path can push primaries out the "
            "(un-widened) side before they finish traversing the "
            "sensitive volume, biasing deposited-energy results low -- "
            "see Documentation/PANDA_MASTER_DESIGN's incident-angle "
            "section. Treat results at this angle with caution, or "
            "reduce incidentAngle / widen sensitiveXY."
        );
    }

    // Surrounding volume matches the sensitive volume's material (bulk
    // substrate the junction is fabricated in). Grow it automatically
    // if the dead/sensitive stack itself is larger than the requested
    // surrounding size -- otherwise a large-sensitive-volume preset
    // (e.g. 5000 um) wouldn't fit inside a 1 mm default surrounding
    // box, which Geant4 would reject as an invalid mother/daughter
    // overlap.
    G4double totalThickness =
        fEffectiveDeadThickness + fEffectiveSensitiveThickness;

    G4double surroundingXY =
        std::max(fSurroundingXY, 1.2 * std::max(fSensitiveXY, fDeadXY));

    G4double surroundingThickness =
        std::max(fSurroundingThickness, 1.2 * totalThickness);

    // SurroundingVolume's front (beam-entrance) face is truncated flush
    // with the dead layer's own front face -- NOT centered symmetrically
    // about the stack's midpoint like the back/lateral margins still
    // are. A symmetric box left ~10-20% of surroundingThickness worth of
    // solid (sensitive-material) silicon sitting *upstream* of the dead
    // layer, in what every other part of this codebase (PrimaryGenerator
    // Action's own comment included) assumed was vacuum -- confirmed via
    // the ORTEC BU-014-050-100 alpha-spectroscopy datasheet cross-check:
    // a 5.486 MeV alpha through that unintended margin lost ~140 keV it
    // had no business losing before ever reaching the modeled dead
    // layer. The back/lateral margins are left exactly as before (still
    // sized from surroundingThickness/surroundingXY above) since those
    // are what the McNulty nearby-nuclear-reaction validation depends
    // on; only the front face moves.
    G4double deadFrontZ =
        -(fEffectiveSensitiveThickness/2.0 + fEffectiveDeadThickness);

    G4double oldSymmetricBackFaceZ = surroundingThickness/2.0;

    G4double surroundingThicknessZ = oldSymmetricBackFaceZ - deadFrontZ;
    G4double surroundingCenterZ =
        (deadFrontZ + oldSymmetricBackFaceZ)/2.0;

    // Thin vacuum margin around the surrounding volume (Geant4's
    // outermost placed volume convention). World stays centered at the
    // origin, so its half-extent must cover whichever face (front or
    // back) sits further from it -- no longer symmetric now that
    // SurroundingVolume itself is offset.
    G4double worldXY = 1.2 * surroundingXY;
    G4double worldZ =
        2.4 * std::max(std::abs(deadFrontZ), std::abs(oldSymmetricBackFaceZ));

    auto solidWorld =
        new G4Box(
            "World",
            worldXY/2,
            worldXY/2,
            worldZ/2
        );

    auto logicWorld =
        new G4LogicalVolume(
            solidWorld,
            vacuum,
            "World"
        );

    auto physWorld =
        new G4PVPlacement(
            nullptr,
            G4ThreeVector(),
            logicWorld,
            "World",
            nullptr,
            false,
            0
        );

    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Bulk material surrounding the dead+sensitive stack -- see
    // fSurroundingXY/fSurroundingThickness in the header for why this
    // exists (nearby nuclear reactions contributing recoils into the
    // sensitive volume; PANDA under-predicts high-deposited-energy
    // events without it -- see Known Limitations in
    // Documentation/PANDA_MASTER_DESIGN).
    auto solidSurrounding =
        new G4Box(
            "SurroundingVolume",
            surroundingXY/2,
            surroundingXY/2,
            surroundingThicknessZ/2
        );

    auto logicSurrounding =
        new G4LogicalVolume(
            solidSurrounding,
            sensitiveMaterial,
            "SurroundingVolume"
        );

    logicSurrounding->SetVisAttributes(G4VisAttributes::GetInvisible());

    fSurroundingLogical = logicSurrounding;

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, 0, surroundingCenterZ),
        logicSurrounding,
        "SurroundingVolume",
        logicWorld,
        false,
        0
    );

    // Dead layer directly below sensitive volume
    auto solidDead =
        new G4Box(
            "DeadLayer",
            fDeadXY/2,
            fDeadXY/2,
            fEffectiveDeadThickness/2
        );

    auto logicDead =
        new G4LogicalVolume(
            solidDead,
            deadMaterial,
            "DeadLayer"
        );

    fDeadLogical = logicDead;

     auto deadVis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // blue
         deadVis->SetForceSolid(true);
         logicDead->SetVisAttributes(deadVis);

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(
            0,
            0,
            -(fEffectiveSensitiveThickness/2 + fEffectiveDeadThickness/2)
                - surroundingCenterZ
        ),
        logicDead,
        "DeadLayer",
        logicSurrounding,
        false,
        0
    );

    // Sensitive volume
    auto solidSensitive =
        new G4Box(
            "SensitiveVolume",
            fSensitiveXY/2,
            fSensitiveXY/2,
            fEffectiveSensitiveThickness/2
        );


    fSensitiveLogical =
        new G4LogicalVolume(
            solidSensitive,
            sensitiveMaterial,
            "SensitiveVolume"
        );

    auto sensVis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red
        sensVis->SetForceSolid(true);
        fSensitiveLogical->SetVisAttributes(sensVis);

    fSensitiveLogical->SetUserLimits(
        //new G4UserLimits(0.01*um)     //Appropriate for nano-electronics
        //new G4UserLimits(0.50*um)     //Appropriate for large structures
        new G4UserLimits(fStepSize)
    );

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, 0, -surroundingCenterZ),
        fSensitiveLogical,
        "SensitiveVolume",
        logicSurrounding,
        false,
        0
    );

    // Production cuts scaled to the sensitive+dead stack, not left at
    // Geant4's stock ~0.7-1 mm default. That default is fine down to the
    // 1-14 um devices this has actually been run at (a delta ray's real
    // range is already smaller than the device), but becomes a real gap
    // at deep-submicron scale (e.g. an 80 nm-node junction's ~20-100 nm
    // depletion depth): every secondary's energy gets dumped locally via
    // restricted stepping instead of being transported as an explicit
    // track, silently losing whatever delta-ray/charge-sharing range
    // structure is comparable to or larger than the sensitive volume
    // itself.
    //
    // Applied via a dedicated G4Region covering just the sensitive+dead
    // volumes (not the whole geometry) -- Geant4 only uses the finer cut
    // for secondaries actually produced inside that region, so the
    // (often much larger, auto-grown up to several mm -- see
    // fSurroundingXY/Thickness) surrounding volume and world keep the
    // default cut and the normal tracking cost. Scaling this cut globally
    // instead would force fine-grained tracking through that entire
    // surrounding volume too, which is exactly the "1 mm shell made runs
    // impractically slow" problem noted below, but far worse.
    //
    // Cut length = smallest relevant dimension / 10 (a commonly-used
    // rule of thumb for resolving secondary-particle structure within a
    // region -- not an exact science), clamped to [1 nm, 0.7 mm]: the
    // floor matches SteppingAction's LET boundary-artifact cutoff (going
    // finer buys nothing -- Geant4's own lowest tracking energy threshold
    // is reached well before 1 nm in Si anyway), the ceiling means this
    // never coarsens relative to Geant4's stock default for large devices
    // (e.g. the 50-1000 um geometries in Macros/run_ceiling_*.mac) --
    // it only ever sharpens the resolution as the device shrinks.
    G4double smallestDim =
        std::min(
            {fSensitiveXY, fEffectiveSensitiveThickness, fDeadXY, fEffectiveDeadThickness}
        );

    G4double cutLength = smallestDim / 10.0;
    cutLength = std::clamp(cutLength, 1.0 * nm, 0.7 * mm);

    G4cout << "DetectorConstruction::Construct() -- SensitiveRegion "
              "production cut: "
           << cutLength / nm << " nm (smallest dimension "
           << smallestDim / nm << " nm / 10)" << G4endl;

    // Every PANDA macro calls /run/initialize then /run/reinitializeGeometry
    // (see any Macros/*.mac), so Construct() runs twice per invocation --
    // once from initialize, once from the explicit reinit. Region names
    // must be globally unique in G4RegionStore, so unconditionally `new`-ing
    // a region here on the second call collides with the first call's
    // still-registered region (a fatal region-store corruption, not just a
    // harmless warning despite what the G4Exception text says). Reuse the
    // existing region object if Construct() has already run once instead.
    auto sensitiveRegion =
        G4RegionStore::GetInstance()->GetRegion("SensitiveRegion", false);
    if (!sensitiveRegion)
        sensitiveRegion = new G4Region("SensitiveRegion");

    sensitiveRegion->AddRootLogicalVolume(fSensitiveLogical);
    sensitiveRegion->AddRootLogicalVolume(logicDead);

    auto sensitiveCuts = new G4ProductionCuts();
    sensitiveCuts->SetProductionCut(cutLength);
    sensitiveRegion->SetProductionCuts(sensitiveCuts);


    return physWorld;
}

G4LogicalVolume* DetectorConstruction::GetSensitiveLogical() const
{
    return fSensitiveLogical;
}

G4Material* DetectorConstruction::ResolveMaterial(const G4String& name)
{
    auto nist = G4NistManager::Instance();

    if (name == "Si")    return nist->FindOrBuildMaterial("G4_Si");
    if (name == "Ge")    return nist->FindOrBuildMaterial("G4_Ge");
    if (name == "GaAs")  return nist->FindOrBuildMaterial("G4_GALLIUM_ARSENIDE");
    if (name == "SiO2")  return nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    if (name == "Al2O3") return nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
    if (name == "TiO2")  return nist->FindOrBuildMaterial("G4_TITANIUM_DIOXIDE");
    if (name == "Au")    return nist->FindOrBuildMaterial("G4_Au");
    if (name == "W")     return nist->FindOrBuildMaterial("G4_W");
    if (name == "Pb")    return nist->FindOrBuildMaterial("G4_Pb");
    if (name == "Ta")    return nist->FindOrBuildMaterial("G4_Ta");

    // SiC and GaN aren't in Geant4's NIST compound database -- build
    // them from elements. G4NistManager::FindOrBuildMaterial caches by
    // name internally, but "SiC"/"GaN" aren't NIST names it recognizes,
    // so we cache these ourselves to avoid rebuilding (and re-warning
    // Geant4's material table about a duplicate name) on every call.
    if (name == "SiC")
    {
        static G4Material* siliconCarbide = nullptr;
        if (!siliconCarbide)
        {
            siliconCarbide = new G4Material("SiC", 3.21*g/cm3, 2);
            siliconCarbide->AddElement(nist->FindOrBuildElement("Si"), 1);
            siliconCarbide->AddElement(nist->FindOrBuildElement("C"), 1);
        }
        return siliconCarbide;
    }

    if (name == "GaN")
    {
        static G4Material* galliumNitride = nullptr;
        if (!galliumNitride)
        {
            galliumNitride = new G4Material("GaN", 6.15*g/cm3, 2);
            galliumNitride->AddElement(nist->FindOrBuildElement("Ga"), 1);
            galliumNitride->AddElement(nist->FindOrBuildElement("N"), 1);
        }
        return galliumNitride;
    }

    G4Exception(
        "DetectorConstruction::ResolveMaterial()",
        "InvalidMaterial",
        FatalException,
        ("Unknown material name: '" + name +
         "' -- expected one of Si, Ge, GaAs, SiC, GaN, SiO2, Al2O3, "
         "TiO2, Au, W, Pb, Ta.").c_str()
    );
    return nullptr;
}

G4double DetectorConstruction::GetSensitivePairCreationEnergy() const
{
    // Mean energy to create one electron-hole pair. Typical room-
    // temperature bulk literature values -- treat as approximate,
    // same spirit as the "(tune these)" mobility/vsat comment below.
    if (fSensitiveMaterialName == "Si")   return 3.6*eV;
    // Standard cited Ge value: 2.9 eV.
    if (fSensitiveMaterialName == "Ge")   return 2.9*eV;
    if (fSensitiveMaterialName == "GaAs") return 4.2*eV;
    // Direct 4H-SiC p-n diode measurement: 7.83 +/- 0.02 eV (more
    // precise than the older ~7.6 eV figure this used to use).
    if (fSensitiveMaterialName == "SiC")  return 7.83*eV;
    // Direct GaN radiation-detector measurement: 10 +/- 0.5 eV
    // (more authoritative than the older ~8.9 eV figure this used to
    // use -- see commit history).
    if (fSensitiveMaterialName == "GaN")  return 10.0*eV;

    G4Exception(
        "DetectorConstruction::GetSensitivePairCreationEnergy()",
        "InvalidMaterial",
        FatalException,
        ("Unknown /sim/sensitiveMaterial name: '" + fSensitiveMaterialName +
         "' -- expected one of Si, Ge, GaAs, SiC, GaN.").c_str()
    );
    return 3.6*eV;
}

G4double DetectorConstruction::GetSensitiveElectronMobility() const
{
    // Low-field electron mobility. Typical room-temperature bulk
    // literature values (tune these).
    if (fSensitiveMaterialName == "Si")   return 1350.0 * cm2/volt/second;
    // Real Ge detectors are cooled to ~77 K (room-temp thermal leakage
    // current is prohibitive), where mobility is far higher than this
    // 300 K figure (~40,000-70,000 cm^2/V/s vs 3900). Kept at the
    // room-temperature value for consistency with every other material
    // here (none of which model temperature dependence) -- it doesn't
    // change this simulation's result anyway, since mobility*field
    // already vastly exceeds vsat at either value, so drift velocity
    // is saturation-capped regardless.
    if (fSensitiveMaterialName == "Ge")   return 3900.0 * cm2/volt/second;
    if (fSensitiveMaterialName == "GaAs") return 8500.0 * cm2/volt/second;
    if (fSensitiveMaterialName == "SiC")  return 1000.0 * cm2/volt/second;
    if (fSensitiveMaterialName == "GaN")  return 1300.0 * cm2/volt/second;

    G4Exception(
        "DetectorConstruction::GetSensitiveElectronMobility()",
        "InvalidMaterial",
        FatalException,
        ("Unknown /sim/sensitiveMaterial name: '" + fSensitiveMaterialName +
         "' -- expected one of Si, Ge, GaAs, SiC, GaN.").c_str()
    );
    return 1350.0 * cm2/volt/second;
}

G4double DetectorConstruction::GetSensitiveSaturationVelocity() const
{
    // Typical room-temperature bulk literature values (tune these).
    // GaAs's electron velocity-field curve is non-monotonic (Gunn
    // effect / negative differential mobility): it peaks around
    // 2.0e7 cm/s near ~3-4 kV/cm, then droops before re-settling at
    // high field. This simple mu*E-capped-at-vsat model can't
    // reproduce that NDR dip, so 2.0e7 (the commonly tabulated
    // "saturation velocity" figure, e.g. Sze & Ng) is used as the cap
    // -- more realistic than reusing Si's 1.0e7, which was an
    // oversimplification from the first cut of this model.
    if (fSensitiveMaterialName == "Si")   return 1.0e7 * cm/second;
    // Cited 300 K electron saturation velocity: 0.7e7 cm/s (more
    // precise than the older 6.0e6 figure this used to use).
    if (fSensitiveMaterialName == "Ge")   return 7.0e6 * cm/second;
    if (fSensitiveMaterialName == "GaAs") return 2.0e7 * cm/second;
    if (fSensitiveMaterialName == "SiC")  return 2.0e7 * cm/second;
    // GaN, like GaAs, has a non-monotonic velocity-field curve
    // (Gunn effect): peaks around 2.6e7 cm/s near ~145 kV/cm for
    // undoped GaN, then droops at higher field. Same modeling
    // choice as GaAs: use the peak/commonly-tabulated figure as
    // the cap, since this simple model can't reproduce the NDR dip.
    if (fSensitiveMaterialName == "GaN")  return 2.6e7 * cm/second;

    G4Exception(
        "DetectorConstruction::GetSensitiveSaturationVelocity()",
        "InvalidMaterial",
        FatalException,
        ("Unknown /sim/sensitiveMaterial name: '" + fSensitiveMaterialName +
         "' -- expected one of Si, Ge, GaAs, SiC, GaN.").c_str()
    );
    return 1.0e7 * cm/second;
}

void DetectorConstruction::ConstructSDandField()
{
    // Attach the SEE biasing operator to the sensitive and dead-layer
    // volumes. This is called once (serial) or once per worker thread
    // (MT), after Construct() -- the Geant4-recommended, thread-safe
    // place to attach biasing operators (each thread gets its own
    // operator instance, matching Geant4's official biasing examples).
    //
    // With fBiasCrossSectionFactor == 1.0 (the default), the operator
    // still runs but proposes no actual change to the physical cross
    // section, so results should be statistically identical to not
    // having biasing installed at all. This is the required sanity
    // check before trusting any boosted-factor run -- see the
    // messenger command help text and SEEBiasingOperator.hh for why.
    //
    // The particle to bias is whatever /sim/particle selected (read
    // via fParticleName, set before /run/initialize triggers this
    // method). PANDA.cc must have wrapped that same particle's
    // hadronic inelastic process via PhysicsBias() at physics-
    // construction time -- if it hasn't (e.g. a particle name not in
    // that list), SEEBiasingOperator::StartRun() prints an explicit
    // WARNING and biasing is silently inert for this run.
    // Heavy-ion primaries: unlike proton/neutron/alpha/deuteron/triton/
    // He3 (built into Geant4's particle table from startup), a specific
    // heavy ion doesn't exist as a G4ParticleDefinition until explicitly
    // created via G4IonTable::GetIon(Z,A,...) -- and that call requires
    // GenericIon's own process manager to already be built (it clones
    // that onto the new ion). ConstructSDandField() runs too early for
    // this: confirmed empirically (G4IonTable::CreateIon() logs PART105
    // "GenericIon is not ready" if attempted here) -- geometry/SD
    // construction (this method) happens before physics construction in
    // Geant4's own init sequence, for master and worker threads alike.
    // The earliest point that's guaranteed safe is
    // SEEBiasingOperator::StartRun(), triggered on the Idle->GeomClosed
    // state transition, which happens only after
    // G4RunManagerKernel::RunInitialization() has already built physics
    // tables (see its source -- physicsInitialized is asserted true on
    // entry, and BuildPhysicsTables() runs before the state transition
    // that fires StartRun()). So for a known heavy-ion name, resolution
    // is deferred: primaryDef stays null and (Z,A) is passed through
    // instead, resolved later by StartRun() -- see SEEBiasingOperator.cc.
    // By the time PrimaryGeneratorAction::GeneratePrimaries() runs for
    // the first event (necessarily after StartRun()), the ion already
    // exists, so its own FindParticle(fParticleName) needs no changes.
    struct IonSpec { const char* name; G4int z; G4int a; };
    static const IonSpec kHeavyIons[] = {
        {"C12",   6,  12},
        {"F19",   9,  19},
        {"Cl35",  17, 35},
        {"Ni58",  28, 58},
        {"I127",  53, 127},
        {"Au197", 79, 197},
    };

    G4ParticleDefinition* primaryDef = nullptr;
    G4int pendingIonZ = 0;
    G4int pendingIonA = 0;

    for (const auto& ion : kHeavyIons)
    {
        if (fParticleName == ion.name)
        {
            pendingIonZ = ion.z;
            pendingIonA = ion.a;
            break;
        }
    }

    if (pendingIonZ == 0)
    {
        primaryDef =
            G4ParticleTable::GetParticleTable()->FindParticle(fParticleName);

        if (!primaryDef)
        {
            G4Exception(
                "DetectorConstruction::ConstructSDandField()",
                "InvalidParticle",
                FatalException,
                ("Unknown /sim/particle name: '" + fParticleName +
                 "' -- not found in G4ParticleTable.").c_str()
            );
        }
    }

    // Geant4 allows only one G4VBiasingOperator attached per logical
    // volume, so biasing both the primary species and secondary
    // neutrons in the same volumes requires ONE operator instance
    // covering both, not two separate instances (an earlier version of
    // this code tried two instances; AttachTo() silently refused the
    // second one for each volume, logging "can not be attached ...
    // already used by another operator" and leaving secondary-neutron
    // biasing completely inert). See SEEBiasingOperator.hh.
    std::vector<SEEBiasingOperator::SpeciesBias> speciesBiases;
    speciesBiases.push_back(
        {primaryDef, fBiasCrossSectionFactor, false, pendingIonZ, pendingIonA}
    );

    // Secondary neutrons (e.g. produced by a proton primary's own
    // reactions) are NOT biased by the primary entry above -- it only
    // ever biases tracks matching fParticleName. Skip adding a second
    // entry when the primary itself is neutron: that entry already
    // biases every neutron track, primary and secondary alike, so a
    // duplicate entry for the same species/process would be redundant.
    if (fParticleName != "neutron" && fSecondaryNeutronBiasFactor != 1.0)
    {
        G4ParticleDefinition* neutronDef =
            G4ParticleTable::GetParticleTable()->FindParticle("neutron");

        speciesBiases.push_back(
            {neutronDef, fSecondaryNeutronBiasFactor, true}
        );
    }

    auto* biasingOperator = new SEEBiasingOperator(speciesBiases);

    if (fSensitiveLogical)
        biasingOperator->AttachTo(fSensitiveLogical);

    if (fDeadLogical)
        biasingOperator->AttachTo(fDeadLogical);

    // The surrounding volume needs biasing too -- otherwise the rare
    // nearby reactions it exists to capture would be just as
    // statistically hard to sample as everywhere else was before
    // biasing existed, defeating the purpose of adding it.
    if (fSurroundingLogical)
        biasingOperator->AttachTo(fSurroundingLogical);
}
