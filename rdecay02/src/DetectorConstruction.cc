//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4Box.hh"


#include "G4Colour.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : fVisAttributes()
{
  fTargetLength      = 0.1*cm; 
  fTargetRadius      = 0.5*cm;
  fTargetRadiusActive = 2 * mm;
  fDetectorLength    = 2.14*cm; 
  fDetectorRadius = 2.57*cm;
  fWorldLength = 0.3 * m;
  fSideAirThickness = 1. * mm;
  fDistanceFromGeToWindow1 = 5.22 * mm;
  fWindowThickness = 0.6 * mm;
  // TODO - check for 2, 5, 10, 15 and 25 cm
  fdistanceFromTarToDet = 2. * cm;
  DefineMaterials();
    
  fDetectorMessenger = new DetectorMessenger();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
    delete fDetectorMessenger;
    for (auto visAttributes : fVisAttributes) {
        delete visAttributes;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // build materials
  //
  /*fDetectorMater = 
  new G4Material("Germanium", 32, 72.61*g/mole, 5.323*g/cm3);*/
  auto nistManager = G4NistManager::Instance();
  // Ge
  G4double z = 32.0;                        // Atomic number of Germanium
  G4double a = 72.63 * g / mole;              // Atomic mass of Germanium
  G4double densityGe = 5.323 * g / cm3;         // Density of Germanium
  G4Element* Ge = new G4Element("Germanium", "Ge", z, a);

  // Create the material Germanium
  G4double temperature = 77 * kelvin;
  G4double pressure = 1. * atmosphere;
  fDetectorMater = new G4Material("Germanium", densityGe, 1, kStateSolid, temperature, pressure);
  fDetectorMater->AddElement(Ge, 1);
  //fDetectorMater = nistManager->FindOrBuildMaterial("G4_Ge");
  // Al
  Al = nistManager->FindOrBuildMaterial("G4_Al");
  
  //Carbon 
  G4double carbon_z = 6.0;                     // Atomic number of Carbon
  G4double carbon_a = 12.011 * g / mole;       // Atomic mass of Carbon
  G4double carbon_density = 2.267 * g / cm3;   // Density of Carbon in solid form
  G4Element* elC = new G4Element("Carbon", "C", carbon_z, carbon_a);
  carbon = new G4Material("CustomCarbon", carbon_density, 1, kStateSolid, temperature, pressure);
  carbon->AddElement(elC, 1);

  G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",   "O", 8, 16.00*g/mole);
  // Air
  G4int ncomponents; G4double fractionmass;      
  G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2,
                      kStateGas, temperature, pressure);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);
  fWorldMater = Air20;

  // Na22
  G4double density_air = 0.97 * g / cm3;
  G4double vacuum_pressure = 1.0e-9 * atmosphere;
  fTargetMater = new G4Material("Na22Material", 11, 22.0 * g / mole, density_air,kStateSolid,temperature, vacuum_pressure);
 
  // Tungsten
  G4double tungsten_density = 19.25 * g / cm3;
  G4double a_tungsten = 183.84 * g / mole;
  G4double z_tungsten = 74;
  G4Element* Tungsten = new G4Element("Tungsten", "W", z_tungsten, a_tungsten);
  fTungstenMater = new G4Material("SolidTungsten", tungsten_density, 1, kStateSolid, temperature, pressure);
  fTungstenMater->AddElement(Tungsten, 1);

  //Kapton 
  G4double kaptonDensity = 1.42 * g / cm3;
  fKaptonMater = new G4Material("Kapton", kaptonDensity, 3);
  fKaptonMater->AddElement(nistManager->FindOrBuildElement("C"), 5);
  fKaptonMater->AddElement(nistManager->FindOrBuildElement("H"), 10);
  fKaptonMater->AddElement(nistManager->FindOrBuildElement("O"), 5);

  // Other Options for Disk Material

  G4Material* matPE = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Material* matPP = nistManager->FindOrBuildMaterial("G4_POLYPROPYLENE");
  G4Material* matBe = nistManager->FindOrBuildMaterial("G4_Be");
  G4Material* matLi = nistManager->FindOrBuildMaterial("G4_Li");
  G4Material* matGraphite = nistManager->FindOrBuildMaterial("G4_GRAPHITE");
  G4Material* matPMMA = nistManager->FindOrBuildMaterial("G4_PLEXIGLASS"); // acrylic
  G4Material* matPTFE = nistManager->FindOrBuildMaterial("G4_TEFLON");
  G4Material* matAl = nistManager->FindOrBuildMaterial("G4_Al");
  G4Material* matSi = nistManager->FindOrBuildMaterial("G4_Si");

  // TODO - CHANGE THE DISK MATERIAL HERE
  const G4String diskChoice = "Kapton";

  if (diskChoice == "Kapton")            fDiskMater = fKaptonMater;
  else if (diskChoice == "G4_POLYETHYLENE")   fDiskMater = matPE;
  else if (diskChoice == "G4_POLYPROPYLENE")  fDiskMater = matPP;
  else if (diskChoice == "G4_Be")             fDiskMater = matBe;
  else if (diskChoice == "G4_Li")             fDiskMater = matLi;
  else if (diskChoice == "G4_GRAPHITE")       fDiskMater = matGraphite;
  else if (diskChoice == "G4_PLEXIGLASS")     fDiskMater = matPMMA;
  else if (diskChoice == "G4_TEFLON")         fDiskMater = matPTFE;
  else if (diskChoice == "G4_Al")             fDiskMater = matAl;
  else if (diskChoice == "G4_Si")             fDiskMater = matSi;
  else {
      G4Exception("DetectorConstruction::DefineMaterials()",
          "BadDiskMaterial", JustWarning,
          ("Unknown diskChoice: " + diskChoice + " — falling back to Kapton").c_str());
      fDiskMater = fKaptonMater;
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // World
  BuildWorld();
  BuildField();
  // Target
  BuildSource();
  
// Detectors
  BuildDetector(true); 
  BuildDetector(false);

  // Tungsten Cones
  BuildTungstenCones();

  // Kapton Disks
  BuildKaptonDisks();

// -------------------------------
  PrintParameters();
  
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::BuildWorld() {
    auto worldSolid
        = new G4Box("worldBox",             // name
            fWorldLength, fWorldLength, fWorldLength);   // dimentions
    fworldLogical
        = new G4LogicalVolume(worldSolid, //shape
            fWorldMater,    //material
            "fworldLogical");    // name
    fPhysiWorld
        = new G4PVPlacement(0,  // no rotation
            G4ThreeVector(),    // at (0,0,0)
            fworldLogical,       //logical volume
            "worldPhysical",    //name
            0,                  //mother volume
            false,              //no boolean operation
            0);

    G4double transparency = 0.5;

    auto worldvisAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, transparency));
    worldvisAttributes->SetVisibility(true);
    fworldLogical->SetVisAttributes(worldvisAttributes);
    fVisAttributes.push_back(worldvisAttributes);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::BuildSource() {
    G4Tubs*
        sTarget = new G4Tubs("Target",                                   //name
            0., fTargetRadiusActive, 0.5 * fTargetLength, 0., twopi); //dimensions


    fLogicTarget = new G4LogicalVolume(sTarget,           //shape
        fTargetMater,              //material
        "Target");                 //name

    new G4PVPlacement(0,                         //no rotation
        G4ThreeVector(),             //at (0,0,0)
        fLogicTarget,                //logical volume
        "Target",                    //name
        fworldLogical,                      //mother  volume
        false,                       //no boolean operation
        0);                          //copy number
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::BuildDetector(G4bool isPositiveZ)
{
    G4double zPos = (isPositiveZ ? +1 : -1) * fdistanceFromTarToDet;

    auto carbon_tube = new G4Tubs("CarbonTube", fDetectorRadius + fSideAirThickness,
        fDetectorRadius + fSideAirThickness + fWindowThickness,
        0.5 * (fDetectorLength + fWindowThickness + fDistanceFromGeToWindow1),
        0., twopi);

    auto lcarbon_tube = new G4LogicalVolume(carbon_tube, carbon, "lCarbonTube");
    new G4PVPlacement(0, G4ThreeVector(0, 0, zPos), lcarbon_tube, "lCarbonTube", fworldLogical, false, 0);

    auto air_tube = new G4Tubs("AirTube", fDetectorRadius,
        fDetectorRadius + fSideAirThickness,
        0.5 * (fDetectorLength + fDistanceFromGeToWindow1),
        0., twopi);

    auto lair_tube = new G4LogicalVolume(air_tube, fWorldMater, "lAirTube");
    new G4PVPlacement(0, G4ThreeVector(0, 0, zPos), lair_tube, "lAirTube", fworldLogical, false, 0);

    G4Tubs* sDetector = new G4Tubs("Detector",
        0., fDetectorRadius, 0.5 * fDetectorLength, 0., twopi);

    const char* detLVname = isPositiveZ ? "Detector1" : "Detector2";
    const char* detPVname = isPositiveZ ? "Detector1_phys" : "Detector2_phys";

    auto logicDetector = new G4LogicalVolume(sDetector, fDetectorMater, detLVname);
    new G4PVPlacement(0, G4ThreeVector(0, 0, zPos),
        logicDetector, detPVname, fworldLogical, false, 0);

    // Store pointers
    if (isPositiveZ) fLogicDetector1 = logicDetector;
    else             fLogicDetector2 = logicDetector;

    // Visualization
    G4double transparency = 0.5;
    auto visAttr = new G4VisAttributes(G4Colour(0.8888, 0.0, 0.0, transparency));
    logicDetector->SetVisAttributes(visAttr);
    fVisAttributes.push_back(visAttr);

    auto carbonColor = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8, transparency));
    lcarbon_tube->SetVisAttributes(carbonColor);
    fVisAttributes.push_back(carbonColor);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::BuildTungstenCones() {
    // ---------------- geometry knobs (fractions of your existing sizes) -------------
    const G4double kWallFrac = 0.9;   // inner radius = kWallFrac * outer (thicker wall -> smaller value)
    const G4double kGapDetFrac = 0.10;   // clearance to detector face as fraction of fDetectorLength
    const G4double kGapTarFrac = 0.00;   // clearance to target face as fraction of fTargetLength

	const G4double kDiskThickFrac = 1.00;   // TODOD change here to 0 or 1 - wihtout or with tungeten disk
    const G4double kDiskOuterFrac = 1.00;   // outer radius relative to fTargetRadius
    const G4double kDiskPosFrac = 0.25;   // position between target face and detector face (0..1)

    // ------------------------- derived distances ------------------------------------
    const G4double gapDet = kGapDetFrac * fDetectorLength; // near detector face (both sides)
    const G4double gapTar = kGapTarFrac * fTargetLength;   // near target face  (both sides)

    // Free space from target face (+z) to detector front face (+z)
    G4double freeSpan = fdistanceFromTarToDet - 0.5 * fDetectorLength - 0.5 * fTargetLength;
    if (freeSpan < 0.) freeSpan = 0.;

    // Cones occupy remaining space after gaps; each cone gets half
    G4double usable = freeSpan - (gapDet + gapTar);
    if (usable < 0.) usable = 0.;
    G4double h = 0.5 * usable;               // half-length for G4Cons
    if (h <= 0.) h = 0.1 * mm;               // minimal non-zero length

    // Centers of the cones (distance from world origin)
    const G4double zCenter = 0.5 * fTargetLength + gapTar + h;

    // ------------------------------ cone radii --------------------------------------
    const G4double rmax1 = fTargetRadius;                  // outer at target side
    const G4double rmin1 = kWallFrac * rmax1;              // inner (wall)
    const G4double rmax2 = fDetectorRadius;                // outer at detector side
    const G4double rmin2 = kWallFrac * rmax2;              // inner (wall)

    // ------------------------------ build +Z cone -----------------------------------
    // For G4Cons: (rmin1,rmax1) are at local -Z face; (rmin2,rmax2) at local +Z face.
    // Placed at +z without rotation, the local -Z face is nearer the target, which is correct.
    auto cone1 = new G4Cons("Cone1", rmin1, rmax1, rmin2, rmax2, h, 0., twopi);
    fLogicCone1 = new G4LogicalVolume(cone1, fTungstenMater, "LogicCone1");
    new G4PVPlacement(
        nullptr, G4ThreeVector(0, 0, +zCenter),
        fLogicCone1, "PhysCone1", fworldLogical, false, 0);

    // ------------------------------ build -Z cone (ROTATED 180°) --------------------
    // At -z, we must flip the cone so its target-side (small end defined at local -Z face)
    // also faces the target (toward +Z in global for this placement).
    auto cone2 = new G4Cons("Cone2", rmin1, rmax1, rmin2, rmax2, h, 0., twopi);
    fLogicCone2 = new G4LogicalVolume(cone2, fTungstenMater, "LogicCone2");

    auto rot180 = new G4RotationMatrix();
    rot180->rotateX(180. * deg); // flip so target-side is toward the origin

    new G4PVPlacement(
        rot180, G4ThreeVector(0, 0, -zCenter),
        fLogicCone2, "PhysCone2", fworldLogical, false, 0);

    // visuals for cones
    {
        G4double alpha = 0.5;
        auto visCone = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, alpha));
        fLogicCone1->SetVisAttributes(visCone);
        fLogicCone2->SetVisAttributes(visCone);
        fVisAttributes.push_back(visCone);
    }

    // ------------------- optional flat tungsten disks (symmetric) -------------------
	const G4double diskThickness = kDiskThickFrac * 0.05 * mm; // TODO change here the Tungsten disk thickness
    if (diskThickness > 0.) {
        const G4double halfDisk = 0.5 * diskThickness;

        // place between target face and detector face, at a fraction of that span
        const G4double zRing = 0.75 * fTargetLength;   // Kapton disk position
        const G4double smallGap = 0.1 * mm;            // tiny separation
        const G4double zDisk = zRing + smallGap;

        const G4double rDiskOuter = kDiskOuterFrac * fTargetRadius;
        const G4double rDiskInner = 0.0; // solid disk

        // +Z disk
        auto sDisk1 = new G4Tubs("TungstenDisk1", rDiskInner, rDiskOuter, halfDisk, 0., twopi);
        fLogicTungstenTube1 = new G4LogicalVolume(sDisk1, fTungstenMater, "lTungstenTube1");
        new G4PVPlacement(
            nullptr, G4ThreeVector(0, 0, +zDisk),
            fLogicTungstenTube1, "lTungstenTube1", fworldLogical, false, 0);

        // -Z disk
        auto sDisk2 = new G4Tubs("TungstenDisk2", rDiskInner, rDiskOuter, halfDisk, 0., twopi);
        fLogicTungstenTube2 = new G4LogicalVolume(sDisk2, fTungstenMater, "lTungstenTube2");
        new G4PVPlacement(
            nullptr, G4ThreeVector(0, 0, -zDisk),
            fLogicTungstenTube2, "lTungstenTube2", fworldLogical, false, 0);

        // visuals for disks
        auto visDisk = new G4VisAttributes(G4Colour(0.2, 0.2, 0.6, 0.6));
        fLogicTungstenTube1->SetVisAttributes(visDisk);
        fLogicTungstenTube2->SetVisAttributes(visDisk);
        fVisAttributes.push_back(visDisk);
    }
    else {
        fLogicTungstenTube1 = nullptr;
        fLogicTungstenTube2 = nullptr;
    }
}


void DetectorConstruction::BuildKaptonDisks() {
    // ---- geometry controls ----
    const G4bool   kAsAnnulus = false;         // ring (true) vs solid disk (false)
	const G4double kThickness = 0.4 * mm;    // TODO change Kapton disk thickness here
    const G4double kRingWidth = 2.0 * mm;     // radial material width outside the hole
    const G4double kMargin = 1.20;         // 20% safety margin on the LOS hole

    // ring z position (you used 0.75 * fTargetLength)
    const G4double zRing = 0.75 * fTargetLength;   // > 0 (forward). We'll place symmetric +-zRing.

    // detector geometry we already have
    const G4double Rdet = fDetectorRadius;
    const G4double zDetFace = fdistanceFromTarToDet - 0.5 * fDetectorLength; // +z front face distance

    // --- compute the hole radius that exactly shadows the detector face ---
    // scale = (distance of ring plane)/(distance of detector face)
    const G4double scale = (zDetFace > 0.) ? (zRing / zDetFace) : 0.0;
    G4double rHole = kMargin * Rdet * scale;                // clear aperture at ring plane
    rHole = std::max(0.2 * mm, rHole);                        // clamp to >=0.2 mm

    // pick an outer radius. Use at least the old target radius; larger is fine.
    const G4double rOuter = std::max(fTargetRadius, rHole + kRingWidth);

    // choose inner radius (0 for solid disk, rHole for annulus)
    const G4double rInner = kAsAnnulus ? rHole : 0.0;

    // ---- build shapes ----
    auto s1 = new G4Tubs("Kapton1", rInner, rOuter, 0.5 * kThickness, 0., twopi);
    fLogicKapton_tube1 = new G4LogicalVolume(s1, fDiskMater, "lKapton_tube1");
    new G4PVPlacement(nullptr, { 0,0,+zRing }, fLogicKapton_tube1, "lKapton_tube1", fworldLogical, false, 0);

    auto s2 = new G4Tubs("Kapton2", rInner, rOuter, 0.5 * kThickness, 0., twopi);
    fLogicKapton_tube2 = new G4LogicalVolume(s2, fDiskMater, "lKapton_tube2");
    new G4PVPlacement(nullptr, { 0,0,-zRing }, fLogicKapton_tube2, "lKapton_tube2", fworldLogical, false, 0);

    // ---- visuals ----
    auto vis = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0, 0.5));
    fLogicKapton_tube1->SetVisAttributes(vis);
    fLogicKapton_tube2->SetVisAttributes(vis);
    fVisAttributes.push_back(vis);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n Target : Length = " << G4BestUnit(fTargetLength,"Length")
         << " Radius = " << G4BestUnit(fTargetRadius,"Length")  
         << " Material = " << fTargetMater->GetName();
  G4cout << "\n Detector : Length = " << G4BestUnit(fDetectorLength,"Length")
         << " Tickness = " << G4BestUnit(fDetectorRadius,"Length")  
         << " Material = " << fDetectorMater->GetName() << G4endl;          
  G4cout << "\n" << fTargetMater << "\n" << fDetectorMater << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fTargetMater = pttoMaterial;
    if(fLogicTarget) { fLogicTarget->SetMaterial(fTargetMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetTargetMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fDetectorMater = pttoMaterial;
    if(fLogicDetector1) { fLogicDetector1->SetMaterial(fDetectorMater); }
    if (fLogicDetector2) { fLogicDetector2->SetMaterial(fDetectorMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

void DetectorConstruction::BuildField() {
	// TODO - change if to use magenetic field
    const G4bool   kIsMagneticField = false;
    if (!kIsMagneticField) {
        return;
	}
    auto magField = new G4UniformMagField(G4ThreeVector(0.35 * tesla, 0.35 * tesla, 0.)); // tweak 0.1–0.3 T
    auto fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField); // simple default chord finder
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetRadius(G4double value)
{
  fTargetRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetLength(G4double value)
{
  fTargetLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorThickness(G4double value)
{
  fDetectorRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorLength(G4double value)
{
  fDetectorLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetLength() const
{
  return fTargetLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetRadius() const
{
  return fTargetRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetTargetMaterial() const
{
  return fTargetMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicTarget() const
{
  return fLogicTarget;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorLength() const
{
  return fDetectorLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorThickness() const
{
  return fDetectorRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetDetectorMaterial() const
{
  return fDetectorMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicDetector1() const
{
  return fLogicDetector1;
}

G4LogicalVolume* DetectorConstruction::GetLogicDetector2() const
{
    return fLogicDetector2;
}

G4LogicalVolume* DetectorConstruction::GetLogicDisk1() const
{
    return fLogicTungstenTube1;
}

G4LogicalVolume* DetectorConstruction::GetLogicDisk2() const
{
    return fLogicTungstenTube2;
}

G4LogicalVolume* DetectorConstruction::GetLogicKapton1() const {
    return fLogicKapton_tube1;
}
G4LogicalVolume* DetectorConstruction::GetLogicKapton2() const {
    return fLogicKapton_tube2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
