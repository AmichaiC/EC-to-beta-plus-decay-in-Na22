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
  // TODO - checkk for 2, 5, 10, 15 and 25 cm
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
  G4Element* elC = new G4Element("Carbon", "C", z, a);
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

    auto logicDetector = new G4LogicalVolume(sDetector, fDetectorMater, "Detector");
    new G4PVPlacement(0, G4ThreeVector(0, 0, zPos), logicDetector, "Detector", fworldLogical, false, 0);

    // Store pointers
    if (isPositiveZ) {
        fLogicDetector1 = logicDetector;
    }
    else {
        fLogicDetector2 = logicDetector;
    }

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
    G4RotationMatrix* rotation180 = new G4RotationMatrix();
    rotation180->rotateX(180.0 * deg);
    G4double rmax1 = fTargetRadius;                    // Outer radius at the base
    G4double rmin1 = 0.95 * rmax1;                               // Inner radius at the base
    G4double rmax2 = fDetectorRadius;                  // Outer radius at the top
    G4double rmin2 = 0.95 * rmax2;                              // Inner radius at the top
    G4double distanceFromDetToCone = 0.2 * cm;
    G4double distanceFromTargetToCone = fTargetLength;
    G4double h = 0.5 * (fdistanceFromTarToDet - fTargetLength - 0.5 * fDetectorLength - 2 * distanceFromDetToCone);

    // Cone to Detector 1
    G4Cons* cone1 = new G4Cons("Cone1", rmin1, rmax1, rmin2, rmax2, h, 0., twopi);

    G4LogicalVolume* logicCone1 = new G4LogicalVolume(cone1, fTungstenMater, "LogicCone1");
    new G4PVPlacement(0, G4ThreeVector(0, 0, h + distanceFromTargetToCone), logicCone1, "PhysCone1", fworldLogical, false, 0);

    // Cone to Detector 2
    G4Cons* cone2 = new G4Cons("Cone2", rmin1, rmax1, rmin2, rmax2, h, 0., twopi);

    G4LogicalVolume* logicCone2 = new G4LogicalVolume(cone2, fTungstenMater, "LogicCone2");
    new G4PVPlacement(rotation180, G4ThreeVector(0, 0, -h - distanceFromTargetToCone), logicCone2, "PhysCone2", fworldLogical, false, 0);


    // TODO - CHANGE HERE - do 0.25 mm, 0.05
    G4double TungstenDiskThickness = 0.25 * mm;
    // disk for cone 1
    auto tungsten_tube1
        = new G4Tubs("TungstenTube1", 0, fTargetRadius
            , 0.5 * TungstenDiskThickness, 0., twopi);
    fLogicTungstenTube1
        = new G4LogicalVolume(tungsten_tube1, fTungstenMater, "lTungstenTube1");
    new G4PVPlacement(0, G4ThreeVector(0, 0, distanceFromTargetToCone + 0.5 * fTargetLength + 0.5 * TungstenDiskThickness), fLogicTungstenTube1, "lTungstenTube1", fworldLogical, false, 0);


    // disk for cone 2
    auto tungsten_tube2
        = new G4Tubs("TungstenTube2", 0, fTargetRadius
            , 0.5 * TungstenDiskThickness, 0., twopi);
    fLogicTungstenTube2
        = new G4LogicalVolume(tungsten_tube2, fTungstenMater, "lTungstenTube2");
    new G4PVPlacement(0, G4ThreeVector(0, 0, -distanceFromTargetToCone - 0.5 * fTargetLength - 0.5 * TungstenDiskThickness), fLogicTungstenTube2, "lTungstenTube2", fworldLogical, false, 0);

    // VisAttributes for Cones
    G4double transparency = 0.5;
    auto coneVisAttr1 = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, transparency)); // blue
    logicCone1->SetVisAttributes(coneVisAttr1);
    fLogicTungstenTube1->SetVisAttributes(coneVisAttr1);
    fVisAttributes.push_back(coneVisAttr1);

    auto coneVisAttr2 = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, transparency)); // blue
    logicCone2->SetVisAttributes(coneVisAttr2);
    fLogicTungstenTube2->SetVisAttributes(coneVisAttr2);
    fVisAttributes.push_back(coneVisAttr2);
}

void DetectorConstruction::BuildKaptonDisks() {
    G4double distanceFromTargetToDisk = 0.75 * fTargetLength;
   // todo was before 0.005 mm, 1. , 0.1
    G4double KaptonDiskThickness = 0.1 * mm;

    //Kapton Tube 1
    auto Kapton_tube1
        = new G4Tubs("KaptonTube1", 0, fTargetRadius
            , 0.5 * KaptonDiskThickness, 0., twopi);
    fLogicKapton_tube1
        = new G4LogicalVolume(Kapton_tube1, fKaptonMater, "lKapton_tube1");
    new G4PVPlacement(0, G4ThreeVector(0, 0, distanceFromTargetToDisk), fLogicKapton_tube1, "lKapton_tube1", fworldLogical, false, 0);

    //Kapton Tube 2
    auto Kapton_tube2
        = new G4Tubs("KaptonTube2", 0, fTargetRadius
            , 0.5 * KaptonDiskThickness, 0., twopi);
    fLogicKapton_tube2
        = new G4LogicalVolume(Kapton_tube2, fKaptonMater, "lKapton_tube2");
    new G4PVPlacement(0, G4ThreeVector(0, 0, -distanceFromTargetToDisk), fLogicKapton_tube2, "lKapton_tube2", fworldLogical, false, 0);


    // VisAttributes for Disks
    G4double transparency = 0.5;
    auto KaptonDiskVisAttr1 = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0, transparency)); // orange
    fLogicKapton_tube1->SetVisAttributes(KaptonDiskVisAttr1);
    fLogicKapton_tube2->SetVisAttributes(KaptonDiskVisAttr1);
    fVisAttributes.push_back(KaptonDiskVisAttr1);
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
