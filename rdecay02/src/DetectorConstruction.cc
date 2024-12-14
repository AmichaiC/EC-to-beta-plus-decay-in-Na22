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
  fDetectorLength    = 2.14*cm; 
  fDetectorRadius = 2.57*cm;
  fWorldLength = 0.3 * m;
  fSideAirThickness = 1. * mm;
  fDistanceFromGeToWindow1 = 5.22 * mm;
  fWindowThickness = 0.6 * mm;
  DefineMaterials();
    
  fDetectorMessenger = new DetectorMessenger(this);
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
  //
  G4int ncomponents; G4double fractionmass;      
  G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2,
                      kStateGas, temperature, pressure);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);
  //
  fWorldMater = Air20;

  G4double density_air = 0.97 * g / cm3;
  fTargetMater = new G4Material("Na22Material", 11, 22.0 * g / mole, density_air,kStateSolid,temperature, pressure);
 
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
  auto worldSolid
      = new G4Box("worldBox",             // name
          fWorldLength, fWorldLength, fWorldLength);   // dimentions
  auto worldLogical
      = new G4LogicalVolume(worldSolid, //shape
          fWorldMater,    //material
          "worldLogical");    // name
  fPhysiWorld
      = new G4PVPlacement(0,  // no rotation
          G4ThreeVector(),    // at (0,0,0)
          worldLogical,       //logical volume
          "worldPhysical",    //name
          0,                  //mother volume
          false,              //no boolean operation
          0);
                            
  // Target
  //
  G4Tubs* 
  sTarget = new G4Tubs("Target",                                   //name
                  0., fTargetRadius, 0.5*fTargetLength, 0.,twopi); //dimensions


  fLogicTarget = new G4LogicalVolume(sTarget,           //shape
                             fTargetMater,              //material
                             "Target");                 //name
                               
           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(),             //at (0,0,0)
                           fLogicTarget,                //logical volume
                           "Target",                    //name
                            worldLogical,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number
    // Detector 1
    G4double distanceFromTarToDet = 5. * cm;
    auto carbon_tube1
        = new G4Tubs("CarbonTube1", fDetectorRadius+ fSideAirThickness,
            fDetectorRadius + fSideAirThickness + fWindowThickness
            , 0.5*(fDetectorLength+fWindowThickness+fDistanceFromGeToWindow1)
            , 0., twopi);

    auto lcarbon_tube1
        = new G4LogicalVolume(carbon_tube1, carbon, "lCarbonTube1");
    new G4PVPlacement(0, G4ThreeVector(0, 0, distanceFromTarToDet), lcarbon_tube1, "lCarbonTube1", worldLogical, false, 0);

    auto air_tube1
        = new G4Tubs("AirTube1", fDetectorRadius,
            fDetectorRadius + fSideAirThickness
            , 0.5 * (fDetectorLength + fDistanceFromGeToWindow1)
            , 0., twopi);

    auto lair_tube1
        = new G4LogicalVolume(air_tube1, fWorldMater, "lairTube1");
    new G4PVPlacement(0, G4ThreeVector(0, 0, distanceFromTarToDet), lair_tube1, "lairTube1", worldLogical, false, 0);

  G4Tubs* sDetector1 = new G4Tubs("Detector1",  
      0., fDetectorRadius, 0.5*fDetectorLength, 0.,twopi);


  fLogicDetector1 = new G4LogicalVolume(sDetector1,       //shape
                             fDetectorMater,            //material
                             "Detector1");               //name
                               
           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(0,0,distanceFromTarToDet),     // location    
                           fLogicDetector1,              //logical volume
                           "Detector1",                  //name
                           worldLogical,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

    // Detector 2
           auto carbon_tube2
               = new G4Tubs("CarbonTube2", fDetectorRadius + fSideAirThickness,
                   fDetectorRadius + fSideAirThickness + fWindowThickness
                   , 0.5 * (fDetectorLength + fWindowThickness + fDistanceFromGeToWindow1)
                   , 0., twopi);

           auto lcarbon_tube2
               = new G4LogicalVolume(carbon_tube1, carbon, "lCarbonTube2");
           new G4PVPlacement(0, G4ThreeVector(0, 0, -distanceFromTarToDet), lcarbon_tube2, "lCarbonTube2", worldLogical, false, 0);

           auto air_tube2
               = new G4Tubs("AirTube2", fDetectorRadius,
                   fDetectorRadius + fSideAirThickness
                   , 0.5 * (fDetectorLength + fDistanceFromGeToWindow1)
                   , 0., twopi);

           auto lair_tube2
               = new G4LogicalVolume(air_tube2, fWorldMater, "lairTube2");
           new G4PVPlacement(0, G4ThreeVector(0, 0, -distanceFromTarToDet), lair_tube2, "lairTube2", worldLogical, false, 0);

    G4Tubs* sDetector2 = new G4Tubs("Detector2",
        0., fDetectorRadius, 0.5 * fDetectorLength, 0., twopi);


    fLogicDetector2 = new G4LogicalVolume(sDetector2,       //shape
        fDetectorMater,            //material
        "Detector2");               //name

    new G4PVPlacement(0,                         //no rotation
        G4ThreeVector(0, 0, -distanceFromTarToDet),     // location    
        fLogicDetector2,              //logical volume
        "Detector2",                  //name
        worldLogical,                      //mother  volume
        false,                       //no boolean operation
        0);                          //copy number

 


 // visualization attributes ------------------------------------------------
           G4double transparency = 0.5;

    auto worldvisAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, transparency));
    worldvisAttributes->SetVisibility(true);
    worldLogical->SetVisAttributes(worldvisAttributes);
    fVisAttributes.push_back(worldvisAttributes);

    auto det1visAttributes = new G4VisAttributes(G4Colour(0.8888, 0.0, 0.0, transparency));
    fLogicDetector1->SetVisAttributes(det1visAttributes);
    fVisAttributes.push_back(det1visAttributes);

    auto det2visAttributes = new G4VisAttributes(G4Colour(0.8888, 0.0, 0.0, transparency));
    fLogicDetector2->SetVisAttributes(det2visAttributes);
    fVisAttributes.push_back(det2visAttributes);

    auto carbonColor = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8, transparency)); // silver
    lcarbon_tube1->SetVisAttributes(carbonColor);
    fVisAttributes.push_back(carbonColor);

    lcarbon_tube2->SetVisAttributes(carbonColor);
    fVisAttributes.push_back(carbonColor);

// -------------------------------
  PrintParameters();
  
  //always return the root volume
  //
  return fPhysiWorld;
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

G4double DetectorConstruction::GetTargetLength()
{
  return fTargetLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetRadius()
{
  return fTargetRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetTargetMaterial()
{
  return fTargetMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicTarget()
{
  return fLogicTarget;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorLength()
{
  return fDetectorLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorThickness()
{
  return fDetectorRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetDetectorMaterial()
{
  return fDetectorMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicDetector1()
{
  return fLogicDetector1;
}

G4LogicalVolume* DetectorConstruction::GetLogicDetector2()
{
    return fLogicDetector2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
