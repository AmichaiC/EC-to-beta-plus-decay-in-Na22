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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <vector>
#include "G4VisAttributes.hh"


class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction() override;

  public:
  
    G4VPhysicalVolume* Construct() override;
    
    void SetTargetLength (G4double value);
    void SetTargetRadius (G4double value);
    void SetTargetMaterial (G4String);
    
    void SetDetectorLength(G4double value);           
    void SetDetectorThickness(G4double value);  
    void SetDetectorMaterial(G4String);               
                   
    void PrintParameters();
    
  public:
      
    G4double GetTargetLength();
    G4double GetTargetRadius();
    G4Material* GetTargetMaterial();       
    G4LogicalVolume* GetLogicTarget();
    
    G4double GetDetectorLength();
    G4double GetDetectorThickness();
    G4Material* GetDetectorMaterial(); 
    G4LogicalVolume* GetLogicDetector1();
    G4LogicalVolume* GetLogicDetector2();  
    G4LogicalVolume* GetLogicDisk1();
    G4LogicalVolume* GetLogicDisk2();
                       
  private:
  
    G4double           fTargetLength = 0.; 
    G4double           fTargetRadius = 0.;
    G4double           fTargetRadiusActive = 0.;
    G4Material*        fTargetMater  = nullptr;
    G4LogicalVolume*   fLogicTarget  = nullptr;
                 
    G4double           fDetectorLength    = 0.;
    G4double           fDetectorRadius = 0.;
    G4double           fSideAirThickness = 0.;
    G4double           fDistanceFromGeToWindow1 = 0.;
    G4double           fWindowThickness = 0.;
    G4double           fdistanceFromTarToDet = 0.;
    G4Material*        fDetectorMater = nullptr;
    G4Material*        Al = nullptr;
    G4Material*         carbon = nullptr;
    G4Material* fTungstenMater = nullptr;
    G4Material* fKaptonMater = nullptr;
    G4LogicalVolume*   fLogicDetector1 = nullptr;
    G4LogicalVolume*   fLogicDetector2 = nullptr;
    G4LogicalVolume* fworldLogical = nullptr;
    G4LogicalVolume* fLogicTungstenTube1 = nullptr;
    G4LogicalVolume* fLogicTungstenTube2 = nullptr;
    G4LogicalVolume* fLogicKapton_tube1 = nullptr;
    G4LogicalVolume* fLogicKapton_tube2 = nullptr;
    
               
    G4double           fWorldLength = 0.;
    G4Material*        fWorldMater  = nullptr;     
    G4VPhysicalVolume* fPhysiWorld  = nullptr;
                
    DetectorMessenger* fDetectorMessenger = nullptr;

    std::vector<G4VisAttributes*> fVisAttributes;

  private:
    
    void               DefineMaterials();
    G4VPhysicalVolume* ConstructVolumes(); 
    void BuildWorld();
    void BuildSource();
    void BuildDetector1();
    void BuildDetector2();
    void BuildTungstenCones();
    void BuildKaptonDisks();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

