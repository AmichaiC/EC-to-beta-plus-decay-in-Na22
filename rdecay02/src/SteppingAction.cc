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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Positron.hh"
#include "MyTrackInfo.hh"
                          
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* event)
: fDetector(det), fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    // Get pre- and post-step logical volumes
    G4LogicalVolume* preStepVolume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    G4LogicalVolume* postStepVolume = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    // Determine volume ID for process counting (if needed)
    G4int iVol = 0;
    if (preStepVolume == fDetector->GetLogicTarget())      iVol = 1;
    else if (preStepVolume == fDetector->GetLogicDetector1()) iVol = 2;
    else if (preStepVolume == fDetector->GetLogicDetector2()) iVol = 3;
    // iVol = 0 if in another volume

    // Count processes regardless (as in your original code)
    const G4StepPoint* endPoint = aStep->GetPostStepPoint();
    const G4VProcess* process = endPoint->GetProcessDefinedStep();
    run->CountProcesses(process, iVol);

    // Energy deposit processing (unchanged)
    G4double edepStep = aStep->GetTotalEnergyDeposit();
    if (edepStep <= 0.) return;
    G4double time = aStep->GetPreStepPoint()->GetGlobalTime();
    G4double weight = aStep->GetPreStepPoint()->GetWeight();
    fEventAction->AddEdep(iVol, edepStep, time, weight);

    G4int id = 2;
    analysisManager->FillNtupleDColumn(id, 0, edepStep);
    analysisManager->FillNtupleDColumn(id, 1, time / s);
    analysisManager->FillNtupleDColumn(id, 2, weight);
    analysisManager->AddNtupleRow(id);

    // Part that doesn't work for counting the number of positrons that were detected
    // in the detectors and created in the target 

    //// Get the particle definition and name
    //const G4ParticleDefinition* particle = aStep->GetTrack()->GetParticleDefinition();
    //G4String particleName = particle->GetParticleName();

    //// Retrieve the creation volume from the track's touchable handle
    //G4VPhysicalVolume* creationPhysVolume = aStep->GetTrack()->GetTouchableHandle()->GetVolume(0);
    //if (!creationPhysVolume) return;
    //G4LogicalVolume* creationVolume = creationPhysVolume->GetLogicalVolume();
    //if (!creationVolume) return;

    // //Count positrons reaching a detector from the target 
    // //We want to count a positron only once when it first leaves the target and enters a detector.
    //if (particleName == "e+" && creationVolume == fDetector->GetLogicTarget()) {
    //    // Retrieve the custom user track information
    //    MyTrackInfo* trackInfo = dynamic_cast<MyTrackInfo*>(aStep->GetTrack()->GetUserInformation());
    //    if (!trackInfo) {
    //        trackInfo = new MyTrackInfo();
    //        // Use const_cast because GetTrack() returns a const pointer;
    //        // we want to set the user information once.
    //        const_cast<G4Track*>(aStep->GetTrack())->SetUserInformation(trackInfo);
    //    }
    //    // Only count if this positron hasn't been counted already.
    //    if (!trackInfo->IsCounted()) {
    //        // To avoid multiple counts during many steps, we update only if the step
    //        // is a boundary crossing (i.e. the pre-step volume is different from the post-step volume)
    //        if (preStepVolume != postStepVolume) {
    //            // Now check if the positron has entered a detector volume
    //            if (postStepVolume == fDetector->GetLogicDetector1()) {
    //                run->IncrementPositronReachingDetector1();
    //                trackInfo->SetCounted(true);
    //            }
    //            else if (postStepVolume == fDetector->GetLogicDetector2()) {
    //                run->IncrementPositronReachingDetector2();
    //                trackInfo->SetCounted(true);
    //            }
    //        }
    //    }
    //}
}


