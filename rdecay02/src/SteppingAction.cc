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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "MyTrackInfo.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction(EventAction* event)
    :fEventAction(event)
{
    fDetector = const_cast<DetectorConstruction*>(
        static_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    // Get the current run and analysis manager
    Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    auto analysisManager = G4AnalysisManager::Instance();

    // Pre- and post-step volumes
    auto preStepVolume = aStep->GetPreStepPoint()->GetTouchableHandle()
        ->GetVolume()->GetLogicalVolume();
    auto postStepVolume = aStep->GetPostStepPoint()->GetTouchableHandle()
        ->GetVolume()->GetLogicalVolume();

    // Volume ID for energy deposition
    G4int iVol = 0;
    if (preStepVolume == fDetector->GetLogicTarget())         iVol = 1;
    else if (preStepVolume == fDetector->GetLogicDetector1()) iVol = 2;
    else if (preStepVolume == fDetector->GetLogicDetector2()) iVol = 3;

    // Count all processes by volume
    auto process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
    G4String volName = preStepVolume->GetName();
    run->CountProcesses(process, volName);

    // Energy deposition and ntuple filling
    G4double edep = aStep->GetTotalEnergyDeposit();
    if (edep > 0.) {
        G4double time = aStep->GetPreStepPoint()->GetGlobalTime();
        G4double weight = aStep->GetPreStepPoint()->GetWeight();
        fEventAction->AddEdep(iVol, edep, time, weight);

        G4int id = 2;
        analysisManager->FillNtupleDColumn(id, 0, edep);
        analysisManager->FillNtupleDColumn(id, 1, time / s);
        analysisManager->FillNtupleDColumn(id, 2, weight);
        analysisManager->AddNtupleRow(id);
    }

    // Particle and track
    auto track = aStep->GetTrack();
    auto particle = track->GetParticleDefinition();
    G4String name = particle->GetParticleName();

    // Retrieve user info with birth volume
    auto trackInfo = dynamic_cast<MyTrackInfo*>(track->GetUserInformation());
    if (!trackInfo) return;
    auto birthVolume = trackInfo->GetBirthVolume();
    if (!birthVolume) return;

    // Only count annihilations of e+ born in the Target
    if (name == "e+"
        && birthVolume == fDetector->GetLogicTarget()
        && process
        && process->GetProcessName() == "annihil")
    {
        if (!trackInfo->IsAnnihilationCounted()) {
            run->IncrementUniqueAnnihilationCount();
            run->IncrementAnnihilationByVolumeFromSource(volName);
            trackInfo->SetAnnihilationCounted(true);
        }
    }

    // code for if 1274 keV photon is absorbs in the disks or getting through - doesn't work, program collapse

    // only 1274 keV gammas
    //if (name == "gamma") {
    //    /*HandleGammaStep(aStep, run, preStepVolume, postStepVolume, birthVolume);*/
    //    G4double E = track->GetKineticEnergy();

    //    // Only 1274 keV ±10 keV born in target?
    //    if (std::abs(E - 1.274 * MeV) >= 10 * keV ||
    //        birthVolume != fDetector->GetLogicTarget())
    //    {
    //        return;
    //    }

    //    // Only count when crossing a boundary
    //    if (preStepVolume == postStepVolume) return;

    //    // Forward path
    //    if ((preStepVolume == fDetector->GetLogicKapton1() || postStepVolume == fDetector->GetLogicKapton1()) &&
    //        !trackInfo->IsCountedKapton1())
    //    {
    //        run->AddEnterKapton1();
    //        trackInfo->SetCountedKapton1();
    //    }

    //    if ((preStepVolume == fDetector->GetLogicDisk1() || postStepVolume == fDetector->GetLogicDisk1()) &&
    //        !trackInfo->IsCountedDisk1())
    //    {
    //        run->AddEnterDisk1();
    //        trackInfo->SetCountedDisk1();
    //    }

    //    if (postStepVolume == fDetector->GetLogicDetector1() &&
    //        !trackInfo->IsCountedDet1())
    //    {
    //        run->AddReachDet1();
    //        trackInfo->SetCountedDet1();
    //    }

    //    // Backward path
    //    if ((preStepVolume == fDetector->GetLogicKapton2() || postStepVolume == fDetector->GetLogicKapton2()) &&
    //        !trackInfo->IsCountedKapton2())
    //    {
    //        run->AddEnterKapton2();
    //        trackInfo->SetCountedKapton2();
    //    }

    //    if ((preStepVolume == fDetector->GetLogicDisk2() || postStepVolume == fDetector->GetLogicDisk2()) &&
    //        !trackInfo->IsCountedDisk2())
    //    {
    //        run->AddEnterDisk2();
    //        trackInfo->SetCountedDisk2();
    //    }

    //    if (postStepVolume == fDetector->GetLogicDetector2() &&
    //        !trackInfo->IsCountedDet2())
    //    {
    //        run->AddReachDet2();
    //        trackInfo->SetCountedDet2();
    //    }
    //}

}


void SteppingAction::HandleGammaStep(const G4Step* aStep,
    Run* run,
    G4LogicalVolume* preVol,
    G4LogicalVolume* postVol,
    G4LogicalVolume* birthVol) const
{
    auto track = aStep->GetTrack();
    
}
