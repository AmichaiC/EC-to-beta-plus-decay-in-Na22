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

namespace {
    inline bool InWindow(G4double E, G4double center, G4double halfWidth) {
        return std::abs(E - center) <= halfWidth;
    }
    inline bool entersVolume(G4LogicalVolume* pre, G4LogicalVolume* post, G4LogicalVolume* V) {
        return (V && post == V && pre != V);
    }
}

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

    // Pre- and post-step volumes (guard post at boundaries)
    auto preStepVolume = aStep->GetPreStepPoint()->GetTouchableHandle()
        ->GetVolume()->GetLogicalVolume();

    G4LogicalVolume* postStepVolume = nullptr;
    {
        auto postH = aStep->GetPostStepPoint()->GetTouchableHandle();
        if (postH && postH->GetVolume()) {
            postStepVolume = postH->GetVolume()->GetLogicalVolume();
        }
    }

    // Volume ID for energy deposition
    G4int iVol = 0;
    if (preStepVolume == fDetector->GetLogicTarget())         iVol = 1;
    else if (preStepVolume == fDetector->GetLogicDetector1()) iVol = 2;
    else if (preStepVolume == fDetector->GetLogicDetector2()) iVol = 3;

    // Count all processes by volume (keep as you had)
    auto process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
    G4String volName = preStepVolume->GetName();
    run->CountProcesses(process, volName);

    // Energy deposition and ntuple filling (your original T3 writes)
    G4double edep = aStep->GetTotalEnergyDeposit();
    if (edep > 0.) {
        G4double time = aStep->GetPreStepPoint()->GetGlobalTime();
        G4double weight = aStep->GetPreStepPoint()->GetWeight();
        fEventAction->AddEdep(iVol, edep, time, weight);

        G4int id = 2; // T3
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

    // ------------------------------
    // Photon blocking / reaching counters (minimal overhead)
    // ------------------------------
    // Only proceed for gammas and for boundary crossings
    // ------------------------------
// Simple photon flow counters
// ------------------------------
    // ------------------------------
// Simple photon flow counters with per-track enter/exit pairing
// ------------------------------
    if (name == "gamma") {
        // only count at geometry boundaries
        if (aStep->GetPostStepPoint()->GetStepStatus() != fGeomBoundary) return;

        const G4double E = aStep->GetTrack()->GetKineticEnergy();
        const G4double halfWin = 10.0 * keV;
        const G4bool is511 = std::abs(E - 0.511 * MeV) <= halfWin;
        const G4bool is1274 = std::abs(E - 1.274 * MeV) <= halfWin;
        if (!(is511 || is1274)) return;

        // origin filter (Na-22 chain): born in target OR from annihilation
        bool fromNa22Chain = false;
        if (birthVolume == fDetector->GetLogicTarget()) fromNa22Chain = true;
        else {
            auto cr = aStep->GetTrack()->GetCreatorProcess();
            if (cr && cr->GetProcessName() == "annihil") fromNa22Chain = true;
        }
        if (!fromNa22Chain) return;

        // track info (we need it for inside/outside bookkeeping)
        auto trackInfo = dynamic_cast<MyTrackInfo*>(aStep->GetTrack()->GetUserInformation());
        if (!trackInfo) {
            trackInfo = new MyTrackInfo();
            const_cast<G4Track*>(aStep->GetTrack())->SetUserInformation(trackInfo);
            trackInfo->SetBirthVolume(birthVolume); // keep your existing behavior
        }

        auto run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

        auto preLV = aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume();
        G4LogicalVolume* postLV = nullptr;
        if (auto postH = aStep->GetPostStepPoint()->GetTouchableHandle(); postH && postH->GetVolume())
            postLV = postH->GetVolume()->GetLogicalVolume();

        // helper: count enter/exit only with proper pairing using MyTrackInfo
        auto countFor = [&](G4LogicalVolume* V) {
            if (!V) return;
            const G4String vname = V->GetName();

            // ENTER: post==V && pre!=V
            if (postLV == V && preLV != V) {
                run->AddEnterGamma(vname, is511, is1274);
                trackInfo->MarkEnter(V);
            }
            // EXIT: pre==V && post!=V  --> only if we had marked as inside
            if (preLV == V && postLV != V) {
                if (trackInfo->WasInside(V)) {
                    run->AddExitGamma(vname, is511, is1274);
                    trackInfo->MarkExit(V);
                }
                // else: photon was born inside V (no prior enter) -> do not count exit
            }
            };

        // Count for the parts you care about:
        countFor(fDetector->GetLogicKapton1());
        countFor(fDetector->GetLogicKapton2());
        countFor(fDetector->GetLogicDisk1());   // tungsten near +Z
        countFor(fDetector->GetLogicDisk2());   // tungsten near -Z

        // Optionally include detectors:
        // countFor(fDetector->GetLogicDetector1());
        // countFor(fDetector->GetLogicDetector2());
    }


}



void SteppingAction::HandleGammaStep(const G4Step* aStep,
    Run* run,
    G4LogicalVolume* preVol,
    G4LogicalVolume* postVol,
    G4LogicalVolume* birthVol) const
{
    auto track = aStep->GetTrack();
    
}
