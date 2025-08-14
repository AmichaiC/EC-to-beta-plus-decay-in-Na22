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
/// \file Run.hh
/// \brief Definition of the Run class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include <map>

class DetectorConstruction;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run();
   ~Run() override = default;

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);         
    void CountProcesses(const G4VProcess* process, const G4String& volumeName);
    void ParticleCount(G4String, G4double, G4int); 
    void AddEdep (G4double edep1, G4double edep2, G4double edep3);
    void Run::PrintProcessFrequency();

    void IncrementPositronReachingDetector1() { fPositronReachingDetector1++; }
    void IncrementPositronReachingDetector2() { fPositronReachingDetector2++; }
    G4int GetPositronReachingDetector1() const { return fPositronReachingDetector1; }
    G4int GetPositronReachingDetector2() const { return fPositronReachingDetector2; }
    void IncrementUniqueAnnihilationCount() { fUniqueAnnihilationCount++; }
    G4int GetUniqueAnnihilationCount() const { return fUniqueAnnihilationCount; }
    void IncrementAnnihilationByVolumeFromSource(G4String);
                          
    void Merge(const G4Run*) override;
    void EndOfRun();
    void WriteActivity(G4int); 

    // Detector 1
    void AddEnterDisk1() { ++fEnterDisk1; }
    void AddReachDet1() { ++fReachDet1; }
    void AddEnterKapton1() { ++fEnterKapton1; }
    G4int GetEnterDisk1() const { return fEnterDisk1; }
    G4int GetReachDet1()  const { return fReachDet1; }
    G4int GetEnterKapton1() const { return fEnterKapton1; }


    // Detector 2
    void AddEnterDisk2() { ++fEnterDisk2; }
    void AddReachDet2() { ++fReachDet2; }
    void AddEnterKapton2() { ++fEnterKapton2; }
    G4int GetEnterKapton2() const { return fEnterKapton2; }
    G4int GetEnterDisk2() const { return fEnterDisk2; }
    G4int GetReachDet2()  const { return fReachDet2; }

    void AddEnterGamma(const G4String& volName, G4bool is511, G4bool is1274);
    void AddExitGamma(const G4String& volName, G4bool is511, G4bool is1274);


   
  private:
    struct ParticleData {
     ParticleData()
       : fCount(0), fEmean(0.), fEmin(0.), fEmax(0.) {}
     ParticleData(G4int count, G4double ekin, G4double emin, G4double emax)
       : fCount(count), fEmean(ekin), fEmin(emin), fEmax(emax) {}
     G4int     fCount;
     G4double  fEmean;
     G4double  fEmin;
     G4double  fEmax;
    };
     
  private:
    DetectorConstruction* fDetector = nullptr;
    G4ParticleDefinition* fParticle = nullptr;
    G4double              fEkin = 0.;
    
    G4double fEdepTarget = 0., fEdepTarget2 = 0.;
    G4double fEdepDetect11 = 0., fEdepDetect12 = 0., fEdepDetect21 = 0., fEdepDetect22 = 0.;

    G4int fPositronReachingDetector1 = 0;
    G4int fPositronReachingDetector2 = 0;
    G4int fUniqueAnnihilationCount = 0;

    // Tungsten disk counters
    // Detector 1 (positive Z)
    G4int fEnterDisk1 = 0;
    G4int fReachDet1 = 0;
    // Detector 2 (negative Z)
    G4int fEnterDisk2 = 0;
    G4int fReachDet2 = 0;

    // Kapton disk counters
    G4int fEnterKapton1 = 0, fEnterKapton2 = 0;

    std::map<G4String,G4int>        fProcCounter1;
    std::map<G4String,G4int>        fProcCounter2;   
    std::map<G4String, G4int>        fProcCounter3;
    std::map<G4String,ParticleData> fParticleDataMap1;                    
    std::map<G4String,ParticleData> fParticleDataMap2;
    std::map<G4String, ParticleData> fParticleDataMap3;

    std::map<G4String, std::map<G4String, G4int>> fProcCounterByVolume;

    // maps the number of annihilations in each volume by the positrons emitted by the source
    std::map<G4String, G4int> fAnnihilationByVolumeFromSource; 

    struct SimpleFlow { long long in511 = 0, out511 = 0, in1274 = 0, out1274 = 0; };
    std::map<G4String, SimpleFlow> fSimpleFlow;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

