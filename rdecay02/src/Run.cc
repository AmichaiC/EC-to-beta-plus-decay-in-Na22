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
/// \file Run.cc
/// \brief Implementation of the Run class
//
// 

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4ProcessTable.hh"
#include "G4Radioactivation.hh"
#include "G4TwoVector.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "MyTrackInfo.hh"
#include "G4RunManager.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run()
{ 
   fDetector = const_cast<DetectorConstruction*>(
        static_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(const G4VProcess* process, const G4String& volumeName)
{
  if (!process) return;
    G4String procName = process->GetProcessName();
    // Update the nested map: if volumeName not present, it will be default-constructed.
    fProcCounterByVolume[volumeName][procName]++;
}


void Run::IncrementAnnihilationByVolumeFromSource(G4String volumeName)
{
    // Increment the count of annihilations for positrons
    // created in the target, in the given volume.
    fAnnihilationByVolumeFromSource[volumeName]++;
}
                  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleCount(G4String name, G4double Ekin, G4int iVol)
{
  if (iVol == 1) {
   std::map<G4String, ParticleData>::iterator it = fParticleDataMap1.find(name);
   if ( it == fParticleDataMap1.end()) {
     fParticleDataMap1[name] = ParticleData(1, Ekin, Ekin, Ekin);
   }
   else {
     ParticleData& data = it->second;
     data.fCount++;
     data.fEmean += Ekin;
     //update min max
     G4double emin = data.fEmin;
     if (Ekin < emin) data.fEmin = Ekin;
     G4double emax = data.fEmax;
     if (Ekin > emax) data.fEmax = Ekin; 
   }   
  }
  
  if (iVol == 2) {
   std::map<G4String, ParticleData>::iterator it = fParticleDataMap2.find(name);
   if ( it == fParticleDataMap2.end()) {
     fParticleDataMap2[name] = ParticleData(1, Ekin, Ekin, Ekin);
   }
   else {
     ParticleData& data = it->second;
     data.fCount++;
     data.fEmean += Ekin;
     //update min max
     G4double emin = data.fEmin;
     if (Ekin < emin) data.fEmin = Ekin;
     G4double emax = data.fEmax;
     if (Ekin > emax) data.fEmax = Ekin; 
   }   
  } 

  if (iVol == 3) {
      std::map<G4String, ParticleData>::iterator it = fParticleDataMap3.find(name);
      if (it == fParticleDataMap3.end()) {
          fParticleDataMap3[name] = ParticleData(1, Ekin, Ekin, Ekin);
      }
      else {
          ParticleData& data = it->second;
          data.fCount++;
          data.fEmean += Ekin;
          //update min max
          G4double emin = data.fEmin;
          if (Ekin < emin) data.fEmin = Ekin;
          G4double emax = data.fEmax;
          if (Ekin > emax) data.fEmax = Ekin;
      }
  }
}
                 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// TODO CHANGED edep3 back to edep 2 in the method.
void Run::AddEdep(G4double edep1, G4double edep2, G4double edep3)
{ 
  fEdepTarget  += edep1;
  fEdepTarget2 += edep1*edep1;
  fEdepDetect11  += edep2;
  fEdepDetect12 += edep2*edep2; 
  fEdepDetect21 += edep3;
  fEdepDetect22 += edep3 * edep3;
}
                 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  for (const auto& kv : static_cast<const Run*>(run)->fSimpleFlow) {
      auto& dst = fSimpleFlow[kv.first];
      dst.in511 += kv.second.in511;
      dst.out511 += kv.second.out511;
      dst.in1274 += kv.second.in1274;
      dst.out1274 += kv.second.out1274;
  }
  //primary particle info
  //
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;
  
  // accumulate sums
  //
  fEdepTarget   += localRun->fEdepTarget;  
  fEdepTarget2  += localRun->fEdepTarget2;
  fEdepDetect11   += localRun->fEdepDetect11;  
  fEdepDetect12  += localRun->fEdepDetect12;  
  fEdepDetect21 += localRun->fEdepDetect21;
  fEdepDetect22 += localRun->fEdepDetect22;
      
  //map: processes count in target
  
  std::map<G4String,G4int>::const_iterator itp1;
  for ( itp1 = localRun->fProcCounter1.begin();
        itp1 != localRun->fProcCounter1.end(); ++itp1 ) {

    G4String procName = itp1->first;
    G4int localCount = itp1->second;
    if ( fProcCounter1.find(procName) == fProcCounter1.end()) {
      fProcCounter1[procName] = localCount;
    }
    else {
      fProcCounter1[procName] += localCount;
    }  
  }
  
  //map: processes count in detectors
  
  std::map<G4String,G4int>::const_iterator itp2;
  for ( itp2 = localRun->fProcCounter2.begin();
        itp2 != localRun->fProcCounter2.end(); ++itp2 ) {

    G4String procName = itp2->first;
    G4int localCount = itp2->second;
    if ( fProcCounter2.find(procName) == fProcCounter2.end()) {
      fProcCounter2[procName] = localCount;
    }
    else {
      fProcCounter2[procName] += localCount;
    }  
  }

  std::map<G4String, G4int>::const_iterator itp3;
  for (itp3 = localRun->fProcCounter3.begin();
      itp3 != localRun->fProcCounter3.end(); ++itp3) {

      G4String procName = itp3->first;
      G4int localCount = itp3->second;
      if (fProcCounter3.find(procName) == fProcCounter3.end()) {
          fProcCounter3[procName] = localCount;
      }
      else {
          fProcCounter3[procName] += localCount;
      }
  }
    
  //map: created particles in target   
  std::map<G4String,ParticleData>::const_iterator itc;
  for (itc = localRun->fParticleDataMap1.begin(); 
       itc != localRun->fParticleDataMap1.end(); ++itc) {
    
    G4String name = itc->first;
    const ParticleData& localData = itc->second;   
    if ( fParticleDataMap1.find(name) == fParticleDataMap1.end()) {
      fParticleDataMap1[name]
       = ParticleData(localData.fCount, 
                      localData.fEmean, 
                      localData.fEmin, 
                      localData.fEmax);
    }
    else {
      ParticleData& data = fParticleDataMap1[name];   
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin) data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax) data.fEmax = emax; 
    }   
  }
  
  //map: created particle in detectors       
  std::map<G4String,ParticleData>::const_iterator itn1;
  for (itn1 = localRun->fParticleDataMap2.begin(); 
       itn1 != localRun->fParticleDataMap2.end(); ++itn1) {
    
    G4String name = itn1->first;
    const ParticleData& localData = itn1->second;   
    if ( fParticleDataMap2.find(name) == fParticleDataMap2.end()) {
      fParticleDataMap2[name]
       = ParticleData(localData.fCount, 
                      localData.fEmean, 
                      localData.fEmin, 
                      localData.fEmax);
    }
    else {
      ParticleData& data = fParticleDataMap2[name];   
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin) data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax) data.fEmax = emax; 
    }   
  }

  std::map<G4String, ParticleData>::const_iterator itn2;
  for (itn2 = localRun->fParticleDataMap3.begin();
      itn2 != localRun->fParticleDataMap3.end(); ++itn2) {

      G4String name = itn2->first;
      const ParticleData& localData = itn2->second;
      if (fParticleDataMap3.find(name) == fParticleDataMap3.end()) {
          fParticleDataMap3[name]
              = ParticleData(localData.fCount,
                  localData.fEmean,
                  localData.fEmin,
                  localData.fEmax);
      }
      else {
          ParticleData& data = fParticleDataMap3[name];
          data.fCount += localData.fCount;
          data.fEmean += localData.fEmean;
          G4double emin = localData.fEmin;
          if (emin < data.fEmin) data.fEmin = emin;
          G4double emax = localData.fEmax;
          if (emax > data.fEmax) data.fEmax = emax;
      }

      // Merge fProcCounterByVolume 
      for (const auto& volPair : localRun->fProcCounterByVolume) {
          const G4String& volumeName = volPair.first;
          const std::map<G4String, G4int>& localProcMap = volPair.second;
          for (const auto& procPair : localProcMap) {
              const G4String& procName = procPair.first;
              G4int localCount = procPair.second;
              // Add the count to the corresponding entry in our master map.
              fProcCounterByVolume[volumeName][procName] += localCount;
          }
      }
  }

  // Merge annihilation counts from each thread
  for (const auto& kv : localRun->fAnnihilationByVolumeFromSource) {
      fAnnihilationByVolumeFromSource[kv.first] += kv.second;
  }
  // ——— Merge our per-detector 1274 keV photon counters ———
  fEnterDisk1 += localRun->fEnterDisk1;
  fReachDet1 += localRun->fReachDet1;
  fEnterDisk2 += localRun->fEnterDisk2;
  fReachDet2 += localRun->fReachDet2;
  fEnterKapton1 += localRun->fEnterKapton1;
  fEnterKapton2 += localRun->fEnterKapton2;

  G4Run::Merge(run); 
} 

void Run::PrintProcessFrequency()
{
    G4cout << "\n Process calls frequency by volume:" << G4endl;
    for (const auto& volPair : fProcCounterByVolume) {
        G4String volName = volPair.first;
        G4cout << "\n Volume: " << volName << G4endl;
        for (const auto& procPair : volPair.second) {
            G4String procName = procPair.first;
            G4int count = procPair.second;
            G4cout << "   " << procName << " = " << count << G4endl;
        }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun() 
{
  G4int prec = 5, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);
  
  // run condition
  //   
  G4String Particle = fParticle->GetParticleName();    
  G4cout << "\n The run is " << numberOfEvent << " "<< Particle << " of "
         << G4BestUnit(fEkin,"Energy") << " through : ";
          
  G4cout << "\n Target   : Length = " 
         << G4BestUnit(fDetector->GetTargetLength(),"Length")
         << " Radius    = " 
         << G4BestUnit(fDetector->GetTargetRadius(),"Length")  
         << " Material = " 
         << fDetector->GetTargetMaterial()->GetName();
  G4cout << "\n Detector : Length = " 
         << G4BestUnit(fDetector->GetDetectorLength(),"Length")
         << " Thickness = " 
         << G4BestUnit(fDetector->GetDetectorThickness(),"Length")  
         << " Material = " 
         << fDetector->GetDetectorMaterial()->GetName() << G4endl;

  if (numberOfEvent == 0) { G4cout.precision(dfprec);   return;}
  
  // compute mean Energy deposited and rms in target
  //
  G4int TotNbofEvents = numberOfEvent;
  fEdepTarget /= TotNbofEvents; fEdepTarget2 /= TotNbofEvents;
  G4double rmsEdep = fEdepTarget2 - fEdepTarget*fEdepTarget;
  if (rmsEdep>0.) rmsEdep = std::sqrt(rmsEdep);
  else            rmsEdep = 0.;
  
  G4cout << "\n Mean energy deposit in target,   in time window = "
         << G4BestUnit(fEdepTarget,"Energy") << ";  rms = "
         << G4BestUnit(rmsEdep,    "Energy") 
         << G4endl;

  // compute mean Energy deposited and rms in detectors
  //
  fEdepDetect11 /= TotNbofEvents; fEdepDetect12 /= TotNbofEvents;
  rmsEdep = fEdepDetect12 - fEdepDetect11*fEdepDetect11;
  if (rmsEdep>0.) rmsEdep = std::sqrt(rmsEdep);
  else            rmsEdep = 0.;
  
  G4cout << " Mean energy deposit in detector 1, in time window = "
         << G4BestUnit(fEdepDetect11,"Energy") << ";  rms = "
         << G4BestUnit(rmsEdep,    "Energy") 
         << G4endl;

  fEdepDetect21 /= TotNbofEvents; fEdepDetect22 /= TotNbofEvents;
  rmsEdep = fEdepDetect22 - fEdepDetect21 * fEdepDetect21;
  if (rmsEdep > 0.) rmsEdep = std::sqrt(rmsEdep);
  else            rmsEdep = 0.;

  G4cout << " Mean energy deposit in detector 2, in time window = "
      << G4BestUnit(fEdepDetect21, "Energy") << ";  rms = "
      << G4BestUnit(rmsEdep, "Energy")
      << G4endl;

  // frequency of processes in target
  //
  G4cout << "\n Process calls frequency in target :" << G4endl;
  G4int index = 0;
  std::map<G4String,G4int>::iterator it1;    
  for (it1 = fProcCounter1.begin(); it1 != fProcCounter1.end(); it1++) {
     G4String procName = it1->first;
     G4int    count    = it1->second;
     G4String space = " "; if (++index%3 == 0) space = "\n";
     G4cout << " " << std::setw(20) << procName << "="<< std::setw(7) << count
            << space;
  }
  G4cout << G4endl;
  
  // frequency of processes in detectors
  // 1
  G4cout << "\n Process calls frequency in detector 1:" << G4endl;
  index = 0;
  std::map<G4String,G4int>::iterator it2;    
  for (it2 = fProcCounter2.begin(); it2 != fProcCounter2.end(); it2++) {
     G4String procName = it2->first;
     G4int    count    = it2->second;
     G4String space = " "; if (++index%3 == 0) space = "\n";
     G4cout << " " << std::setw(20) << procName << "="<< std::setw(7) << count
            << space;
  }
  G4cout << G4endl;
  // 2
  G4cout << "\n Process calls frequency in detector 2:" << G4endl;
  index = 0;
  std::map<G4String, G4int>::iterator it3;
  for (it3 = fProcCounter3.begin(); it3 != fProcCounter3.end(); it3++) {
      G4String procName = it3->first;
      G4int    count = it3->second;
      G4String space = " "; if (++index % 3 == 0) space = "\n";
      G4cout << " " << std::setw(20) << procName << "=" << std::setw(7) << count
          << space;
  }
  G4cout << G4endl;

  // TODO - NEW PRINTING OF FREQUENCY MIGHT NEED TO DELETE THE ONES BEFORE
  PrintProcessFrequency();
    
  // particles count in target
  //
  G4cout << "\n List of generated particles in target:" << G4endl;
     
  std::map<G4String,ParticleData>::iterator itn1;               
  for (itn1 = fParticleDataMap1.begin(); itn1 != fParticleDataMap1.end(); itn1++) {
     G4String name = itn1->first;
     ParticleData data = itn1->second;
     G4int count = data.fCount;
     G4double eMean = data.fEmean/count;
     G4double eMin = data.fEmin;
     G4double eMax = data.fEmax;    
         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ")" << G4endl;           
 }

 // particles count in detectors
 // 1
 G4cout << "\n List of generated particles in detector 1:" << G4endl;
     
 std::map<G4String,ParticleData>::iterator itn2;               
 for (itn2 = fParticleDataMap2.begin(); itn2 != fParticleDataMap2.end(); itn2++) {
    G4String name = itn2->first;
    ParticleData data = itn2->second;
    G4int count = data.fCount;
    G4double eMean = data.fEmean/count;
    G4double eMin = data.fEmin;
    G4double eMax = data.fEmax;
         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") << ")" << G4endl; 
  }
  G4cout << G4endl;

  //2 
  G4cout << "\n List of generated particles in detector 2:" << G4endl;

  std::map<G4String, ParticleData>::iterator itn3;
  for (itn3 = fParticleDataMap3.begin(); itn3 != fParticleDataMap3.end(); itn3++) {
      G4String name = itn3->first;
      ParticleData data = itn3->second;
      G4int count = data.fCount;
      G4double eMean = data.fEmean / count;
      G4double eMin = data.fEmin;
      G4double eMax = data.fEmax;

      G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
          << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
          << "\t( " << G4BestUnit(eMin, "Energy")
          << " --> " << G4BestUnit(eMax, "Energy") << ")" << G4endl;
  }
  G4cout << G4endl;

  G4cout << "\n=========================================" << G4endl;
  G4cout << " Positron Statistics: " << G4endl;
  G4cout << "-----------------------------------------" << G4endl;
  G4cout << " Positrons Reaching Detector 1 from target: " << fPositronReachingDetector1 << G4endl;
  G4cout << " Positrons Reaching Detector 2 from target: " << fPositronReachingDetector2 << G4endl;
  G4cout << "=========================================" << G4endl;


  // ——— 1274 keV photon blocking summary ———
  G4cout << "\n--- 1274 keV gamma blocking ---\n";

  // Detector 1 (forward)
  G4int e1 = GetEnterDisk1(), r1 = GetReachDet1();
  G4cout
      << " Det1: entered disk = " << e1
      << ", reached det = " << r1;
  if (e1 > 0) {
      G4double f1 = 1.0 - double(r1) / double(e1);
      G4cout << ", blocked = " << f1 * 100 << " %";
  }
  G4cout << G4endl;

  // Detector 2 (backward)
  G4int e2 = GetEnterDisk2(), r2 = GetReachDet2();
  G4cout
      << " Det2: entered disk = " << e2
      << ", reached det = " << r2;
  if (e2 > 0) {
      G4double f2 = 1.0 - double(r2) / double(e2);
      G4cout << ", blocked = " << f2 * 100 << " %";
  }
  G4cout << "\n-----------------------------\n" << G4endl;

  // Kapton summary - 1274 keV photons crossing
  {
      G4int k1 = GetEnterKapton1();
      G4int k2 = GetEnterKapton2();
      G4cout << "\n--- Kapton disk crossings ---\n"
          << " Kapton1 (forward): " << k1 << G4endl
          << " Kapton2 (backward): " << k2 << "\n"
          << "------------------------------\n" << G4endl;
  }

  // --- Print annihilations by volume ---
  G4cout << "\n Annihilations of target-born e+ by volume:" << G4endl;
  for (const auto& kv : fAnnihilationByVolumeFromSource) {
      G4cout << "  Volume: "
          << std::setw(12) << kv.first
          << " -> " << kv.second << " annihilations"
          << G4endl;
  }

  G4cout << "\n=== Photon flow per volume (+-10 keV windows) ===\n";
  G4cout << std::left << std::setw(22) << "Volume"
      << std::right << std::setw(12) << "in_511"
      << std::setw(12) << "out_511"
      << std::setw(12) << "abs_511"
      << std::setw(12) << "in_1274"
      << std::setw(12) << "out_1274"
      << std::setw(12) << "abs_1274" << G4endl;

  for (const auto& kv : fSimpleFlow) {
      const auto& name = kv.first;
      const auto& c = kv.second;
      G4cout << std::left << std::setw(22) << name
          << std::right << std::setw(12) << c.in511
          << std::setw(12) << c.out511
          << std::setw(12) << (c.in511 - c.out511)
          << std::setw(12) << c.in1274
          << std::setw(12) << c.out1274
          << std::setw(12) << (c.in1274 - c.out1274) << G4endl;
  }
  G4cout << "===============================================\n";

 
  // activities in VR mode
  //
  WriteActivity(numberOfEvent);
 
  //remove all contents in fProcCounter, fCount 
  fProcCounter1.clear();
  fProcCounter2.clear();
  fProcCounter3.clear();
  fParticleDataMap1.clear();    
  fParticleDataMap2.clear();
  fParticleDataMap3.clear();
                          
  //restore default format         
  G4cout.precision(dfprec);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEnterGamma(const G4String& volName, G4bool is511, G4bool is1274) {
    auto& c = fSimpleFlow[volName];
    if (is511)  ++c.in511;
    if (is1274) ++c.in1274;
}

void Run::AddExitGamma(const G4String& volName, G4bool is511, G4bool is1274) {
    auto& c = fSimpleFlow[volName];
    if (is511)  ++c.out511;
    if (is1274) ++c.out1274;
}

void Run::WriteActivity(G4int nevent)
{
 G4ProcessTable *pTable = G4ProcessTable::GetProcessTable();
 G4Radioactivation* rDecay = (G4Radioactivation *)
         pTable->FindProcess("Radioactivation", "GenericIon");
   
 // output the induced radioactivities (in VR mode only)
 //
 if ((rDecay == 0) || (rDecay->IsAnalogueMonteCarlo())) return;
 
 G4String fileName = G4AnalysisManager::Instance()->GetFileName() + ".activity";
 std::ofstream outfile (fileName, std::ios::out );
 
 std::vector<G4RadioactivityTable*> theTables =
                              rDecay->GetTheRadioactivityTables();

 for (size_t i = 0 ; i < theTables.size(); i++) {
    G4double rate, error;
    outfile << "Radioactivities in decay window no. " << i << G4endl;
    outfile << "Z \tA \tE \tActivity (decays/window) \tError (decays/window) "
            << G4endl;

    map<G4ThreeVector,G4TwoVector> *aMap = theTables[i]->GetTheMap();
    map<G4ThreeVector,G4TwoVector>::iterator iter;
    for (iter=aMap->begin(); iter != aMap->end(); iter++) {
       rate = iter->second.x()/nevent;
       error = std::sqrt(iter->second.y())/nevent;
       if (rate < 0.) rate = 0.;                // statically it can be < 0.
       outfile << iter->first.x() <<"\t"<< iter->first.y() <<"\t"
               << iter->first.z() << "\t" << rate <<"\t" << error << G4endl;
    }
    outfile << G4endl;
 }
 outfile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
