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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  
  fParticleGun->SetParticleEnergy(0*eV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,0.));


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// TODO - change to decay or calibration
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // Part for Calibration of 511 and 1274 keV photons
    auto gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(gamma);
	fParticleGun->SetParticleEnergy(1275 * keV); // TODO change to 511 or 1274
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
    G4double cosTheta = 2.0 * G4UniformRand() - 1.0; 
    G4double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    G4double phi = 2.0 * CLHEP::pi * G4UniformRand();
    G4ThreeVector direction(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    //simulation of one e+

    //fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("e+"));
    //fParticleGun->SetParticleEnergy(1. * MeV);  // Slow positron
    //fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));  // Optional
    //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));  // Optional
    //fParticleGun->GeneratePrimaryVertex(anEvent);


    // Part for Na-22 decay

    //auto ion = G4IonTable::GetIonTable()->GetIon(11, 22, 0. * keV);
    //fParticleGun->SetParticleDefinition(ion);
    //fParticleGun->SetParticleCharge(0. * eplus);
    //fParticleGun->SetParticleEnergy(0. * eV);                           // at rest
    //fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, 0));
    //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1)); // dir irrelevant at 0 eV
    //fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

