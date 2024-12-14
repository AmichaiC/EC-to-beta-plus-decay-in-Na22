// File: plotHisto.C
// Purpose: Plot histograms and ntuples generated in Geant4 simulation

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"

void plotHisto() {
    // Open the ROOT file containing histograms
    TFile* file = TFile::Open("C:\\Users\\amich\\Study\\rdecay02\\build\\Release\\rdecay02.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open the ROOT file!" << std::endl;
        return;
    }

    // Retrieve histograms from the file
    TH1D* h11 = (TH1D*)file->Get("H11");
    TH1D* h17 = (TH1D*)file->Get("H17");
    TH1D* h18 = (TH1D*)file->Get("H18");

    if (!h11 || !h17 || !h18) {
        std::cerr << "Error: One or more histograms could not be found!" << std::endl;
        file->Close();
        return;
    }

    // Create a canvas to draw histograms
    TCanvas* c1 = new TCanvas("c1", "Histograms", 800, 600);
    c1->Divide(2, 2); // Divide canvas into sub-pads

    // Style settings
    gStyle->SetOptStat(0); // Disable the statistics box

    // Draw H11 histogram
    c1->cd(1);
    h11->SetLineColor(kBlue);
    h11->SetLineWidth(2);
    h11->SetTitle("H11: Energy Deposit");
    h11->GetXaxis()->SetTitle("Energy (MeV)");
    h11->GetYaxis()->SetTitle("Counts");
    h11->Draw();

    // Draw H17 histogram
    c1->cd(2);
    h17->SetLineColor(kRed);
    h17->SetLineWidth(2);
    h17->SetTitle("H17: Energy Deposit");
    h17->GetXaxis()->SetTitle("Energy (MeV)");
    h17->GetYaxis()->SetTitle("Counts");
    h17->Draw();

    // Draw H18 histogram
    c1->cd(3);
    h18->SetLineColor(kGreen);
    h18->SetLineWidth(2);
    h18->SetTitle("H18: Energy Deposit");
    h18->GetXaxis()->SetTitle("Energy (MeV)");
    h18->GetYaxis()->SetTitle("Counts");
    h18->Draw();

    // Save the canvas to a file
    c1->SaveAs("histograms.png");

    // Clean up
    file->Close();
}
