#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <iostream>

void plot_gamma_histogram(const std::string &filename) {
    // Open the ROOT file
    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Get the histogram from the file
    TH1F *hist = dynamic_cast<TH1F*>(file->Get("gamma_energy_histogram"));
    if (!hist) {
        std::cerr << "Error: Histogram 'gamma_energy_histogram' not found in " << filename << std::endl;
        file->Close();
        return;
    }

    // Create a canvas to draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Gamma Energy Histogram", 800, 600);
    hist->Draw();

    // Save the canvas as a PNG file
    canvas->SaveAs("gamma_energy_histogram.png");

    // Clean up
    delete canvas;
    file->Close();
}

int main() {
    // Your code for plotting histograms
    return 0;
}
