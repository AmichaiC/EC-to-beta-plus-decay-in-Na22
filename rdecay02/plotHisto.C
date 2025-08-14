// plotHisto.C  — minimal ROOT macro, no arguments needed
#include <iostream>
#include <fstream>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"

// ---- SETTINGS (edit these two lines if needed) ----
static const char* kRootPath =
"C:/Users/amich/Study/rdecay02/build/Release/rdecay02.root";
static const bool  kWriteCSV = true;   // set false if you don't want CSVs
// ---------------------------------------------------

static void writeCSV(TH1* h, const std::string& path) {
    if (!h) return;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) { std::cout << "CSV: cannot open " << path << "\n"; return; }
    out << "BinCenter,Count\n";
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        out << h->GetBinCenter(i) << "," << h->GetBinContent(i) << "\n";
}

void plotHisto() {
    std::cout << "Opening: " << kRootPath << "\n";
    TFile* file = TFile::Open(kRootPath, "READ");
    if (!file || file->IsZombie()) { std::cerr << "Error: cannot open file\n"; return; }

    // where to save outputs
    std::string rootDir = kRootPath;
    for (auto& ch : rootDir) if (ch == '\\') ch = '/';
    auto pos = rootDir.find_last_of('/');
    std::string outDir = (pos == std::string::npos) ? "." : rootDir.substr(0, pos);
    gSystem->mkdir(outDir.c_str(), kTRUE);

    // nice style
    gStyle->SetOptStat(0);

    // list histograms we might have (it’s OK if some are missing)
    const char* names[] = { "H11","H16","H17","H18","H19" };
    const char* titles[] = { "H11","H16","H17","H18","H19" };

    for (int i = 0; i < 5; ++i) {
        TH1* h = nullptr;
        file->GetObject(names[i], h);
        if (!h) { std::cout << "Note: histogram \"" << names[i] << "\" not found (skipping)\n"; continue; }

        // draw to its own canvas
        std::string cname = std::string("c_") + names[i];
        TCanvas* c = new TCanvas(cname.c_str(), titles[i], 900, 600);
        h->SetLineWidth(2);
        h->Draw("HIST");
        c->Update();

        // save PNG
        std::string png = outDir + "/" + std::string(names[i]) + ".png";
        c->SaveAs(png.c_str());
        std::cout << "Saved " << png << "\n";

        // optional CSV
        if (kWriteCSV) {
            std::string csv = outDir + "/" + std::string(names[i]) + "_data.csv";
            writeCSV(h, csv);
            std::cout << "Wrote " << csv << "\n";
        }
    }

    file->Close();
    std::cout << "Done.\n";
}
