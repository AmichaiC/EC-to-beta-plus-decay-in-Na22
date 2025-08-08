#ifndef MYTRACKINFO_HH
#define MYTRACKINFO_HH

#include "G4VUserTrackInformation.hh"
#include "G4LogicalVolume.hh"

class MyTrackInfo : public G4VUserTrackInformation {
public:
    // Initialize flags and birth volume pointer
    MyTrackInfo()
        : fCounted(false)
        , fAnnihilationCounted(false)
        , fBirthVolume(nullptr) {}
    virtual ~MyTrackInfo() {}

    // Detector crossing count
    void SetCounted(G4bool flag) { fCounted = flag; }
    G4bool IsCounted() const { return fCounted; }

    // Unique annihilation counting
    void SetAnnihilationCounted(G4bool flag) { fAnnihilationCounted = flag; }
    G4bool IsAnnihilationCounted() const { return fAnnihilationCounted; }

    // Record and retrieve the volume where the track was born
    void SetBirthVolume(G4LogicalVolume* vol) { fBirthVolume = vol; }
    G4LogicalVolume* GetBirthVolume() const { return fBirthVolume; }

    G4bool IsCountedKapton1() const { return fCountedKapton1; }
    void SetCountedKapton1(G4bool val = true) { fCountedKapton1 = val; }

    G4bool IsCountedDisk1() const { return fCountedDisk1; }
    void SetCountedDisk1(G4bool val = true) { fCountedDisk1 = val; }

    G4bool IsCountedDet1() const { return fCountedDet1; }
    void SetCountedDet1(G4bool val = true) { fCountedDet1 = val; }

    // --- Negative side ---
    G4bool IsCountedKapton2() const { return fCountedKapton2; }
    void SetCountedKapton2(G4bool val = true) { fCountedKapton2 = val; }

    G4bool IsCountedDisk2() const { return fCountedDisk2; }
    void SetCountedDisk2(G4bool val = true) { fCountedDisk2 = val; }

    G4bool IsCountedDet2() const { return fCountedDet2; }
    void SetCountedDet2(G4bool val = true) { fCountedDet2 = val; }

private:
    G4bool              fCounted;               // flag for detector crossing
    G4bool              fAnnihilationCounted;   // flag for counting annihilation once
    G4LogicalVolume* fBirthVolume;           // birth volume of this track
    G4bool fCountedKapton1 = false;
    G4bool fCountedDisk1 = false;
    G4bool fCountedDet1 = false;
    G4bool fCountedKapton2 = false;
    G4bool fCountedDisk2 = false;
    G4bool fCountedDet2 = false;
};

#endif // MYTRACKINFO_HH
