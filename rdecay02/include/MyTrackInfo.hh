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

private:
    G4bool              fCounted;               // flag for detector crossing
    G4bool              fAnnihilationCounted;   // flag for counting annihilation once
    G4LogicalVolume* fBirthVolume;           // birth volume of this track
};

#endif // MYTRACKINFO_HH
