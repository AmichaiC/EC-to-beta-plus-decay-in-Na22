#ifndef MyTrackInfo_h
#define MyTrackInfo_h 1

#include "G4VUserTrackInformation.hh"

class MyTrackInfo : public G4VUserTrackInformation {
public:
    MyTrackInfo() : fCounted(false) {}
    virtual ~MyTrackInfo() {}

    void SetCounted(G4bool flag) { fCounted = flag; }
    G4bool IsCounted() const { return fCounted; }

private:
    G4bool fCounted; // flag indicating whether the track has been counted as reaching a detector
};

#endif
