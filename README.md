System: Na22 source in the middle, from both sides Ge detectors with tungsten cones connecting the source to the detectors, and in between from both sides kapton and tungsten disks to annihilate the positrons that are emitted.
Changing the source - 3 options: Na22, Photon (511 or 1274 keV) or positrons in PrimaryGeneratorAction.cc in the method GeneratePrimaries, uncomment the needed part and comment out the rest.

Changing the distance from the source to the detectors and the thickness of the Kapton and Tungsten disks in DetectorConstruction.cc all marked in TODO comment.
Distance - fdistanceFromTarToDet in line 75.
Tungsten - diskThickness in line 390
Kapton - kThickness in line 432 
It is also possible to change the disks to annulus and to add magnetic fields - both options not used in the simulation but with changing Boolean values or inner radius it can be changed.

The code is compatible with running on different threads, although in my computer it didn't show significant improvement which probably indicates a bug which I didn't find.

In order to visualize the code use vis.mac.
In order to run the code use run.mac.

The code creates root files which can be used with the file plotHisto.C to show the distribution of energy of the photons detected in graphs and in excel files of both detectors. H11 is for detector 1, and the second is the same in H19.

.exe file is being created with those excel files in build\Release.

results from the simulation can be seen under 5Kapton with the the thickness of the Kapton and tungsten including excel files and text files which has also the results of the simulation - Particles detected in different parts of the system, events and more...
