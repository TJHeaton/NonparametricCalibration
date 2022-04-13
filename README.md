# NonparametricCalibration
Code for "Non-parametric calibration of multiple related radiocarbon determinations and their calendar age summarisation"


13th April 2022

This folder contain code to perform non-parametric calibration and subsequent summarisation of radiocarbon dates

The code here provides code as used for JRSS C submission. All files need to be placed in a directory as organised here. 

The essential functions for users who wish to calibrate their own datasets are:

1) NPCalibrationWalker.R
2) NPCalibrationNeal.R

In these files, one should be able to calibrate any set of samples, and also estimate the joint calendar age density, by changing x and xsig at the start of each file. I would suggest trying the Walker version first as it is slightly faster if you have a huge number of samples to jointly calibrate. The current files calibrate the data from Kerr et al. 

I have added a small progress bar to the code to let a user know it is running (and how far along it is) and have also tried to create some nice plotting functions which show the summarised joint density, and also the individual posterior calendar age for any individual determination. You can access the other variables fairly easily too - hopefully they are sensibly named so they are self-explanatory.

Any feedback on the code is welcome. I intend to try and integrate it with OxCal (https://c14.arch.ox.ac.uk/oxcal.html) so that it is widely available to the radiocarbon community.


---- Simulation Study - IF YOU WISH TO RUN THEM THEY ARE SLOW

The simulation study code can be found in the three files:

1) SimulationStudySingleCluster_Revisedv1.R
2) SimulationStudyThreeCluster_Revisedv1.R
3) SimulationStudyUniform_Revisedv1.R

These generate output in a subdirectory called SimOutput/ so you will need to create this first. The code will print to screen the values in Table 1 but you can then also create the boxplots using the code in:

BoxplotImprovements.R

---- Paper Examples (Known Truth)

These can be recreated by the files:

1) TestDensityEstimateNormal_Revisedv1.R
2) TestDensityEstimateUniform_Revisedv1.R
3) TestDensityMexicoPopulation.R
4) TestDensityMexicoPopulationShifted.R

---- Paper Examples (Real-Life)

Use the code in:

1) RealDatasetExamples.R
2) RealExamplesWalkerPostProcessingFinal.R









