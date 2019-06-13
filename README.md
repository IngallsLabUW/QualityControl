# QualityControl

This code performs a quality-control check on output from [Skyline](https://skyline.gs.washington.edu/labkey/project/home/software/Skyline/begin.view). 

The repository is split into two sections: quality control for a Thermo Q Exactive HF (Orbitrap) and a Waters Xevo TQ-S (triple quadrupole) mass spectrometer.

**Quality control for QE data**
This code requires the use ofthe command line. An example of how to use it is below:
'''R
$ Rscript QE_Targeted_QC.R <filepath of Skyline output file> <filepath of blank matcher file>
'''

You will then be guided through the necessary steps to set parameters for the quality control.

**Quality control for TQS data**

You will also need a "Master list" that notes which transition should be used for quantification.  

Please cite the following paper when using this code:

# Setup
Run the samples should be run in the following manner for the quality control and [B-MIS normalization] (https://github.com/IngallsLabUW/B-MIS-normalization):

* Standards run (all mixed) at least once at the beginning and end of the run, labeled:  Date_Std_AdditionalID (e.g. 161018_Std_FirstStandardinH20)
* Standards run (in representative matrix, all mixed) at least once the beginning and end of the run, labeled Date_Std_AdditionalID (e.g. 161018_Std_FirstStandardinMatrix)
*  Blanks run (preferably method/filter blanks) at least once labeled: Date_Blk_AdditionalID (e.g. 161018_Blk_FirstBlank)
* A pooled sample run at least three times throughout the run, labeled: Date_Poo_AdditionalID_Rep (e.g. 161018_Poo_PooledSample_1)
* Samples, labeled: Date_Smp_AdditionalID_Rep (e.g. 161018_Std_FirstSample_BioRep1)

# Software
Both the QE and the TQS code is written in R. The QE quality control code requires the use of the command line. 

It requires the following packages available on CRAN: 
* plyr
* reshape2

## Example to test
Input files for testing code are: 
*  HILIC_MasterList_Example.csv
* ExampleSkylineOutput.csv

Output should match the file
 *  QC_outputExampleSkylineOutput.csv



#  Acknowledgements
This code was vetted by Katherine Heal.



