# QualityControl
Quality control for TQS data

This code performs a quality-control check on output from [Skyline] ( https://skyline.gs.washington.edu/labkey/project/home/software/Skyline/begin.view) .  You will also need "Master list" that notes which transition should be used for quantification.  The output is a .csv of sample ID, compound ID, any flags and vetted area.

Please cite the following paper when using this code:



# Software
This code is written in R .

It requires the following packages available on CRAN: 
* plyr
* reshape2

# Example to test
Input files for testing code are: 
*  HILIC_MasterList_Example.csv
* ExampleSkylineOutput.csv

Output should match the file
 *  QC_outputExampleSkylineOutput.csv

#  Acknowledgements
This code was vetted by Katherine Heal.



