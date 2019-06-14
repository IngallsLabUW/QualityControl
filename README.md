# Quality Control

This code performs a user-defined quality-control check on output from [Skyline](https://skyline.ms/project/home/software/Skyline/begin.view).

The repository is split into two sections: quality control for targeted and untargeted metabolomic analysis. Within those sections, choose between code for Thermo Q Exactive HF (Orbitrap) and a Waters Xevo TQ-S (triple quadrupole) mass spectrometers.

#### Software

Both the QE and the TQS code is written in R. The QE quality control code requires the use of the command line.

#### Setup

Run the samples should be run in the following manner for the quality control and [B-MIS normalization](https://github.com/IngallsLabUW/B-MIS-normalization):

* Standards run (all mixed) at least once at the beginning and end of the run, labeled: Date_Std_AdditionalID (e.g. 161018_Std_FirstStandardinH20)
* Standards run (in representative matrix, all mixed) at least once the beginning and end of the run, labeled Date_Std_AdditionalID (e.g. 161018_Std_FirstStandardinMatrix)
* Blanks run (preferably method/filter blanks) at least once labeled: Date_Blk_AdditionalID (e.g. 161018_Blk_FirstBlank)
* A pooled sample run at least three times throughout the run, labeled: Date_Poo_AdditionalID_Rep (e.g. 161018_Poo_PooledSample_1)
* Samples, labeled: Date_Smp_AdditionalID_Rep (e.g. 161018_Std_FirstSample_BioRep1)

#### Required Input Files

QE | TQS
------------ | -------------
CSV of skyline output | *in progress*
CSV for identifying blank replicates| *in progress*

*This repository is still a work in progress! We are working on completing the code for untargeted analysis. Thanks for your patience.*

***

## Targeted Data Quality Control

### Q Exactive HF

QE quality control uses the command line. An example of how to use the targeted quality control for the Qexactive is below: 
```shell
$ Rscript QE_Targeted_QC.R <filename of Skyline output> <filename of blank matching csv>
```
After hitting enter, you will then be guided through the necessary steps to set parameters for the quality control.

**Output**

The output will be saved as "QEQC_ + original skyline output filename" to your working directory. For example, HILICpos.csv would be saved as QEQC_HILICpos.csv.

**Required packages available on CRAN**

```R
library(rlist)
library(tidyverse)
```

***

### TQ-S Quality Control

*This section is still in progress. We appreciate your patience!*

***


## Untargeted Data Quality Control 

*This section is still in progress. We appreciate your patience!*


### Acknowledgements
This code was vetted by Katherine Heal.
Please cite the following paper when using this code:




