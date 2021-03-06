# Quality Control

This code performs a user-defined quality-control check on output from [Skyline](https://skyline.ms/project/home/software/Skyline/begin.view).

The repository is split into two sections: targeted and untargeted metabolomic analysis. Within those sections, choose between code for Thermo Q Exactive HF (Orbitrap) and a Waters Xevo TQ-S (triple quadrupole) mass spectrometers.

#### Software

Both the QE and the TQS code is written in R. The QE quality control code requires the use of the command line.

#### LCMS Setup

Samples should be run in the following manner for the quality control and [B-MIS normalization](https://github.com/IngallsLabUW/B-MIS-normalization):

* Please label all samples in the following manner: Date_RunType_AdditionalID (e.g. 161018_Std_FirstStandardinH20). RunType refers to whether the sample is a standard (Std), sample (Smp), pooled (poo), or blank (blk).
* Standards run (all mixed) at least once at the beginning and end of the run.
* Standards run (in representative matrix, all mixed) at least once the beginning and end of the run. Example label: 161019_Std_FirstStandardinMatrix
* Blanks run (preferably method/filter blanks) at least once. Example label: 161018_Blk_FirstBlank
* A pooled sample run at least three times throughout the run. Example label:161018_Poo_PooledSample_1
* Samples. Example label: Date_Smp_AdditionalID_Rep

*This repository is still a work in progress! Thanks for your patience.*

## Targeted Data Quality Control

#### Required Input Files

QE | TQS
------------ | -------------
CSV of skyline output | CSV of skyline output
CSV pairing blank replicates with their appropriate samples | CSV master list of compounds for matching traces.

For an example of the blank/samples csv file, please see Targeted/QE_QC/datafiles/Samps_With_Blanks.csv.
For an example of a master list csv file, please see Targeted/TQS_QC/datafiles/HILIC_MasterList_Example.csv.

Please note that for the TQS, the master list should have the following columns: *Group, Compound Name, Parent, Daughter, Retention Time, Parent Charge, Daughter Charge, Notes, Function Number, Cone_voltage, Collision_energy, Quan Trace, 2nd trace*


***

### Q Exactive HF

QE quality control uses the command line. An example of how to use the targeted quality control for the Qexactive is below: 
```shell
$ Rscript QE_Targeted_QC.R <filename of Skyline output> <filename of matching blanks-to-samples csv>
```
After hitting enter, you will then be guided through the necessary steps to set parameters for the quality control.

**Output**

The output will be saved as "QEQC_ + original skyline output filename" to your working directory. For example, HILICpos.csv would be saved as QEQC_HILICpos.csv.

**Required packages available on CRAN**

```R
library(plyr)
library(rlist)
library(tidyverse)
```

***

### TQ-S Quality Control

TQS quality control is written as a normal R program, requiring manual input of file paths for the Skyline output and the appropriate master list. Please note that this manual input will be changed in a future release. 

**Output**

Similar to the QE Quality Control, this output will be saved as "TQSQC_ + original skyline output filename" to your working directory. HILICpos.csv would be saved as TQSQC_HILICpos.csv.

Remember to set your numerical output options to undo the default scientific notation setting. This will make your values much easier to read.

**Required packages available on CRAN**

```R
library(plyr)
library(reshape2)
library(tidyverse)
options(scipen=999)
```
*This section is still in progress. We appreciate your patience!*

***


## Untargeted Data Quality Control 

*This section is still in progress. We are working on completing the code for untargeted analysis. We appreciate your patience!*


### Acknowledgements
This code was vetted by Katherine Heal.





