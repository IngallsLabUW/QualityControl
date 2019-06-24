## First draft of a new TQS Targeted QC
# June 24th, 2019
# Regina Lionheart

library(reshape2)
library(tidyverse)

# Import files
# TODO (rlionheart): Decide on best import practice. Upload csv via shiny app (https://shiny.rstudio.com/gallery/file-upload.html)?
areas.raw <- read.csv("./Targeted/TQS_QC/datafiles/ExampleSkylineOutput_TQS.csv", header = TRUE)
master <- read.csv("./Targeted/TQS_QC/datafiles/HILIC_MasterList_Example.csv")

## Set the parameters for the QC -----------
max.height <- 1.0e8
min.height <- 1000
RT.flex <- 0.4
IR.flex <- 0.3
blk.thresh <- 0.3
SN.thresh <- 4

## ID run types
run.type <- tolower(str_extract(areas.raw$Replicate.Name, "(?<=_)[^_]+(?=_)"))

## Change variables classes as needed.
areas.raw$Replicate.Name <- as.character(areas.raw$Replicate.Name)
areas.raw$Precursor.Ion.Name <- as.factor(areas.raw$Precursor.Ion.Name)
areas.raw$Area <- as.numeric(areas.raw$Area)
areas.raw$Retention.Time <- as.numeric(areas.raw$Retention.Time)
areas.raw$Background <- as.numeric(areas.raw$Background)
areas.raw$Height <- as.numeric(areas.raw$Height)

## Add sample type column 
areas.raw$Sample.Type <- as.factor(run.type)

# Retention Time  ----------------------------------------------------

# For each precursor ion name, get the retention time range.

RT.range_test <- areas.raw %>%
     select(Precursor.Ion.Name, Retention.Time, Sample.Type)  %>%
     group_by(Sample.Type) %>%
     filter(Sample.Type == "std") %>%
     group_by(Precursor.Ion.Name) %>%
     summarize(
          minRT = min(Retention.Time),
          maxRT = max(Retention.Time)
     )

###########################################
# TODO (rlionheart): This is leftover from OG tqs QC
stdRows <- areas.raw$Sample.Type == "std"
blkRows <- areas.raw$Sample.Type == "blk"

areas.split<-split(areas.raw, areas.raw$Sample.Type)

#type.options <- names(areas.split)

cmpd.blk.list<- split(areas.split[["blk"]],
                      areas.split[["blk"]]$Precursor.Ion.Name)


##########################################


     





