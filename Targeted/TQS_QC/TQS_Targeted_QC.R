## First draft of a new TQS Targeted QC
# June 24th, 2019
# Regina Lionheart

library(plyr)
library(reshape2)
library(tidyverse)
options(scipen=999)


# TODO (rlionheart): Add style guide and move functions to top.
# TODO (kheal, aboysen): Would it be nice to have all these tables accessible for download (ie, the RT.Range.Table, Area.Table, whatever)? Helpful to see an overview of the run?
# Alternately, I could include a section that will print a custom summary of the run with whatever is interesting.
# TODO (rlionheart, kheal, aboysen): Why is the IR table giving different data than the original skyline TQS QC?

# Import files
input_file <- "./Targeted/TQS_QC/datafiles/ExampleSkylineOutput_TQS.csv"
areas.raw <- read.csv("./Targeted/TQS_QC/datafiles/ExampleSkylineOutput_TQS.csv", row.names = NULL, header = TRUE) %>%
        select(-X)
master <- read.csv("./Targeted/TQS_QC/datafiles/HILIC_MasterList_Example.csv") %>%
        rename(Second.Trace = X2nd.trace)

## Set the parameters for the QC -----------
max.height <- 1.0e8
min.height <- 1000
area.min <- 1000
RT.flex <- 0.4
IR.flex <- 0.3
blk.thresh <- 0.3
SN.min <- 4

## ID run types
run.type <- tolower(str_extract(areas.raw$Replicate.Name, "(?<=_)[^_]+(?=_)"))

## Change variables classes as needed, add Sample.Type column
# TODO (rlionheart): add documentation
TransformVariables <- function(skyline.output) {
     before <- lapply(skyline.output, class)
     cat("Original class variables ", "\n")
     print(paste(colnames(skyline.output), ":", before))
     
     skyline.output <- skyline.output %>%
          mutate(Replicate.Name     = suppressWarnings(as.character(Replicate.Name))) %>%
          mutate(Precursor.Ion.Name = suppressWarnings(as.factor(Precursor.Ion.Name))) %>%
          mutate(Area               = suppressWarnings(as.numeric(as.character(Area)))) %>%
          mutate(Retention.Time     = suppressWarnings(as.numeric(as.character(Retention.Time)))) %>%
          mutate(Background         = suppressWarnings(as.numeric(as.character(Background)))) %>%
          mutate(Height             = suppressWarnings(as.numeric(as.character(Height)))) %>%
          mutate(Sample.Type        = suppressWarnings(as.factor(run.type)))
     
     after <- lapply(skyline.output, class)
     cat("New class variables ", "\n")
     print(paste(colnames(skyline.output), ":", after))
     
     return(skyline.output)
}
areas.transformed <- TransformVariables(areas.raw)


# Ion Ratio Ranges Table  ----------------------------------------------------
# TODO (rlionheart): add documentation to CheckFragments function. 
# TODO (rlionheart): add explanation of what is happening to find the IR range table.
# TODO (rlionheart): Document and clean up function

CheckFragments <- function(areas.transformed) { 

        fragment.check <- areas.transformed %>%
        filter(Sample.Type == "std") %>%
        select(Replicate.Name, Precursor.Ion.Name, Precursor.Mz, Product.Mz, Sample.Type)

        fragment.unique <- fragment.unique <-unique(fragment.check %>% select(Precursor.Ion.Name, Precursor.Mz, Product.Mz))

        fragment.multi.unique <- fragment.unique %>%
                count(Precursor.Ion.Name) %>%
                mutate(Two.Fragments = ifelse((n==1), FALSE, TRUE)) %>%
                select(-n)

        fragments.checked <- fragment.check %>%
                left_join(fragment.multi.unique, by = "Precursor.Ion.Name" ) %>%
                merge(y = master,
                        by.x = c("Precursor.Ion.Name", "Product.Mz"),
                        by.y = c("Compound.Name", "Daughter"),
                        all.x = TRUE) %>%
                select(Replicate.Name, Precursor.Ion.Name, Precursor.Mz, Product.Mz, Two.Fragments, Quan.Trace, Second.Trace) %>%
                mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
                mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
                mutate(QT.Five.Percent = ifelse((Two.Fragments == TRUE & Quan.Trace == TRUE), 0.05 * Product.Mz, NA)) %>%
                mutate(Significant.Size = QT.Five.Percent < Product.Mz) %>%
                mutate(Std.Type = str_match(Replicate.Name, "\\d+_\\w+_(\\w+)_\\d+")[,2]) %>%
                group_by(Std.Type, Precursor.Ion.Name, Product.Mz) %>%
                summarise_all(first) %>%
                arrange(Precursor.Ion.Name)
        
        return(fragments.checked)
}

fragments.checked <- CheckFragments(areas.transformed)
IR.Table <- fragments.checked %>%
        group_by(Precursor.Ion.Name, Std.Type) %>%
        mutate(IR.Ratio = ifelse(TRUE %in% Significant.Size, (Product.Mz[Quan.Trace == TRUE] / Product.Mz[Second.Trace == TRUE]), NA)) %>%
        select(Precursor.Ion.Name, Std.Type, IR.Ratio)  %>%
        unique()

        

## Retention Time  ----------------------------------------------------
RT.Range.Table <- areas.transformed %>%
        select(Precursor.Ion.Name, Retention.Time, Sample.Type) %>%
        filter(Sample.Type == "std") %>%
        group_by(Precursor.Ion.Name) %>%
        mutate(RT.min = min(Retention.Time, na.rm = TRUE)) %>%
        mutate(RT.max = max(Retention.Time, na.rm = TRUE)) %>%
        mutate(RT.Reference = mean(Retention.Time, na.rm = TRUE)) %>%
        select(-Retention.Time) %>%
        unique()


## Blanks ---------------------------------------
## if the area of the blank is more than blk.thresh% of the
## area in the sample, make a flag
Blank.Table <- areas.transformed %>%
        select(Precursor.Ion.Name, Area, Sample.Type) %>%
        filter(Sample.Type == "blk") %>%
        group_by(Precursor.Ion.Name) %>%
        mutate(Blank.max = max(Area, na.rm = TRUE)) %>%
        select(-Area) %>%
        unique()


## Height  ---------------------------------------
Height.Table <- areas.transformed %>%
        select(Replicate.Name, Sample.Type, Precursor.Ion.Name, Height) %>%
        filter(Sample.Type == "smp")

## Area  ---------------------------------------
Area.Table <- areas.transformed %>%
        select(Replicate.Name, Sample.Type, Precursor.Ion.Name, Area) %>%
        filter(Sample.Type == "smp") %>%
        mutate(area.min.Flag = ifelse((Area < area.min), "area.min.Flag", NA))


## Signal to Noise  ---------------------------------------
# TODO (rlionheart): Ditto question for S.N?
SN.Table <- areas.transformed %>%
        select(Replicate.Name, Sample.Type, Precursor.Ion.Name, Area, Background) %>%
        filter(Sample.Type == "smp") %>%
        mutate(Signal.to.Noise = (Area / Background)) %>%
        mutate(SN.Flag = ifelse((Signal.to.Noise < SN.min), "SN.Flag", NA))


## Join datasets  ---------------------------------------

## Construct final output
# TODO (rlionheart): Add IR info & flags
# IR.ok[1,] <- IR.range[1,]*(1-IR.flex)
# IR.ok[2,] <- IR.range[2,]*(1+IR.flex)
# 
# IR.flags.added <- areas.transformed %>%
#         filter(Sample.Type %in% c("smp", "poo")) %>%
#         left_join(IR.Table, by = "Precursor.Ion.Name") 
####
        

RT.flags.added <- areas.transformed %>%
        filter(Sample.Type %in% c("smp", "poo")) %>%
        left_join(RT.Range.Table, by = "Precursor.Ion.Name") %>%
        mutate(RT.Flag = ifelse((abs(Retention.Time - RT.Reference) > RT.flex), "RT.Flag", NA)) %>%
        select(Replicate.Name, Precursor.Ion.Name:Background, Height, RT.Flag) 
        
Blank.flags.added <- RT.flags.added %>%
        left_join(Blank.Table, by = "Precursor.Ion.Name") %>%
        mutate(Blank.Reference = mean(Area, na.rm = TRUE)) %>%
        mutate(blank.Flag = ifelse((Area / Blank.Reference) < blk.thresh, "blank.Flag", NA)) %>%
        select(Replicate.Name:RT.Flag, blank.Flag)

Height.flags.added <- Blank.flags.added %>%
        mutate(height.min.Flag = ifelse((Height < min.height), "height.min.Flag", NA)) %>%
        mutate(overloaded.Flag = ifelse((Height > max.height), "overloaded.Flag", NA)) 

Area.flags.added <- Height.flags.added %>%
        mutate(area.min.Flag = ifelse((Area < area.min), "area.min.Flag", NA))

SN.flags.added <- Area.flags.added %>%
        mutate(Signal.to.Noise = (Area / Background)) %>%
        mutate(SN.Flag = ifelse((Signal.to.Noise < SN.min), "SN.Flag", NA)) %>%
        select(Replicate.Name:area.min.Flag, SN.Flag)
        
final.table <- SN.flags.added %>%
        mutate(all.Flags      = paste(RT.Flag, blank.Flag, height.min.Flag, overloaded.Flag, area.min.Flag, SN.Flag, sep = ", ")) %>%
        mutate(all.Flags      = as.character(all.Flags %>% str_remove_all("NA, ") %>%  str_remove_all("NA")))

        
# If there are any standards, add those to the bottom of the dataset.
Stds.test <- grepl("_Std_", areas.raw$Replicate.Name)

if (any(Stds.test == TRUE)) {
        print("There are standards in this run. Joining standard samples to the bottom of the dataset!", quote = FALSE)
        standards <- areas.transformed[grep("Std", areas.transformed$Replicate.Name), ]
        final.table <- rbind.fill(final.table, standards)
} else {
        print("No standards exist in this set.")
}


# Print to file with comments and new name!`    `
con <- file(paste("TQSQC_", basename(input_file), sep = ""), open = "wt")
writeLines(paste("Hello! Welcome to the world of TQS Quality Control! ",
                 "Maximum height for peaks: ", max.height, ". ",
                 "Minimum height for a real peak: ", min.height, ". ",
                 "Minimum area for a real peak: ", area.min, ". ",
                 "RT flexibility: ", RT.flex, ". ",
                 "Ion ratio (IR) flexibility: ", IR.flex, ". ",
                 "Blank can be this fraction of a sample: ", blk.thresh, ". ",
                 "S/N ratio: " , SN.min, ". ",
                 "Processed on: ", Sys.time(), ". ",
                 sep = ""), con)
write.csv(final.table, con, row.names = FALSE)
close(con)





