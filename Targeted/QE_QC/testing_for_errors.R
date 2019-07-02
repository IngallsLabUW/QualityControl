library(rlist)
library(plyr)
library(stringr)
library(tidyverse)



original.output <- read.csv("./Targeted/QE_QC/datafiles/HILICPos_DiatomBatch1_SkylineExport.csv", header = TRUE) %>%
     select(-Protein.Name, -Protein) 
blank.matcher <- read.csv("./Targeted/QE_QC/datafiles/NSamps_with_Blanks.csv", header = TRUE)

area.min <- 1000
RT.flex <- 0.4
blank.ratio.max <- 0.2
SN.min <- 4
ppm.flex <- 7


TransformVariables <- function(skyline.output) {
     # Transforms variable classes and renames troublesome columns. 
     # Prints variables classes before and after to ensure success.
     #
     # Args 
     #  skyline.output: Raw output file from Skyline. 
     #
     # Returns
     #   skyline.output: File of the same name, with chosen columns
     #   changed to numeric classes, and Precursor.Ion.Name renamed
     #   to the clearer Mass.Feature.
     #
     before <- lapply(skyline.output, class)
     cat("Original class variables ", "\n")
     print(paste(colnames(skyline.output), ":", before))
     
     skyline.output <- skyline.output %>%
          mutate(Replicate.Name = suppressWarnings(as.character(Replicate.Name))) %>%
          mutate(Retention.Time = suppressWarnings(as.numeric(as.character(Retention.Time)))) %>%
          mutate(Area           = suppressWarnings(as.numeric(as.character(Area)))) %>%
          mutate(Background     = suppressWarnings(as.numeric(as.character(Background)))) %>%
          mutate(Mass.Error.PPM = suppressWarnings(as.numeric(as.character(Mass.Error.PPM)))) %>%
          rename(Mass.Feature   = Precursor.Ion.Name)
     
     after <- lapply(skyline.output, class)
     cat("New class variables ", "\n")
     print(paste(colnames(skyline.output), ":", after))
     
     return(skyline.output)
}
skyline.variables.transformed <- TransformVariables(original.output)

CreateFirstFlags <- function(skyline.output, area.min, SN.min, ppm.flex) {
     # Creates a dataset containing signal to noise (SN), Parts per million (ppm), 
     # and Area minimum (area.min) flags.
     # 
     # Args
     #   skyline.output: Output file from machine to be modified.
     #   area.min: User-entered value for minimum area to be accepted.
     #   SN.min: User-entered value for acceptable signal to noise ratio.
     #   ppm.flex: User-entered value for parts per million cutoff.
     #
     # Returns
     #   Dataset containing columns of flags corresponding to values
     #   that fall outside of the user-entered parameters for signal to noise,
     #   parts per million, and minimum area.
     #
     first.flags <- skyline.output %>%
          mutate(SN.Flag       = ifelse(((Area / Background) < SN.min), "SN.Flag", NA)) %>%
          mutate(ppm.Flag      = ifelse((Mass.Error.PPM > ppm.flex), "ppm.Flag", NA)) %>%
          mutate(area.min.Flag = ifelse((Area < area.min), "area.min.Flag", NA))
     
     return(first.flags)
} # Area, SN, ppm flags
firstflagsmade <- CreateFirstFlags(skyline.variables.transformed, area.min, SN.min, ppm.flex)

CreateRTFlags <- function(skyline.output) {
        ############################################################
        # CHECK OPTIONAL STD.TAGS FUNCTIONALITY#
        ############################################################
        # Creates a dataset with Retention Time (RT) reference flags.
        #
        # Args
        #   skyline.output: Output file from machine to be modified.
        #
        # Returns
        #   Dataset containing columns of flags corresponding to values
        #   falling outside of user-entered parameters when compared 
        #   to a standard or template.
        #
        retention.time.flags <- skyline.output %>%
                select(Mass.Feature, Retention.Time) %>%
                group_by(Mass.Feature) %>%
                summarise(RT.Reference = mean((Retention.Time), na.rm = TRUE))
        
        return(retention.time.flags)
}
rtflagsmade <- CreateRTFlags(skyline.variables.transformed)

CreateBlankFlags <- function(skyline.output, blank.matcher) {
        # Creates a dataset with blank flags.
        #
        # Args
        #   skyline.output: Output file from machine to be modified.
        #   blank.matcher: File with matches which samples go with which blanks.
        #
        # Returns
        #   Dataset containing columns of flags marking which replicates
        #   are blanks.
        #
        blank.flags <- skyline.output %>%
                filter(Replicate.Name %in% blank.matcher$Blank.Name) %>%
                rename(Blank.Name = Replicate.Name,
                       Blank.Area = Area) %>%
                select(Blank.Name, Mass.Feature, Blank.Area) %>%
                  left_join(blank.matcher) %>% select(-Blank.Name)

        return(blank.flags)
}
blkflagsmade <- CreateBlankFlags(skyline.variables.transformed, blank.matcher)



# Combine signal to noise, ppm, area min, and RT flags into one dataset.
first.join <- firstflagsmade %>%
        left_join(rtflagsmade) %>%
        mutate(Retention.Time = as.numeric(as.character(Retention.Time))) %>%
        mutate(RT.Flag        = ifelse((abs(Retention.Time - RT.Reference) > RT.flex), "RT.Flag", NA))

# Join blank flags to the dataset.
second.join <- first.join %>%
        left_join(blkflagsmade) %>%
        mutate(blank.Flag = ifelse((Area / Blank.Area) < blank.ratio.max, "blank.Flag", NA))

# Finally, combine all the flags.
last.join <- second.join %>%
        mutate(all.Flags      = paste(SN.Flag, ppm.Flag, area.min.Flag, RT.Flag, blank.Flag, sep = ", ")) %>%
        mutate(all.Flags      = as.character(all.Flags %>% str_remove_all("NA, ") %>%  str_remove_all("NA"))) %>%
        mutate(Area.with.QC   = ifelse(str_detect(all.Flags, "Flag"), NA, Area))
        

# If there are any standards, add those to the bottom of the dataset.
Stds.test <- grepl("_Std_", original.output$Replicate.Name)

if (any(Stds.test = TRUE)) {
        print("There are standards in this run. Joining standard samples to the bottom of the dataset!", quote = FALSE)
        standards <- skyline.variables.transformed[grep("Std", skyline.variables.transformed$Replicate.Name), ]
        last.join <- rbind.fill(last.join, standards)
} else {
        print("No standards exist in this set.")
}

toporiginal <- head(original.output, 20L)
toplast <- head(last.join, 20L)
testbind <- rbind.fill(toporiginal, toplast)


# AM I CRAZY --------------------------------------------------------------

Has.Standards <- read.csv("./Targeted/QE_QC/datafiles/Natalie_ODZ_DepthProfile_pos_Report.csv", header = TRUE) %>%
        select(-Protein.Name, -Protein) 
No.Standards <- read.csv("./Targeted/QE_QC/datafiles/HILICPos_DiatomBatch1_SkylineExport.csv", header = TRUE) %>%
        select(-Protein.Name, -Protein) 

HS.test <- grepl("_Std_", Has.Standards$Replicate.Name)
NS.test <- grepl("_Std_", No.Standards$Replicate.Name)

which(HS.test, arr.ind = TRUE)
which(NS.test, arr.ind = TRUE)

any(HS.test == TRUE)
any(NS.test == TRUE)

