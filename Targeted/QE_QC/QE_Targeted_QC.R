library(plyr)
library(rlist)
library(tidyverse)

# This script is intended to be used for .csv files output from Skyline (https://skyline.ms/project/home/begin.view?), and an accompanying
# csv that matches the samples of your run to their appropriate blanks. See included file in Github repository for reference.

# Run this script using the command line. When in the proper directory, enter Rscript <this script> <skyline output csv> <blank sample matching csv>
# The script will then guide you through the steps to obtaining your results, filtered for quality control!

PromptToContinue <- function() {
  # Asks the user if they are ready to continue to the next step. 
  repeat {
    cat("Ready to move on? (Y/N): " )
    answer <- tolower(readLines("stdin", n = 1))
    if (answer == "y")  {
      break
    } else if (answer == "n") {
      stop("Not ready to move on.")
    }
  }
}
PromptForStdTags <- function() {
  # Asks the user to set tag matches by entering sample replicate names. 
  # These tags are used to find standards for the creation of reliable retention times.
  std.tags <- list()
  repeat {
    cat("Sample tag for reliable retention time (leave blank if finished): ")
    tag <- readLines("stdin", n = 1);
    if (tag == "") {
      return(std.tags)
    }
    std.tags <- list.append(std.tags, tag)
  }
}
CheckBlankMatcher <- function(blank.matcher) {
  # Takes a blank matcher file and separates any multi-value variable
  # columns into their own row.
  #
  # Args:
  #   blank.matcher: CSV entered by user to match samples with 
  #   appropriate blanks.
  #
  # Returns:
  #   blank.matcher: new CSV with any duplicate values separated
  #   into their own rows.
  #
  blank.matcher <- do.call("rbind", Map(data.frame, 
    Blank.Name=strsplit(as.character(blank.matcher$Blank.Name), ","), 
    Replicate.Name=(blank.matcher$Replicate.Name))
  )
  blank.matcher <- blank.matcher[c(2, 1)]
  
  return(blank.matcher)
}
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
IdentifyRunTypes <- function(skyline.output) {
  # Identify run types and return each unique value present in the Skyline output.
  #
  # Args
  #   skyline.output: Raw output file from Skyline.
  #
  # Returns
  #   run.type: list of labels identifying the run types, isolated from Replicate.Name.
  #   Options conssist of samples (smp), pooled (poo), standards (std), and blanks (blk).
  #
  run.type <- tolower(str_extract(skyline.output$Replicate.Name, "(?<=_)[^_]+(?=_)"))
  cat("The replicate types in this run are:", "\n")
  print(toString(unique(run.type)))
} 

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
    filter(Replicate.Name %in% checked.blanks$Replicate.Name) %>%
    mutate(SN.Flag       = ifelse(((Area / Background) < SN.min), "SN.Flag", NA)) %>%
    mutate(ppm.Flag      = ifelse(abs(Mass.Error.PPM) > ppm.flex, "ppm.Flag", NA)) %>%
    mutate(area.min.Flag = ifelse((Area < area.min), "area.min.Flag", NA))
  
  return(first.flags)
}
CreateRTFlags <- function(skyline.output, std.tags) {
  # Creates a dataset with Retention Time (RT) flags.
  #
  # Args
  #   skyline.output: Output file from machine to be modified.
  #   std.tags: User-entered tag matches from replicate names.
  #
  # Returns
  #   Dataset containing columns of flags corresponding to values
  #   falling outside of user-entered parameters when compared 
  #   to a standard or template.
  #
  retention.time.flags <- skyline.output %>%
    filter(Replicate.Name %in% std.tags) %>%
    select(Mass.Feature, Retention.Time) %>%
    group_by(Mass.Feature) %>%
    summarise(RT.Reference = mean((Retention.Time), na.rm = TRUE))
  
  return(retention.time.flags)
}
CreateBlankFlags <- function(skyline.output, checked.blanks) {
  # Creates a dataset with blank flags.
  #
  # Args
  #   skyline.output: Output file from machine to be modified.
  #   checked.blanks: File with matches which samples go with which blanks.
  #
  # Returns
  #   Dataset containing columns of flags marking which replicates 
  #   are blanks.
  #
  blank.flags <- skyline.output %>%
    filter(Replicate.Name %in% checked.blanks$Blank.Name) %>%
    rename(Blank.Name = Replicate.Name,
           Blank.Area = Area) %>%
    select(Blank.Name, Mass.Feature, Blank.Area) %>%
    left_join(checked.blanks, by = "Blank.Name") %>% 
    select(-Blank.Name) %>%
    arrange(desc(Blank.Area)) %>%
    group_by(Mass.Feature, Replicate.Name) %>% 
    filter(row_number() == 1)
  
  return(blank.flags)
}


## Processing----------------------------------------------

# Main function for command line input.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Please enter a csv of Skyline output and sample blanks to compare.")
}

# Assign files to variables within script.
skyline.output <- read.csv(file = args[1], header = TRUE)
blank.matcher  <- read.csv(file = args[2], header = TRUE)
input.file     <- basename(args[1])

# Set parameters for quality control processing.
cat("Pick the minimum height to be counted as a 'real' peak (QE suggestion: HILIC - 1000, Cyano - 5000): " )
area.min        <- as.double(readLines("stdin", n = 1))
cat("Pick retention time (RT) flexibility (QE suggestion: +/- 0.4 min for HILIC, +/- 0.2 min for Cyano): ")
RT.flex         <- as.double(readLines("stdin", n = 1))
cat("Pick signal size comparison between sample and blank to merit inclusion (QE suggestion: +/- 0.2): ")
blank.ratio.max <- as.double(readLines("stdin", n = 1))
cat("Pick acceptable signal to noise ratio value. Note: broader peaks create more background noise(QE suggestion: 5 for Cyano, 4 for HILIC.): ")
SN.min          <- as.double(readLines("stdin", n = 1))
cat("Pick an absolute value for a cutoff for parts per million (ppm) (QE suggestion: 7): ")
ppm.flex        <- as.double(readLines("stdin", n = 1))

cat("Your minimum height value is:", area.min, "\n")
cat("Your retention time flexibility value is:", RT.flex, "\n")
cat("Your blank threshold value is:", blank.ratio.max, "\n")
cat("Your signal to noise ratio value is:", SN.min, "\n")
cat("Your parts per million value is:", ppm.flex, "\n")
PromptToContinue()

# Create standard tags for matching which samples go with which blanks. 
# Double check with user that the tags are correct.
std.tags <- PromptForStdTags()
printed.tags <- toString(std.tags)
cat("Your tags for sample matching are:", printed.tags, "\n")
PromptToContinue()

# Reorder blank matcher csv to ensure there is one value per variable.
#testing <- read.csv("./Targeted/QE_QC/datafiles/Samps_with_Blanks2.csv", header = TRUE)
checked.blanks <- CheckBlankMatcher(blank.matcher)

# Drop unnecessary colums, transform and rename variables as needed.
skyline.columns.dropped <- skyline.output %>%
  select(-Protein.Name, -Protein) 
skyline.classes.transformed <- TransformVariables(skyline.columns.dropped)
skyline.runtypes.identified <- IdentifyRunTypes(skyline.output)

# Create datasets for different flag types.
SNPPMAM.flags <- CreateFirstFlags(skyline.classes.transformed, area.min, SN.min, ppm.flex)
RT.flags      <- CreateRTFlags(skyline.classes.transformed, std.tags)
blank.flags   <- CreateBlankFlags(skyline.classes.transformed, checked.blanks)

# Joining datasets---------------------------------------

# Combine signal to noise, ppm, area min, and RT flags into one dataset.
first.join <- SNPPMAM.flags %>%
  left_join(RT.flags) %>%
  mutate(Retention.Time = as.numeric(as.character(Retention.Time))) %>%
  mutate(RT.Flag        = ifelse((abs(Retention.Time - RT.Reference) > RT.flex), "RT.Flag", NA))

# Join blank flags to the dataset.
second.join <- first.join %>%
  left_join(blank.flags) %>%
  mutate(blank.Flag = ifelse((Area / Blank.Area) < blank.ratio.max, "blank.Flag", NA))

# Finally, combine all the flags and throw out any peak with a flag.
last.join <- second.join %>%
  mutate(all.Flags      = paste(SN.Flag, ppm.Flag, area.min.Flag, RT.Flag, blank.Flag, sep = ", ")) %>%
  mutate(all.Flags      = as.character(all.Flags %>% str_remove_all("NA, ") %>%  str_remove_all("NA"))) %>%
  mutate(Area.with.QC   = ifelse(str_detect(all.Flags, "Flag"), NA, Area)) 

# If there are any standards, add those to the bottom of the dataset.
Stds.test <- grepl("_Std_", skyline.output$Replicate.Name)

if (any(Stds.test == TRUE)) {
  print("There are standards in this run. Joining standard samples to the bottom of the dataset!", quote = FALSE)
  standards <- skyline.classes.transformed[grep("Std", skyline.classes.transformed$Replicate.Name), ]
  last.join <- rbind.fill(last.join, standards)
  } else {
  print("No standards exist in this set.")
}

# Print to file with comments and new name!
con <- file(paste("QEQC_", input.file, sep = ""), open = "wt")
writeLines(paste("Hello! Welcome to the world of QE Quality Control! ",
                 "Minimum area for a real peak: ", area.min, ". ",
                 "RT flexibility: ", RT.flex, ". ",
                 "Blank can be this fraction of a sample: ", blank.ratio.max, ". ",
                 "S/N ratio: " , SN.min, ". ",
                 "Parts per million flexibility: ", ppm.flex, ". ",
                 "Sample tag matches:", printed.tags
                 , ". ",
                 "Processed on: ", Sys.time(), ". ",
                 sep = ""), con)
write.csv(last.join, con, row.names = FALSE)
close(con)
