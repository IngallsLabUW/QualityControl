library(plyr)
library(reshape2)
library(stringr)
library(tidyverse)
options(scipen=999)

TransformVariables <- function(skyline.output) {
  # Transforms variable classes and renames troublesome columns. 
  # Prints variables classes before and after to ensure success.
  #
  # Args 
  #  skyline.output: Raw output file from Skyline. 
  #
  # Returns
  #   skyline.output: File of the same name, with chosen columns
  #   changed to numeric classes
  #
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
CheckStdFragments <- function(areas.transformed) { 
  # Modifies transformed skyline output to prepare for standard ion ratio detection by running several tests.
  # 1. Isolate standards.
  # 2. Isolate unique Product.Mz fragments per compound.
  # 3. Check if each compound has two unique fragments.
  # 4. If it does, identify which fragment is quantitative and which is secondary by comparing them to the master compound list.
  # 5. Find 5% of the quantitative fragment.
  # 6. Determine if the secondary trace is > the 5% value from step 5.
  #
  # Args:
  #   areas.transformed: Output from skyline that has had its variables modified to numeric values.
  #
  # Returns:
  #   fragments.checked: Modified data frame with added columns reflecting the above tests.
  #
  fragment.check <- areas.transformed %>%
    filter(Sample.Type == "std") %>%
    select(Replicate.Name, Precursor.Ion.Name, Area, Precursor.Mz, Product.Mz, Sample.Type)
  
  fragment.unique <-unique(fragment.check %>% select(Precursor.Ion.Name, Precursor.Mz, Product.Mz))
  
  fragment.multi.unique <- fragment.unique %>%
    count(Precursor.Ion.Name) %>%
    mutate(Two.Fragments = ifelse((n==1), FALSE, TRUE)) %>%
    select(-n)
  
  fragments.checked <- fragment.check %>%
    left_join(fragment.multi.unique, by = "Precursor.Ion.Name" ) %>%
    left_join(standard.types, by = "Replicate.Name") %>%
    merge(y = master,
          by.x = c("Precursor.Ion.Name", "Product.Mz"),
          by.y = c("Compound.Name", "Daughter"),
          all.x = TRUE) %>%
    select(Replicate.Name, Std.Type, Precursor.Ion.Name, Area, Precursor.Mz, Product.Mz, Two.Fragments, Quan.Trace, Second.Trace) %>%
    mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
    mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
    mutate(QT.Five.Percent = ifelse((Two.Fragments == TRUE & Quan.Trace == TRUE), 0.05 * Product.Mz, NA)) %>%
    mutate(Significant.Size = QT.Five.Percent < Product.Mz) %>%
    group_by(Replicate.Name, Precursor.Ion.Name, Product.Mz) %>% 
    summarise_all(first) %>% 
    arrange(Precursor.Ion.Name) %>%
    select(Replicate.Name, Std.Type, Precursor.Ion.Name, Precursor.Mz, Product.Mz, Area, Two.Fragments:Significant.Size)
  
  return(fragments.checked)
}
CheckSmpFragments <- function(areas.transformed) {
  # Modifies transformed skyline output to prepare for standard ion ratio detection by running several tests.
  # 1. Isolate samples and pooled runs.
  # 2. Isolate unique Product.Mz fragments per compound.
  # 3. Check if each compound has two unique fragments.
  # 4. If it does, identify which fragment is quantitative and which is secondary by comparing them to the master compound list.
  # 5. Find 5% of the quantitative fragment.
  # 6. Determine if the secondary trace is > the 5% value from step 5.
  # 7. Find ion ratio by dividing the area of the quantitative trace by the area of the secondary trace.
  #
  # Args:
  #   areas.transformed: Output from skyline that has had its variables modified to numeric values.
  #
  # Returns:
  #   all.samples.IR: Modified data frame with added columns reflecting the above tests.
  unique.smp.frags <- unique(all.samples %>% select(Precursor.Ion.Name, Precursor.Mz, Product.Mz))
  
  unique.smp.frags2 <- unique.smp.frags %>%
    count(Precursor.Ion.Name) %>%
    mutate(Two.Fragments = ifelse((n==1), FALSE, TRUE)) %>%
    select(-n)
  
  all.samples.IR <- all.samples %>%
    left_join(unique.smp.frags2, by = "Precursor.Ion.Name" ) %>%
    merge(y = master,
          by.x = c("Precursor.Ion.Name", "Product.Mz"),
          by.y = c("Compound.Name", "Daughter"),
          all.x = TRUE) %>%
    select(Replicate.Name, Precursor.Ion.Name, Protein.Name, Precursor.Mz, Product.Mz, Area, Two.Fragments, Quan.Trace, Second.Trace) %>%
    mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
    mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
    mutate(QT.Five.Percent = ifelse((Two.Fragments == TRUE & Quan.Trace == TRUE), 0.05 * Product.Mz, NA)) %>%
    mutate(Significant.Size = QT.Five.Percent < Product.Mz) %>%
    group_by(Precursor.Ion.Name) %>%
    mutate(IR.Ratio = ifelse(TRUE %in% Significant.Size, (Area[Quan.Trace == TRUE] / Area[Second.Trace == TRUE]), NA))
  return(all.samples.IR)
}

# Import files - this import format will be changed when integrated with Shiny
input.file <- "./Targeted/TQS_QC/datafiles/ReRuns_TQS_Transition_Results.csv"
areas.raw  <- read.csv("./Targeted/TQS_QC/datafiles/ReRuns_TQS_Transition_Results.csv", row.names = NULL, header = TRUE) #%>% select(-X)
master     <- read.csv("./Targeted/TQS_QC/datafiles/HILIC_TQS_GBT_MasterList.csv") %>% rename(Second.Trace = X2nd.trace)

standard.types <- read.csv("./Targeted/TQS_QC/datafiles/GBT_filter.csv") %>%
  mutate(Replicate.Name = suppressWarnings(as.character(Replicate.Name))) %>%
  mutate(Std.Type       = suppressWarnings(as.character(Std.Type)))


# Set the parameters for the QC
max.height <- 1.0e8
min.height <- 1000
area.min   <- 1000
RT.flex    <- 0.4
IR.flex    <- 0.3
blk.thresh <- 0.3
SN.min     <- 4

std.tags <- c()
standard.pattern <- "Glycine Betaine" # or custom pattern without the standard types csv.

# ID run types
run.type <- tolower(str_extract(areas.raw$Replicate.Name, "(?<=_)[^_]+(?=_)"))

# Change variables classes as needed, add Sample.Type column
areas.transformed <- TransformVariables(areas.raw)

# Ion Ratio Ranges Table  ----------------------------------------------------
# Find Ion Ratio by dividing the area of the quantitative trace by the area of the secondary trace. 
# Find the minimum and maximum IR to create reference table of IR ranges.
IR.Table <- CheckStdFragments(areas.transformed) %>%
  # filter(Replicate.Name %in% std.tags) %>%
  # mutate(Std.Type2 = sapply(strsplit(Replicate.Name, "_"), "[[", 3)) %>% # An optional way to isolate the mix type, if naming conventions are consistent.
  group_by(Precursor.Ion.Name, Replicate.Name) %>%
  mutate(Std.Ion.Ratio = ifelse(Quan.Trace == TRUE, (Area[Quan.Trace == TRUE]) / (Area[Second.Trace == TRUE]), NA)) %>%
  mutate(Diff.Mix = ifelse((standard.pattern %in% str_extract(Precursor.Ion.Name, standard.pattern)), TRUE, FALSE)) %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(IR.min = ifelse(Diff.Mix == TRUE, suppressWarnings(min(Std.Ion.Ratio[Std.Type == "GBTMix"], na.rm = TRUE)), suppressWarnings(min(Std.Ion.Ratio[Std.Type == "Mix1"], na.rm = TRUE)))) %>%
  mutate(IR.max = ifelse(Diff.Mix == TRUE, suppressWarnings(max(Std.Ion.Ratio[Std.Type == "GBTMix"], na.rm = TRUE)), suppressWarnings(max(Std.Ion.Ratio[Std.Type == "Mix1"], na.rm = TRUE)))) %>%
  #mutate(IR.min = min(Std.Ion.Ratio, na.rm = TRUE)) %>%
  #mutate(IR.max = max(Std.Ion.Ratio, na.rm = TRUE)) %>%
  select(Precursor.Ion.Name, IR.min, IR.max) %>%
  unique()


# Retention Time Table ----------------------------------------------------
# Find the minimum and maximum Retention Times and take the average.
# Use this as a reference table for acceptable Retention Times.
RT.Range.Table <- areas.transformed %>%
  select(Precursor.Ion.Name, Retention.Time, Sample.Type) %>%
  filter(Sample.Type == "std") %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(RT.min = min(Retention.Time, na.rm = TRUE)) %>%
  mutate(RT.max = max(Retention.Time, na.rm = TRUE)) %>%
  mutate(RT.Reference = mean(Retention.Time, na.rm = TRUE)) %>%
  select(-Retention.Time) %>%
  unique()

# Blanks ---------------------------------------
# Isolate the blanks in the sample and add a column 
# with maximum blank for each Precursor ion name.
Blank.Table <- areas.transformed %>%
  merge(y = master,
        by.x = c("Precursor.Ion.Name", "Product.Mz"),
        by.y = c("Compound.Name", "Daughter"),
        all.x = TRUE) %>%
  mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
  mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
  filter(Sample.Type == "blk") %>%
  filter(Quan.Trace == TRUE) %>%
  select(Precursor.Ion.Name, Replicate.Name, Area, Sample.Type) %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(Blank.max = max(Area, na.rm = TRUE)) %>%
  select(-Area) %>%
  unique()

# Height  ---------------------------------------
# Isolate all pooled and sample Heights
Height.Table <- areas.transformed %>%
  select(Replicate.Name, Sample.Type, Precursor.Ion.Name, Precursor.Mz, Product.Mz, Height) %>%
  filter(Sample.Type %in% c("smp", "poo")) 

# Area  ---------------------------------------
# Isolate all pooled and sample Areas.
Area.Table <- areas.transformed %>%
  select(Replicate.Name, Sample.Type, Precursor.Ion.Name, Area) %>%
  filter(Sample.Type %in% c("smp", "poo")) 


# Signal to Noise  ---------------------------------------
# Isolate all pooled and sample runs. Find the Signal to Noise
# by dividing the Background of each run by its Area.
SN.Table <- areas.transformed %>%
  ##
  merge(y = master,
        by.x = c("Precursor.Ion.Name", "Product.Mz"),
        by.y = c("Compound.Name", "Daughter"),
        all.x = TRUE) %>%
  mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
  mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
  filter(Sample.Type %in% c("smp", "poo")) %>%
  filter(Quan.Trace == TRUE) %>%
  select(Replicate.Name, Sample.Type, Precursor.Ion.Name, Area, Background) %>%
  mutate(Signal.to.Noise = (Area / Background)) 

# Construct final output of sample and pooled runs.
all.samples <- areas.transformed %>%
  filter(Sample.Type %in% c("smp", "poo"))

# Ion Ratio Flags  ---------------------------------------
# If the Ion Ratio falls outside of the IR.Table range +/- the
# IR.flex value, add a flag.
IR.flags.added <- CheckSmpFragments(areas.transformed) %>%
  left_join(IR.Table, by = "Precursor.Ion.Name") %>%
  mutate(IR.Flag = ifelse((IR.Ratio < (IR.min - IR.flex) & IR.Ratio > (IR.max + IR.flex)), "IR.Flag", NA)) %>%
  select(Replicate.Name:Area, Quan.Trace:Second.Trace, IR.Flag)

# Retention Time Flags  ---------------------------------------
# If the Retention Time is "RT.flex" further away from the RT.Reference 
# Range from the RT.Range Table, add a flag. 
RT.flags.added <- IR.flags.added %>%
  merge(y = all.samples, all.x = TRUE) %>%
  left_join(RT.Range.Table, by = "Precursor.Ion.Name") %>%
  mutate(RT.Flag = ifelse((Retention.Time >= (RT.max + RT.flex) | Retention.Time <= (RT.min - RT.flex)), "RT.Flag", NA)) %>%
  select(Replicate.Name:Area, Quan.Trace:Second.Trace, IR.Flag, RT.Flag) %>%
  arrange(Precursor.Ion.Name, Product.Mz)

# Blank Flags  ---------------------------------------
# If the Area divided by the Blank.Reference value is
# greater than the set blk.thresh value, add a flag.
Blank.flags.added <- RT.flags.added %>%
  left_join(select(Blank.Table, Blank.max, Precursor.Ion.Name), by = "Precursor.Ion.Name") %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(Blank.Reference = Area * blk.thresh) %>%
  #mutate(blank.Flag = ifelse(((Protein.Name != "Internal Std") & (Area * blk.thresh) < Blank.max), "blank.Flag", NA)) %>%
  mutate(blank.Flag = ifelse(((Protein.Name != "Internal Std") & (Area * blk.thresh) < Blank.max), 
                             "blank.Flag", 
                             ifelse(((Protein.Name == "Internal Std") & (Area * blk.thresh < Blank.max)), "IS.blank.Flag", NA))) %>%
  select(Replicate.Name:RT.Flag, blank.Flag)


# Height Flags  ---------------------------------------
# Add a height.min.flag if the Height falls below the min.height
# value. Add an overloaded flag if the Height falls above the
# max.height value.
Height.flags.added <- Blank.flags.added %>%
  left_join(Height.Table, by = c("Replicate.Name", "Precursor.Ion.Name", "Precursor.Mz", "Product.Mz")) %>%
  mutate(height.min.Flag = ifelse((Height < min.height), "height.min.Flag", NA)) %>%
  mutate(overloaded.Flag = ifelse((Height > max.height), "overloaded.Flag", NA)) %>%
  select(-Height)

# Area Flags  ---------------------------------------
# If the Area is less than the area.min value, add a flag.
Area.flags.added <- Height.flags.added %>%
  mutate(area.min.Flag = ifelse((Area < area.min), "area.min.Flag", NA)) %>%
  select(-Sample.Type)

# Signal to Noise Flags  ---------------------------------------
# If the Signal to Noise ratio is less than the SN.min, add a flag.
SN.flags.added <- Area.flags.added %>%
  left_join(SN.Table, by = c("Replicate.Name", "Precursor.Ion.Name", "Area")) %>%
  mutate(SN.Flag = ifelse((Signal.to.Noise < SN.min), "SN.Flag", NA)) %>%
  select(-c(Sample.Type:Signal.to.Noise))

# All Flags  ---------------------------------------
# Add a column with all flags from the previous steps. 
final.table <- SN.flags.added %>%
  mutate(all.Flags = paste(IR.Flag, RT.Flag, blank.Flag, height.min.Flag, overloaded.Flag, area.min.Flag, SN.Flag, sep = ", ")) %>%
  mutate(all.Flags = as.character(all.Flags %>% str_remove_all("NA, ") %>%  str_remove_all("NA")))

# Remove Secondary trace ---------------------------------------
# Filter rows where Second.Trace == TRUE, keeping only Quan.Trace.
# Remove columns once finished.
final.table <- final.table %>%
  filter(Quan.Trace == TRUE) %>%
  select(-Quan.Trace, -Second.Trace)


# Standards & blank addition  ---------------------------------------
# Test for standars and blanks in the run. Add those standards
# and blanks back into the final table.

Stds.test <- grepl("_Std_", areas.raw$Replicate.Name)
Blks.test <- grepl("_Blk_", areas.raw$Replicate.Name)

if (any(Stds.test == TRUE)) {
  print("There are standards in this run. Joining to the bottom of the dataset!", quote = FALSE)
  ##
  standards <- areas.transformed %>%
    filter(Sample.Type == "std") %>%
    merge(y = master,
          by.x = c("Precursor.Ion.Name", "Product.Mz"),
          by.y = c("Compound.Name", "Daughter"),
          all.x = TRUE) %>%
    mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
    mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
    filter(Quan.Trace == TRUE) %>%
    select(Replicate.Name:Sample.Type)
  final.table <- rbind.fill(final.table, standards)
} else {
  print("No standards exist in this set.")
}

if (any(Blks.test == TRUE)) {
  print("There are blanks in this run. Joining to the bottom of the dataset!", quote = FALSE)
  ##
  blanks <- areas.transformed %>%
    filter(Sample.Type == "blk") %>%
    merge(y = master,
          by.x = c("Precursor.Ion.Name", "Product.Mz"),
          by.y = c("Compound.Name", "Daughter"),
          all.x = TRUE) %>%
    mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
    mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
    filter(Quan.Trace == TRUE) %>%
    select(Replicate.Name:Sample.Type)
  final.table <- rbind.fill(final.table, blanks)
} else {
  print("No blanks exist in this set.")
}

# Rename and save  ---------------------------------------
# Add comments restating the given QC parameters. Save to 
# current working directory with a new name, 
# "TQSQC_<original file name>.csv
con <- file(paste("TQSQC_", basename(input.file), sep = ""), open = "wt")
writeLines(paste("Hello! Welcome to the world of TQS Quality Control! ",
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








