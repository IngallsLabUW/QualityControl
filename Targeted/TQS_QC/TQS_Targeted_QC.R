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
area.min <- 1000
RT.flex <- 0.4
IR.flex <- 0.3
blk.thresh <- 0.3
SN.min <- 4

## ID run types
run.type <- tolower(str_extract(areas.raw$Replicate.Name, "(?<=_)[^_]+(?=_)"))

## Change variables classes as needed, add Sample.Type column
# TODO (rlionheart): drop index columm (why does it show up here but not in the other script?)
TransformVariables <- function(skyline.output) {
     # TODO (rlionheart): add documentation
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


# Retention Time  ----------------------------------------------------

RT.Range.Table <- areas.transformed %>%
        select(Precursor.Ion.Name, Retention.Time, Sample.Type) %>%
        group_by(Sample.Type) %>%
        filter(Sample.Type == "std") %>%
        group_by(Precursor.Ion.Name) %>%
        mutate(RT.Reference = mean(Retention.Time)) %>%
        mutate(RT.min = min(Retention.Time)) %>%
        mutate(RT.max = max(Retention.Time)) %>%
        mutate(RT.Flag = ifelse((abs(Retention.Time - RT.Reference) > RT.flex), "RT.Flag", NA)) 


## Height  ---------------------------------------
Height.Table <- areas.transformed %>%
     select(Precursor.Ion.Name, Height) %>%
     mutate(height.min.Flag = ifelse((Height < min.height), "height.min.Flag", NA)) %>%
     mutate(overloaded.Flag = ifelse((Height > max.height), "overloaded.Flag", NA))

## Area  ---------------------------------------
Area.Table <- areas.transformed %>%
     select(Precursor.Ion.Name, Area) %>%
     mutate(area.min.Flag = ifelse((Area < area.min), "area.min.Flag", NA))


## Signal to Noise  ---------------------------------------
SN.Table <- areas.transformed %>%
     select(Precursor.Ion.Name, Area, Background) %>%
     mutate(SN.Flag = ifelse(((Area / Background) < SN.min), "SN.Flag", NA)) 


# Ion Ratio Ranges Table  ----------------------------------------------------
# for each PIN, check if there is more than one fragment. If yes, determine which is quant and which is sec trace.

compound.list <- areas.raw %>%
     select(Precursor.Ion.Name, Precursor.Mz, Product.Mz, Sample.Type) %>%
     group_by(Precursor.Ion.Name) %>%
     filter(Sample.Type == "std") 

compounds <- as.character(unique(compound.list$Precursor.Ion.Name))
fragments <- unique(compound.list$Product.Mz)

for (i in 1:length(compounds)) {
     if (length(fragments) > 1) {
          print("woohoo")
     } else {
          print("nada")
     }
}


#CheckFragments function here





##########################################
cmpd.std.list<- split(areas.split[["std"]],
                      areas.split[["std"]]$Precursor.Ion.Name)
cmpds <- names(cmpd.std.list)
IR.range <- c()

for(i in 1:length(cmpds)){
     frags <- unique(cmpd.std.list[[i]]$Product.Mz)
     runs <- unique(cmpd.std.list[[i]]$Replicate.Name)
     ratios <-vector()
     
     ## is there more than one fragment?
     if (length(frags)>1){
          ## if yes then figure out which is the quantification trace
          ## and which is the confirmation (sec) trace
          mset<-master[master$Compound.Name==cmpds[[i]],]
          q.trace   <- mset$Daughter[mset$Quan.Trace=="yes"]
          sec.trace <- mset$Daughter[mset$X2nd.trace=="yes"]
          
          for(j in 1:length(runs)){
               ## is trace 2 > 5% of trace  1?
               sset<-cmpd.std.list[[i]][cmpd.std.list[[i]]$Replicate.Name==runs[j],]
               ratio <- sset$Area[sset$Product.Mz == q.trace] / sset$Area[sset$Product.Mz == sec.trace]
               
               ## if yes: get ion ratio for QC
               if(!is.na(ratio) & ratio < 100){
                    ratios <- c(ratios,ratio)
               } else {
                    ## if no: ignore ion ratio
                    ratios <- c(ratios, NA)
               }
          }
     } else {
          ## if no: ignore ion ratio
          ratios <- c(ratios,rep(NA,length(runs)))
     }
     IR.range <- cbind(IR.range,range(ratios,na.rm = T))
}

colnames(IR.range) <- cmpds
###########################################




