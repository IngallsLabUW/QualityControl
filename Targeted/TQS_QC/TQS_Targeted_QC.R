## First draft of a new TQS Targeted QC
# June 24th, 2019
# Regina Lionheart

library(reshape2)
library(tidyverse)

# Import files
# TODO (rlionheart): Decide on best import practice. Upload csv via shiny app (https://shiny.rstudio.com/gallery/file-upload.html)?
# TODO (rlionheart): Maybe rename the X2nd trace column for clarity.
areas.raw <- read.csv("./Targeted/TQS_QC/datafiles/ExampleSkylineOutput_TQS.csv", header = TRUE)
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


# Ion Ratio Ranges Table  ----------------------------------------------------
# CheckFragments function
# for each PIN, check if there is more than one fragment. If yes, determine which is quant and which is second trace.
# Is Second.Trace > 5% of Quant.Trace?

fragment.check <- areas.transformed %>%
        filter(Sample.Type == "std") %>%
        select(Replicate.Name, Precursor.Ion.Name, Precursor.Mz, Product.Mz, Sample.Type)

fragment.unique <-unique(fragment.check %>% select(Precursor.Ion.Name, Precursor.Mz, Product.Mz))

fragment.multi.unique <- fragment.unique %>%
        count(Precursor.Ion.Name)  %>%
        mutate(Two.Fragments = ifelse((n==1), FALSE, TRUE)) %>%
        select(-n)

fragment.check <- fragment.check %>%
        left_join(fragment.multi.unique) 

fragment.check <- fragment.check %>%
        merge(y = master,
              by.x = c("Precursor.Ion.Name", "Product.Mz"),
              by.y = c("Compound.Name", "Daughter"),
              all.x = TRUE) %>%
        select(Precursor.Ion.Name, Precursor.Mz, Product.Mz, Two.Fragments, Quan.Trace, Second.Trace) %>%
        mutate(Second.Trace = if_else(Second.Trace == "", "no" , as.character(Second.Trace))) %>%
        mutate(QT.Five.Percent = ifelse(Quan.Trace == "yes", 0.05 * Product.Mz, NA)) %>%
        # Everthing good above this line
        
        group_by(Precursor.Ion.Name) %>%
        mutate(Significant.Size = QT.Five.Percent > (Product.Mz & Second.Trace == "yes")) %>%
        mutate(IR.Ratio = ifelse(Significant.Size == TRUE, (Product.Mz[Quan.Trace == "yes"] / Product.Mz[Second.Trace == "yes"]), NA))

######################################        
# THIS WORKS SAVE IT
fake.data.frame <- areas.transformed %>%
        filter(Sample.Type == "std") %>%
        filter(Precursor.Ion.Name %in% c("Picolinic Acid", "Ketoglutaric Acid")) %>%
        select(Replicate.Name, Precursor.Ion.Name, Precursor.Mz, Product.Mz, Sample.Type) 
fake.data.frame$Product.Mz[which(fake.data.frame$Product.Mz == 57.0843)] = 101.0140


fake.data.unqiue <- unique(fake.data.frame %>% select(Precursor.Ion.Name, Precursor.Mz, Product.Mz))

fake.data.fragment.ismulti <- fake.data.unqiue %>%
        count(Precursor.Ion.Name)  %>%
        mutate(fragments = ifelse((n==1), FALSE, TRUE)) %>%
        select(-n)

fake.data.frame2 <- fake.data.frame %>%
        left_join(fake.data.fragment.ismulti)
######################################        


## if yes: get ion ratio for QC
if(!is.na(ratio) & ratio < 100){
        ratios <- c(ratios,ratio)
} else {
        ## if no: ignore ion ratio
        ratios <- c(ratios, NA)
}
######################################3
## Retention Time  ----------------------------------------------------
# TODO (rlionheart): Is the filtering of the "std" imperative? does it affect the Blank and RT reference??
RT.Range.Table <- areas.transformed %>%
        select(Precursor.Ion.Name, Retention.Time, Sample.Type) %>%
        #group_by(Sample.Type) %>%
        #filter(Sample.Type == "std") %>%
        group_by(Precursor.Ion.Name) %>%
        mutate(RT.Reference = mean(Retention.Time)) %>%
        mutate(RT.min = min(Retention.Time)) %>%
        mutate(RT.max = max(Retention.Time)) %>%
        mutate(RT.Flag = ifelse((abs(Retention.Time - RT.Reference) > RT.flex), "RT.Flag", NA)) 

RT.Range.Table <- as.data.frame(RT.Range.Table)

## Blanks ---------------------------------------
Blank.Table <- areas.transformed %>%
        select(Precursor.Ion.Name, Area, Sample.Type) %>%
        group_by(Sample.Type) %>%
        # filter(Sample.Type == "blk") %>%
        group_by(Precursor.Ion.Name) %>%
        mutate(Blank.Reference = mean(Area, na.rm = TRUE)) %>%
        mutate(blank.Flag = ifelse((Area / Blank.Reference) < blk.thresh, "blank.Flag", NA))

        
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



first.joins <- areas.transformed %>%
        select(Replicate.Name:Background, Height) %>%
        left_join(Height.Table) %>%
        left_join(Area.Table) %>%
        left_join(SN.Table)


        



##########################################
areas.split <- split(areas.transformed, areas.transformed$Sample.Type)

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
          sec.trace <- mset$Daughter[mset$Second.Trace=="yes"]
          
          for(j in 1:length(runs)){
               ## is trace 2 > 5% of trace 1?
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
     IR.range <- cbind(IR.range, range(ratios, na.rm = T))
}

colnames(IR.range) <- cmpds
###########################################




