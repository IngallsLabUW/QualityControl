## Quality Control for TQS Skyline output

## Note that the csv can have blanks but not characters in 
## the RT, MZ, and Area columns

## Load needed libraries
library("plyr")
library("reshape2")

## Set the parameters for the QC ----------------------------------------
filename <- "ExampleSkylineOutput.csv"
areas.raw <- read.csv(filename ,header=TRUE, comment.char="", as.is=TRUE)
masterfile <- "HILIC_MasterList_Example.csv"
master <- read.csv(masterfile)

## pick overload value (suggestion: 1e8)
max.height <- 1.0e8
## pick the minimum height to be counted as a 'real' peak
## (suggestions: HILIC - 1000, cyano - 100)
min.height <- 1000
## pick RT flexibility 
## (suggeestions: +/- 0.4 min for HILIC, +/- 0.2 min for cyano)
RT.flex <- 0.4
## pick IR flexibility (suggestion: +/- 30%)
IR.flex <- 0.3
## pick how big of a signal the samples need to have in comparison
## to a blank in order to 'count' (suggestion: +/- 20%)
blk.thresh <- 0.3
## pick the value that will be the acceptable signal to noise ratio
## (suggestion maybe 5 for cyano and 4 for HILIC?  (broader peaks make for more background))
SN.thresh <- 4

## ID run types ---------------------------
## Standards (std), Samples (smp), Blanks (blk), Pooled (poo)
t1<-strsplit(areas.raw$Replicate.Name, "[_]")
type <- vector()
for(i in 1:length(t1)){
     this.one = t1[[i]][2]
     type <- c(type, this.one)  
}
type<- tolower(type)

areas.raw$sample.type <- as.factor(type)
areas.raw$Precursor.Ion.Name <- as.factor(areas.raw$Precursor.Ion.Name)

areas.raw$Area <- as.numeric(areas.raw$Area)
areas.raw$Retention.Time <- as.numeric(areas.raw$Retention.Time)
areas.raw$Background <- as.numeric(areas.raw$Background)
areas.raw$Height <- as.numeric(areas.raw$Height)

stdRows <- areas.raw$sample.type == "std"
blkRows <- areas.raw$sample.type == "blk"

areas.split<-split(areas.raw,areas.raw$sample.type)

type.options <- names(areas.split)

cmpd.blk.list<- split(areas.split[["blk"]],
                      areas.split[["blk"]]$Precursor.Ion.Name)

## Check the range of RT and ion ratio in Standards ---------------------
## Range of RTs
RT.range<-sapply(split(areas.split[["std"]]$Retention.Time,areas.split[["std"]]$Precursor.Ion.Name),
                 range,na.rm = T)

## Range of ion ratios
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
          q.trace <- mset$Daughter[mset$Quan.Trace=="yes"]
          sec.trace <- mset$Daughter[mset$X2nd.trace=="yes"]
          
          for(j in 1:length(runs)){
               ## is trace 2 > 5% of trace  1?
               sset<-cmpd.std.list[[i]][cmpd.std.list[[i]]$Replicate.Name==runs[j],]
               ratio <- sset$Area[sset$Product.Mz==q.trace] /
                    sset$Area[sset$Product.Mz==sec.trace]
        
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

## Discard/flag data that is outside these ranges ---------------------- 
## Set up dataframes for getting the relevant info out:
samp.data <- areas.split[["smp"]]
## if there are pooled samples, include them in the 'sample' list
if(any(type.options=="poo")){
     poo.data <- areas.split[["poo"]]
     samp.data <- rbind(samp.data,poo.data)
}
cmpds <- unique(samp.data$Precursor.Ion.Name)
samples <- unique(samp.data$Replicate.Name)

cmpd.samp.dfs <- split(samp.data,samp.data$Precursor.Ion.Name)

## RT range +/- 0.3 min = RT.ok ------------------------------
RT.ok <- RT.range
RT.ok[1,] <- RT.range[1,] - RT.flex
RT.ok[2,] <- RT.range[2,] + RT.flex

# Get the sample RTs
RT.matrix <- c()
for (i in 1:length(cmpds)){
     RT <- sapply(split(cmpd.samp.dfs[[cmpds[i]]]$Retention.Time,
                        cmpd.samp.dfs[[cmpds[i]]]$Replicate.Name),
                  mean, na.rm = T)
     # if they are outside the right range, replace them with NA
     compound.class <- master$Group[master$Compound.Name==as.character(cmpds[i])][1]
     if(compound.class!="Internal Std"){
          for (j in 1:length(RT)){
               if (!is.na(RT[j]) & RT[j]<RT.ok[,cmpds[i]][1]) {
                    RT[j] <- NA
               } else if (!is.na(RT[j]) & RT[j]>RT.ok[,cmpds[i]][2]) {
                    RT[j] <- NA
               }
          }
     }
     # save into RT.matrix
     RT.matrix <- cbind(RT.matrix, RT)
     colnames(RT.matrix)[ncol(RT.matrix)]<-as.character(cmpds[i])
}

#convert short format matrix to long form ----------------------
output <- melt(RT.matrix, value.name = "Retention Time")
colnames(output) <- c("Replicate.Name","Compound.Name","Retention.Time")

# add note if the area was thrown out because of a bad RT
output$Notes <- rep("", nrow(output))
for(i in 1:nrow(output)){
     if(is.na(output$Retention.Time[i])){
          output$Notes[i] <- "Bad RT"
     }
}

## Add areas, heights, and Ion Ratios to output ----------------------
Area <- c()
Pk.Height <- c()
IR <- c()
blank.data <- c()
S.N <- c()
i=1
for (i in 1:length(cmpds)){
     ## find the right trace for that compound
     mset<-master[master$Compound.Name==as.character(cmpds[i]),]
     q.trace <- mset$Daughter[mset$Quan.Trace=="yes"]
     key <- as.character(cmpds[i])
     ssets <- split(cmpd.samp.dfs[[key]],cmpd.samp.dfs[[key]]$Product.Mz)
     areas.quant <- ssets[[as.character(q.trace)]]
     if (length(grep(pattern = "yes",mset$X2nd.trace))>0){
          sec.trace <- mset$Daughter[mset$X2nd.trace=="yes"]
          areas.sec <- ssets[[as.character(sec.trace)]]
          key.info <- merge(areas.quant, areas.sec, 
                       by = c("Precursor.Ion.Name","Replicate.Name",
                              "Precursor.Mz","sample.type"),
                       all = TRUE)
          key.IR <- key.info$Area.x/key.info$Area.y
     } else {
          key.IR <- rep(NA,length(areas.quant$Area))
     }
     # print(paste(cmpds[i],nrow(areas.quant)))
     Area <- c(Area, areas.quant$Area[order(areas.quant$Replicate.Name)])
     
     Pk.Height <- c(Pk.Height, areas.quant$Height[order(areas.quant$Replicate.Name)])
     IR <- c(IR, key.IR)
     S.N <- c(S.N, (areas.quant$Area[order(areas.quant$Replicate.Name)]+
                         areas.quant$Background[order(areas.quant$Replicate.Name)])/
                   areas.quant$Background[order(areas.quant$Replicate.Name)])
     ## Use the q.trace to get the blank data we should use
     bsets <- split(cmpd.blk.list[[key]],cmpd.blk.list[[key]]$Product.Mz)
     ## Find the blank with the highest peak area
     blk.max <- max(bsets[[as.character(q.trace)]]$Area,na.rm = TRUE)
     ## Save out the info for this run and this transition so we can put
     ## it at the end of the output later on
     blank.data <- rbind(blank.data, 
                         bsets[[as.character(q.trace)]][which(bsets[[as.character(q.trace)]]$Area==blk.max),])

}
output$Area <- Area
output$Height <- Pk.Height
output$IR <- IR
output$S.N <- S.N

## Add a column for the unmodified area - we won't touch this during QC
output$rawArea <- output$Area

## if the retention time is bad then delete the area
for (i in 1:nrow(output)){
     if (output$Notes[i]!=""){
          output$Area[i] <- NA
     }
}

## Is it overloaded? ---------------------------------------
## Check peak height
## If it is > specified value (max.height) then keep the area data
## but make a note that it may be overloaded
for (i in 1:nrow(output)){
     if (!is.na(output$Height[i]) & output$Height[i]>max.height){
          note = "overloaded?"
          output$Notes[i]<-paste(output$Notes[i],note,sep="")
     }
}


## Ion ratio range +/- 20% -------------------------------------
IR.ok <- IR.range
IR.ok[1,] <- IR.range[1,]*(1-IR.flex)
IR.ok[2,] <- IR.range[2,]*(1+IR.flex)

## Check the sample Ion Ratios

for (i in 1:nrow(output)){
     cmpd.name <- as.character(unique(output$Compound.Name[i]))
     this.range <- IR.ok[,cmpd.name]
     compound.class <- master$Group[master$Compound.Name==cmpd.name][1]
     
     if (is.na(this.range[1]) | is.infinite(this.range[1])) {
          note <- ""
     } else if (compound.class != "Internal Std") {
          if (is.na(output$IR[i])){
               note <- "No second trace"
          } else if (output$IR[i] < this.range[1] | 
                     output$IR[i] > this.range[2] |
                     is.infinite(output$IR[i])){
               note <- "Bad IR"
          } else {
               note <- ""
          }
     }
     output$Notes[i] <- paste(output$Notes[i], note, sep="") 
}


## Absolute height minimum > 100 (of larger transition/quantification trace) ------
for (i in 1:nrow(output)){
     if (!is.na(output$Height[i]) & output$Height[i]<min.height){
          output$Area[i] <- NA
          output$Notes[i] <- paste(output$Notes[i],"too small",sep="")
     }
}


## Blank Correct -----------------------------
## If the peak area of the sample isn't big enough in comparison
## to that of the largest blank, then get rid of that data

for(j in 1:nrow(output)){
     key <- as.character(output$Compound.Name[j])
     blk.max <- blank.data$Area[blank.data$Precursor.Ion.Name==key]
     compound.class <- master$Group[master$Compound.Name==key][1]
     
     ## if the area of the blank is more than 20% of the
     ## area in the sample, make a note
     if(length(blk.max)!=1){
     } else if (!is.na(output$Area[j]) && 
                (output$Area[j]*blk.thresh)<blk.max &&
                compound.class!="Internal Std"){
          if (!grepl(pattern = "loaded",output$Notes[j])){
               output$Area[j] <- NA
          }
          output$Notes[j] <- paste(output$Notes[j],"comparable to blk",sep="")
     }   
}


## Signal to Noise check -----------------------------------------------
## If S/N is less than 5, throw it out!
for (i in 1:nrow(output)){
     if (!is.na(output$S.N[i]) & output$S.N[i]<SN.thresh){
          output$Notes[i] <- paste(output$Notes[i],"bad S/N",sep="")
          if (!grepl(pattern = "loaded",output$Notes[i])){
               output$Area[i] <- NA
          }
     }
}


## attach blank data to output ---------------------------
blank.data$Compound.Name <- blank.data$Precursor.Ion.Name
blank.data$Notes <- rep("Blank used for comparison",nrow(blank.data))
blank.data$IR <- rep(NA,nrow(blank.data))
blank.data$S.N <- (blank.data$Area+blank.data$Background)/blank.data$Background
blank.data$rawArea <- blank.data$Area

final.output <- rbind(output, blank.data[,colnames(output)])

## Output with: comment -------------------------
## Ion name, area, was a peak removed?
comment.text <- paste("# Hello! welcome to your data! ","Overload height: ", 
                      max.height, ". ", "RT flexibility: ", RT.flex, ". ",
                      "Ion ratio flexibility: ", IR.flex, ". ",
                      "Blank can be this fraction of a sample: ",blk.thresh, ". ", 
                      "S/N threshold: " , SN.thresh, ". ",
                      "Minimum peak height: ", min.height, ". ",
                      "Processed on: ", Sys.time(), sep="")
new.filename <- paste("QC_output",filename,sep="")
con <- file(new.filename, open="wt")
writeLines(paste(comment.text), con)
write.csv(final.output, con)
close(con)

