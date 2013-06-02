normWrapper <- function(galFile,files=dir(pattern="*\\.gpr$"),fileType=".gpr",BGsub.methods="edwards",norm.Methods="loess",scale=FALSE,prenormFilterFun="empty",
                        filterFlaggedSpots=TRUE,filterPercentPresent=FALSE,percentPresentThreshold=0.9,filterRGR=FALSE,filterControls=FALSE,RGRthreshold=c(0.5,2),
                        PDFName="MAD.SD.pdf",verbose=TRUE,hSQindex=5,nameIndex=4,oligoIndex=11,duplicatesFile="duplicates.txt",boxplotname="boxplot.pdf"
){
  #FUNCTION: normWrapper
  #calls all other functions in sequence
  #ARGUMENTS: 
  #     galFile : (string) the galfile used to extract gene names and IDs. 
  #     fileType : (string) “.srr ” or “.gpr”. What file should be used for filtering?
  #     BGsub.methods: (string vector) which background subtraction methods should be used? C
  #     Norm.Methods (string vector) which normalization methods should be used? 
  #     scale (logical) should scale between-array normalization be used? h
  #     PrenormFilterFun (string) What function should be used to assign weights to spots before normal- ization? Choices are:
  #         ”none” : Assigns all spots weight 1.
  #         ”empty”: assigns spots with name ”EMPTY” weight 0, all others 1.
  #         "flagempty” : assigns flagged and ”EMPTY” spots weight 1, alll others 0.
  #         ”control”: assigns control spots and ”EMPTY” spots weight 0, all others 1. 
  #                   NOTE! to use this option, you must provide a .gal file with oligoIDs in column oligoIndex
  #      filterFlaggedSpots (logical) Should flagged spots be filtered out after normalization?
  #      filter PercentPresent(logical)Should spots with fewe rthan percentPresentThreshold spots present 
  #                   across arrays be filtered out after normalization?
  #      percentPresentThreshold (number between 0 and 1) the lower bound of percent of spots present that is allowed in the filtered data.
  #      removeRGR (logical) should spots with red-green ratios between the parameters of RGRthreshold be filtered out?
  #      RGRthreshold (2-element numerical vector) Spots with red-green ratio between RGRthreshold[1] and RGRthreshold[2] are filtered out. 
  #      PDFName (string) The name of the pdf where MAD/SD plots will be deposited.
  #      verbose (logical) should progress reports be printed to the screen?
  #      hSQindex (numerical) column in the .gal file containing SUIDs.
  #      nameIndex (numerical) column in the .gal file containing full gene names.
  #      oligoIndex (numerical) column in the .gal file containing oligoIDs.
  #      duplicatesFile (string) a file contianing a logical vector of which rows in these arrays are 
  #            duplicates. This file makes AvgDupRows run much faster. It will be generated the first time you run 
  #            AvgDupRows in this directory. It is also sufficient to use a duplicate file generated from a run of 
  #            AvgDupRows on a separate, identical print run.
  
  #convert .srr files to .gpr
  files <- SRRtoGPR(files=files,galFile = galFile,hSQindex=hSQindex,nameIndex=nameIndex)
  
  #normalize the .gpr files; store normalized filenames in NormalizedFiles
  NormalizedFiles <- BGsubandNormalize(files = files, BGsub.methods=BGsub.methods,scale=scale,norm.Methods=norm.Methods,
                                       galFile=galFile,verbose=verbose,filterFun=prenormFilterFun,
                                       oligoIndex=oligoIndex,nameIndex=nameIndex,hSQindex=hSQindex) 
  
  #filter normalized data; store filtered filenames in FilteredFiles
  FilteredFiles <- filter(filterFiles=NormalizedFiles,galFile=galFile,fileType=fileType,oligoIndex=oligoIndex,
                          filterControls=filterControls,filterFlaggedSpots=filterFlaggedSpots,filterPercentPresent=filterPercentPresent,
                          percentPresentThreshold=percentPresentThreshold,filterRGR=filterRGR,RGRthreshold=RGRthreshold,verbose=verbose)
  
  # generate MAD vs. SD plots of filtered data
  ReplicateDiffPlot(files=FilteredFiles,PDFName=PDFName)
  
  # average duplicate rows in filtered data; store averaged filenames in averagedFiles
  averagedFiles <- AvgDupRows(files=FilteredFiles,dupFile=duplicatesFile)
  
  #generate boxplots of averaged data
  ArrayBoxPlot(averagedFiles,boxplotname)
  
browser()
}


SRRtoGPR <- function(files=dir(pattern="*\\.srr$"),galFile,hSQindex=5,nameIndex=4){
  #FUNCTION: SRRTOGPR
  #       Converts .srr files to .gpr files and deposits them in the working directory.
  # ARGUMENTS:
  #       files (string vector) : the file names to be converted.
  #       galFile : the gal file that will be used to fill in gene names and SUIDs in the output. 
  #       hSQIndex (numerical) column in the .gal file containing SUIDs.
  #       NameIndex (numerical) : column in the .gal file containing full gene names.
  
  #if there are already .gpr files in the current directory corresponding to each , return them
  gprexists <- sapply(files,function(srr){
    file.exists(paste(substr(srr,start=0,stop=nchar(srr) - 4),".gpr",sep=""))
  })
  if (sum(gprexists) == length(gprexists)) return(dir(pattern="*\\.gpr$"))
  
  #error checking
  if (!file.exists(galFile)) stop("SRRtoGPR: File ",galFile," not found \n")
  if (length(Files) == 0) stop("SRRtoGPR: no srr files found\n")
  
  genenames <- read.delim(galFile,header=TRUE, skip = 55)
  
  #dataHeader: column headers for the data section of the new .gpr file. 
  dataHeader <- c("Block",  "Column",  "Row",  "Name",  "ID",  "X",  "Y",  "Dia.",  "F635 Median","F635 Mean",
                  "F635 SD",  "B635",  "B635 Median",   "B635 Mean",	"B635 SD", "% > B635+1SD",	"% > B635+2SD",	
                  "F635 % Sat.","F532 Median",	"F532 Mean",	"F532 SD",	"B532",	"B532 Median", "B532 Mean",	
                  "B532 SD",	"% > B532+1SD",	"% > B532+2SD",	"F532 % Sat.","Ratio of Medians (635/532)",	
                  "Ratio of Means (635/532)",	"Median of Ratios (635/532)",	"Mean of Ratios (635/532)",	"Ratios SD (635/532)",
                  "Rgn Ratio (635/532)",	"Rgn R? (635/532)",	"F Pixels",	"B Pixels",	"Circularity","Sum of Medians (635/532)",	
                  "Sum of Means (635/532)",	"Log Ratio (635/532)",	"F635 Median - B635",	"F532 Median - B532",
                  "F635 Mean - B635",	"F532 Mean - B532",	"F635 Total Intensity", "F532 Total Intensity",	
                  "SNR 635",	"SNR 532",	"Flags",	"Normalize",	"AutoFlag")
  outfiles <- vector()
  for (i in 1:length(files)){
    cat("converting ",files[i]," to .gpr\n")
    #construct a .gpr name from the .srr name.
    outputfilename <- paste(c(substr(files[i], start = 1, stop = nchar(files[i]) - 3),"gpr"), sep = "", collapse = "")    
    #header: stores the first 32 rows of the file
    header <- as.matrix(read.delim(files[i], header=FALSE))[1:32,1:2]
    #retreive data from the .srr files.
    data <- as.matrix(read.delim(files[i], header=FALSE, comment.char="", skip = skip))
    #output: stores all information (header + data) that will go in the .gpr file
    output <- matrix(nrow = (32 + nrow(data)), ncol = 52, data = "")
    #copy the header into output
    output[1:32,1:2] <- header
    #copy the data section header into line 33 of output
    output[33,] <- dataHeader
    #convert all relevant columns in .srr to corresponding columns in .gpr
    output[34:nrow(output),1] <- data[2:nrow(data),3]
    output[34:nrow(output),2] <- data[2:nrow(data),5]
    output[34:nrow(output),3] <- data[2:nrow(data),4]
    output[34:nrow(output),4] <- as.character(genenames[,nameIndex])
    output[34:nrow(output),5] <- as.character(genenames[,hSQindex])
    output[34:nrow(output),6:7] <- data[2:nrow(data),10:11]
    output[34:nrow(output),8] <- data[2:nrow(data),14]
    output[34:nrow(output),9] <- data[2:nrow(data),19]
    output[34:nrow(output),10] <- data[2:nrow(data),21]
    output[34:nrow(output),11] <- data[2:nrow(data),23]
    output[34:nrow(output),12] <- data[2:nrow(data),43]
    output[34:nrow(output),13] <- data[2:nrow(data),43]
    output[34:nrow(output),14] <- data[2:nrow(data),45]
    output[34:nrow(output),15] <- data[2:nrow(data),47]
    output[34:nrow(output),16] <- data[2:nrow(data),31]
    output[34:nrow(output),17] <- data[2:nrow(data),33]
    output[34:nrow(output),18] <- data[2:nrow(data),27]
    output[34:nrow(output),19] <- data[2:nrow(data),18]
    output[34:nrow(output),20] <- data[2:nrow(data),20]
    output[34:nrow(output),22] <- data[2:nrow(data),42] 
    output[34:nrow(output),23] <- data[2:nrow(data),42]
    output[34:nrow(output),24] <- data[2:nrow(data),44]
    output[34:nrow(output),25] <- data[2:nrow(data),46]
    output[34:nrow(output),26] <- data[2:nrow(data),30]
    output[34:nrow(output),27] <- data[2:nrow(data),32]
    output[34:nrow(output),28] <- data[2:nrow(data),26]
    output[34:nrow(output),29] <- data[2:nrow(data),57]
    output[34:nrow(output),30] <- data[2:nrow(data),61]
    output[34:nrow(output),31] <- data[2:nrow(data),63]
    output[34:nrow(output),32] <- data[2:nrow(data),65]
    output[34:nrow(output),33] <- data[2:nrow(data),67]
    output[34:nrow(output),34] <- data[2:nrow(data),77]
    output[34:nrow(output),35] <- data[2:nrow(data),79]
    output[34:nrow(output),36:37] <- data[2:nrow(data),8:9]
    output[34:nrow(output),38] <- data[2:nrow(data),15]
    output[34:nrow(output),39] <- data[2:nrow(data),69]
    output[34:nrow(output),40] <- data[2:nrow(data),71]
    output[34:nrow(output),41] <- data[2:nrow(data),59]
    output[34:nrow(output),42] <- data[2:nrow(data),35]
    output[34:nrow(output),43] <- data[2:nrow(data),34]
    output[34:nrow(output),44] <- data[2:nrow(data),37]
    output[34:nrow(output),45] <- data[2:nrow(data),36]
    output[34:nrow(output),46] <- data[2:nrow(data),39]
    output[34:nrow(output),47] <- data[2:nrow(data),38]
    output[34:nrow(output),48] <- data[2:nrow(data),41]
    output[34:nrow(output),49] <- data[2:nrow(data),40]
    output[34:nrow(output),50] <- data[2:nrow(data),6]
    output[34:nrow(output),52] <- data[2:nrow(data),17]
    #update the ATF-format column count to reflect .gpr format 
    output[2,2] <- 52
    output[1,1] <-  "ATF"
    #write to output file
    write.table(output, file = outputfilename,sep = "\t", quote = FALSE,row.names = FALSE, col.names = FALSE)
    outfiles <- c(outfiles,outputfilename)
  }
  
  return(outputfilename)
}


BGsubandNormalize <- function(BGsub.methods="edwards",norm.Methods="loess",galFile,verbose=TRUE, files=dir(pattern="*\\.gpr$"),
                              scale=FALSE,filterFun="empty",oligoIndex=11,hSQindex=5,nameIndex=4){
  #ARGUMENTS:
  #     galfile : (string) the galfile used to extract gene names and IDs. 
  #     files (string vector) : which files should be normalized?
  #     BGsub.methods: (string vector) which background subtraction methods should be used? 
  #     norm.Methods (string vector) which normalization methods should be used? 
  #     scale (logical) should scale between-array normalization be used? h
  #     PrenormFilterFun (string) What function should be used to assign weights to spots before normal- ization? Choices are:
  #         ”none” : Assigns all spots weight 1.
  #         ”empty”: assigns spots with name ”EMPTY” weight 0, all others 1.
  #         "flagempty” : assigns flagged and ”EMPTY” spots weight 1, alll others 0.
  #         ”control”: assigns control spots and ”EMPTY” spots weight 0, all others 1. 
  #                  NOTE! to use this option, you must provide a .gal file with oligoIDs in column oligoIndex
  #      hSQIndex (numerical) column in the .gal file containing SUIDs.
  #      NameIndex (numerical) column in the .gal file containing full gene names.
  #      OligoIndex (numerical) column in the .gal file containing oligoIDs.
  require(limma)  
  gal <- read.delim(galFile,header=TRUE, skip = 55)
  
  # genenames : will become the first 2 columns of the output file; SUIDs and gene names.
  genenames <- cbind(as.character(gal[,hSQindex]),as.character(gal[,nameIndex]))
  if (ncol(gal) >= oligoIndex) oligoIDs <- gal[,oligoIndex]
  
  # RGraws : RGlist objects containing red and green fluorescence values for all arays
  if (!(filterFun %in% c("control","empty","flagempty","prefilter","none"))) 
    stop("Choices for filterFun are: control,empty,flagempty,prefilter,none. You provided ", filterFun)
  if (filterFun=="control"){
    RGraws <- read.maimages(files,"genepix",wt.fun=function(x){
      x <- cbind(x,oligoIDs)      
      x$Name <- genenames[2]
      x$ID <- genenames[1]
      ((substr(x$oligoIDs,start=2,stop=2) != "C") & (x$Name != "EMPTY"))
    })
  }else   RGraws <- read.maimages(files,"genepix",wt.fun=stringToFun(filterFun))   
  # If I intend to filter by controls, add the oligoIds into coulumn 53 of .gal files
  
  normFiles <- c()
  #for each background subtraction and normalization method, perform background-subtaction and normalization
  for (b in 1:length(BGsub.methods)){
    for (n in 1:length(norm.Methods)){ #for each normalization and background subtraction combination...
      if ((n != 1)|(b != 1)) rm(outputmatrix)
      if (verbose) cat(c("normalizing :", norm.Methods[n], ", ", BGsub.methods[b], "    scale:", scale, "\n"))
      outputFilename <- paste(c("Prenormfiltered-",filterFun,"_",norm.Methods[n], "_",BGsub.methods[b], ".txt"), sep = "", collapse = "")
      
      # normalize using Limma
      normalized <- normalizeWithinArrays(RGraws, layout = RGraws$printer, method=norm.Methods[n],  
                                          span=0.3, iterations=6, controlspots=NULL, df=5, bc.method=BGsub.methods[b], offset=50)
      # perform scale normalization if necessary
      if (scale) {
        normScaled <- normalizeBetweenArrays(normalized,method="scale",cyclic.method="affy")
        outputmatrix <- cbind(genenames,normScaled$M)
        outputFilename <- paste(c("scale_",outputFilename), sep = "", collapse = "")
      }else outputmatrix <- cbind(genenames,normalized$M)
      # construct output out of gene names and M (red-green ratios) values from Limma
      colnames(outputmatrix)[1:2] <- c("suid","name")
      
      #write output to file
      write.table(outputmatrix,file=outputFilename,append=FALSE,quote=FALSE,sep="\t", row.names = FALSE)
      normFiles <- c(normFiles,outputFilename)
    }
  } 
  return(normFiles)
}


stringToFun <- function(string){
  #stringToFun: convert a user-provided string token into a function for pre-normalization filtering.
  #to  create more pre-norm filtering options, write a funtion that reuturns 1 on spots you 
  # want, and 0 on spots you don't want. 
  # useful tip: FALSE = 0, TRUE = 1.
  
  #"none" : return 1 for all spots
  if (string == "none") return(function(x){return (1)})
  
  #empty: return 0 for spots with name "EMPTY", 1 otherwise.
  if (string == "empty") return (function(x){ as.numeric(x$Name != "EMPTY")})
  
  #"flagempty": return 0 for flagged spots and "EMPTY" spots, 1 otherwise
  if (string == "flagempty" || string == "prefilter") return (flagempty <- function(x){
    okflags <- (x$Flags == 0)
    okName <- (x$Name != "EMPTY")
    as.numeric(okflags)
  })
  
  #"control" : return 0 for control spots and "EMPTY" spots, 1 otherwise.
  # must have oligoIDs in column 53 of .gpr file for this to work.
  if (string == "control")  return(function(x,oligoIDs) ((substr(x[,53],start=2,stop=2) != "C") && (x$Name != "EMPTY"))) 
}


AvgDupRows <- function(files,dupFile="duplicates.txt"){
  #FUNCTION: AvgDupRows
  #      combines rows with identical SUIDs in to a single row containing the average values for all rows of that SUID.
  #ARGUMENTS: 
  #      files (string vector) The names of the files to be filtered. 
  #      duplicatesFile (string) a file contianing a logical vector of which rows in these arrays are duplicates.
  
  #make sure you have Hmisc for aggregate()
  if (!is.element("Hmisc", installed.packages()[,1])) install.packages("Hmisc")
  require(Hmisc)
  averagedFiles <- c()
  
  for (i in 1:length(files)){
    if (i == 1) {
      # find which rows are duplicates. If there is a valid duplicate file, use it; 
      #otherwise, sort data, iterate through all rows and find the duplicates.
      if (!is.null(dupFile) && file.exists(dupFile))  dupVec <- read.delim(dupFile)[,1]
      else{
        data <- read.delim(files[1])
        data <- data[order(data$suid),]
        dupVec <- lapply(1:nrow(data), function(x, vec){
          if ((x %% 500) == 0) cat("finding duplicates in row",x,"\n")
          return(length(unique(vec[(x - 1):(x + 1)][-2])) == length(unique(vec[(x - 1):(x + 1)])))
        }, data$suid)
        write.table(unlist(dupVec),file="duplicates.txt",quote=FALSE,row.names=FALSE)
      }
    }
    cat("\nAveraging duplicate rows in ",files[i])
    data <- read.delim(files[i])
    ordered <- data[order(data$suid),]
    # add an additional column to the data to tell which rows are duplicates
    ordered$dup <- dupVec
    
    #dupDF: dataframe containing only duplicate rows.
    dupDF <- ordered[ordered$dup == TRUE,1:ncol(ordered) - 1]
    
    #uniqueDF: dataframe containing only unique rows.
    uniqueDF <- ordered[ordered$dup == FALSE,1:ncol(ordered) - 1]
    
    #aggregate duplicate rows and find their average.
    suppressWarnings(agg <- aggregate(dupDF,by=list(dupDF$suid),FUN=function(x){
      if (!is.numeric(x)) return (x[1])
      return (mean(na.omit(x)))
    }))
    
    #combine the aggregated duplicates with the originally unique rows.
    agg <- rbind(agg[-1],uniqueDF)
    agg <- agg[order(agg$suid),]
    
    #write to output file.
    outFilename <- paste(c("Avg_",files[i]), sep = "", collapse = "")
    write.table(agg, outFilename, quote = FALSE, sep = "\t",row.names=FALSE)
    averagedFiles <- c(averagedFiles,outFilename)
  }
  return(averagedFiles)
}


ArrayBoxPlot <- function(files,outFilename="Boxplots.pdf"){
  #ARGUMENTS:
  #     files: (string vector) the files to be boxplotted. These should contain 
  #     data matrices of the form returned by BGsubandNormalize, filter, or AvgDupRows; 
  #     that is, they contain gene names and IDs in the first two columns, array names in the first row,
  #     and red-green ratios in the rest of the file.
  #     outFilename: (string) : the .pdf name where boxplots will be deposited.
  if (length(files) == 0) stop("ArrayBoxPlot: specify what files to boxplot")
  pdf(outFilename)
  
  #set parameters to 4 plots per page and smaller font size
  if (length(files) > 2) par(mfrow = c(2,2)) # plots per mage
  par(cex.main = 0.5) #font size 
  for (f in 1:length(files)){
    cat("plotting ", files[f],"\n")
    data <- read.delim(files[f])
    
    #convert data to a list because that's how you produce multiple plots
    datalist <- list()
    for (i in 3:(ncol(data))) datalist[[i - 2]] <- as.numeric(data[,i])
    
    #see ?ploxplot for information about these parameters
    boxplot(datalist, names = substring(tail(colnames(data),-2),5, 7), cex = .5, cex.axis = .7,
            las = 2,main = substr(files[f],start=1,stop=nchar(files[f]) - 3),xlab = "Arrays",
            ylab="log2 ratios", medcex = 0.5, outcex = 0.5)
  }
  dev.off()
}


ReplicateDiffPlot <- function(files,PDFName="MAD.SD.pdf"){
  # ARGUMENTS:
  #       files (string vector): the files to be plotted. These should contain data matrices 
  #       of the form returned by BGsubandNormalize, filter, or AvgDupRows; that is, 
  #       they contain gene names and IDs in the first two columns, array names in the first row, 
  #       and red-green ratios in the rest of the file.
  #       outFilename  (string) : the name of the .pdf where plots will be deposited.
  
  # construct the filename for infofile (number of duplicate spots, standard devation, and MAD)
  infoFileOutfilename <- (paste(c(substr(PDFName, start = 0, stop = nchar(PDFName) - 4),
                                  "info.txt"),sep="",collapse=""))
  colors <- c("black","blue","mediumpurple","green","green4","orange","red","pink","brown")[1:length(files)]
  #error checking
  for (i in 1:length(files)) if (!file.exists(files[i])) stop ("file ", files[i], " does not exist")
  
  pdf(PDFName)
  normData <- (read.delim(files[1],header=TRUE))
  
  #set up info file containing the correct number of rows, column names and row names
  infoFile <- matrix(data = "", ncol = 3*length(files) + 1, nrow = ncol(normData) - 1)
  infoFile[,1] <- c(tail(colnames(normData),-2), "average")
  colnames <- vector ("character",3*length(files))
  colnames[1] <- ""
  for (i in 1:length(files)){
    colnames[i*3 - 1] <- paste(c(files[i], " duplicate MAD"), sep = "", collapse = "")
    colnames[i*3] <- paste(c(files[i], " RGR stdevs"), sep = "",collapse = "")
    colnames[i*3 + 1] <- paste(c(files[i], " npairs"), sep = "", collapse = "")
  }
  colnames(infoFile) <- colnames
  
  for (j in 1:length(files)){ # for each file...
    cat("Plotting ", files[j], "\n")
    #get all normalized data from this file
    normData <- read.delim(files[j],header=TRUE)
    uniqueDiffs <- sort(unique(normData[,1]))
    # remove empty spots
    if (uniqueDiffs[1] == "EMPTY") uniqueDiffs <- tail(uniqueDiffs,-1)
    normData <- normData[normData[,1] != "EMPTY",]
    # for each arrray within this file, calculate the standard deviation of RGR and
    # populate uniqueDiffs with the pairwise differences between duplicate spots.
    stdevs <- allnspots <- medians <- vector()
    for (i in 3:ncol(normData)) { 
      #find the standard deviation of an array (a column in NormData)
      stdev <- sd(na.omit(normData[,i]))
      
      #can't plot NAs
      if (is.na(stdev)) stdev = 0
      stdevs[i - 2] = stdev
      
      #diffs: vector of vectors of differences; each index is a SUID; each index in a vector is a difference between adjacent replicates
      diffs <- tapply(as.numeric(normData[,i]),as.vector(normData[,1]), diff, na.rm=TRUE)
      #arrayDiffs: vector of all differences between duplicate spots for this array
      arrayDiffs <- (abs(na.omit(as.numeric(as.vector(unlist(diffs))))))
      allnspots[i - 2] <- length(arrayDiffs)
      medians[i - 2] <- median(arrayDiffs)
    }
    #plot stardard deviations vs. MAD
    plot(x = stdevs, y = medians, type = "p", main = "Array Quality", xlab = "Stdev of red/green ratio", 
         ylab = "Median Absolute Difference between duplicate spots", ylim = c(0,max(medians) + 0.05), xlim = c(0,max(stdevs) + 0.25), 
         col = colors[j])
    par(new = TRUE)
    infoFile[,j * 3 ] <- c(stdevs, mean(stdevs))
    infoFile[,j * 3 - 1] <- c(medians,mean(medians))
    infoFile[,j * 3 + 1] <- c(allnspots,mean(allnspots))
  }
  legend(x = "topright", legend = files, fill = colors, cex =.7)
  par(new = FALSE)
  write.table(infoFile, infoFileOutfilename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  dev.off()
}


filter <- function(filterFiles,galFile=NULL,fileType=".gpr",oligoIndex=NULL,filterControls=FALSE,filterFlaggedSpots=TRUE,filterPercentPresent=FALSE,
                   percentPresentThreshold=TRUE,filterRGR=FALSE,RGRthreshold=c(0.5,2),verbose=TRUE){
  # ARGUMENTS: 
  #      filterfiles (string vector) the files to be filtered.
  #      fileType (string) : either ".srr" or ".gpr"; the fileType ot use for filtering.
  #      oligoIndex (string) : the index in the gpr file of the column containing oligoIDs; only needed for control filtering.
  #      filterFlaggedSpots (logical) Should flagged spots be filtered out after normalization?
  #      filterPercentPresent(logical)ShouldspotswithfewerthanpercentPresentThresholdspotspresent across arrays be filtered out after normalization?
  #      percentPresentThreshold (number between 0 and 1) the lower bound of percent of spots present that is allowed in the filtered data.
  #      filterRGR (logical) should spots with red-green ratios between the parameters of RGRthreshold be filtered out?
  #      RGRthreshold (2-element numerical vector) Spots with red-green ratio between RGRthreshold[1] and RGRthreshold[2] are filtered out. 
  #      verbose (logical) : should progress reports be printed to the console?
  #error checking
  if (!file.exists(galFile)) stop ("filter: galfile",galFile,"not found")
  if (fileType != ".gpr" && fileType != ".srr") stop("filter fileType can be .srr or .gpr. You provided ", fileType)
  require(limma)
  outputPrefix <- makeFilterPrefix(filterPercentPresent,filterFlaggedSpots,filterRGR,filterControls,percentPresentThreshold,RGRthreshold)  
  
  #if I plan to remove controls, find out which indexes are controls by looking at oligoIDs.
  if (filterControls){
    if (is.null(oligoIndex)) stop("filter: specify which column the oligo ID is in the galfile")
    gal <- read.delim("HOANgal.txt")
    oligoIDs <- gal[56:nrow(gal),oligoIndex]
    controls <- (substr(oligoIDs,start=2,stop=2) == "C")
  }
  
  #assign which indexes in the file to look at based on whether they are .gpr or .srr files.
  if (fileType == ".gpr"){
    FlaggedIndex <- 50
    RGRIndex <- 35
  }
  if (fileType == ".srr"){ 
    FlaggedIndex <- 6
    RGRIndex <- 78
  }
  
  # find out which arrays we are analyzing
  filterfile <- read.delim(filterFiles[1])
  
  arraynames <- lapply(colnames(filterfile),function(colname){paste(c(colname,fileType),sep="",collapse="")})[3:ncol(filterfile)]
  #read GPR data into allGPRdata for use in filtering
  allGPRdata <- lapply(arraynames,function(filename,verbose){ 
    if (verbose) cat ("reading", filename,"\n")
    as.matrix(read.delim(filename,header=FALSE,skip=33))
  },verbose)
  
  #filter each file with a call to filterOneFile
  sapply(filterFiles,filterOneFile,normdata,allGPRdata,filterFlaggedSpots,filterRGR,filterControls,RGRthreshold,
         filterPercentPresent,percentPresentThreshold,verbose,outputPrefix,FlaggedIndex,RGRIndex,controls)
}


filterOneFile<-function(fileName,normdata,allGPRdata,filterFlaggedSpots, filterRGR,filterControls,RGRthreshold,
                        filterPercentPresent,percentPresentThreshold,verbose,outputPrefix,FlaggedIndex,RGRIndex,controls){
  #filter a single file. yay for decomposition.
  if (verbose) cat(paste(c("filtering ",fileName,"\n")))
  outputFilename <- paste(c(outputPrefix, fileName), sep = "", collapse = "")
  
  #normdata: data matrix with unfiltered data, read from user's input files
  normdata <- read.delim(fileName,header=T,comment.char="")  
  filteredNormData <- c()
  
  for(i in 1:length(allGPRdata)){ #for each array...
    #normvector : the column of normdata currently being filtered
    normvector <- normdata[,i + 2]
    GPRdata <- allGPRdata[[i]]
    
    #flaggedFailed, RGRfailed : spots that fail the flag and RGR criteria; will be filtered out later
    if (filterFlaggedSpots) flaggedFailed <- which(as.numeric(GPRdata[,FlaggedIndex])!=0)
    if (filterRGR) {
      RGRgreater <- RGRless <- rep(0,nrow(GPRdata))
      RGRgreater[which(as.numeric(GPRdata[,RGRIndex])>=RGRthreshold[1])] <- 1
      RGRless[which(as.numeric(GPRdata[,RGRIndex]) <=RGRthreshold[2])] <- 1
      RGRfailed <- RGRgreater & RGRless
    }
    #filter the column according to parameters
    if (filterFlaggedSpots)  normvector <- removeRows(flaggedFailed, normvector) 
    if (filterControls)   normvector <- removeRows(controls,normvector)     
    if (filterRGR) normvector <- removeRows(RGRfailed, normvector)
    suppressWarnings(filteredNormData <- cbind(filteredNormData,normvector))
  }
  
  #after filtering all other parameters, go back and filter out rows that fail the percent present threshold.
  if (filterPercentPresent) {
    #calculate which rows fail the percent present threshold
    percentPresentFailed <- percentIsPresent(percentPresentThreshold, filteredNormData)
    for (i in 1:ncol(filteredNormData)){
      if (filterPercentPresent) filteredNormData[,i] <- removeRows(percentPresentFailed, filteredNormData[,i])
    }
  }
  # insert gene names and IDs in front of the filtered data; assign column names and write to file
  filteredNormData <- cbind(as.character(normdata[,1]),as.character(normdata[,2]),filteredNormData)
  
  colnames(filteredNormData) <- colnames(normdata)
  write.table(filteredNormData,outputFilename,col.names=TRUE,row.names=FALSE,quote=F,sep="\t",append=F)
  return(outputFilename)
}

removeRows <- function (parameter, normvector){
  #FUNCTION: removeRows
  # filters out undesired spots from a single column by changing their values to NA.
  #ARGUMENTS:
  #    parameter: (logical vector the same length as normvector) - TRUE when spots should be filtered out,
  #             FALSE when spots should be retained.
  #    normvector: (numerical vector) - a single column (array) of unfiltered data.
  normvector[parameter] <- NA
  return(normvector)
}


percentIsPresent <- function(threshold, normdata){
  #FUNCTION: pctpresentrow
  # calculates whether a row of data should be retained based on whether or not it has at least
  #  threshold% of spots present.
  #ARGUMENTS:
  #    row: (numerical vector) a single row of data (corresponding to a spot across multiple arrays)
  #     threshold (numerical, between 0 and 1): what proportion of data must be present in order to retain it?
  cat ("calculating percent present per row \n")
  failedrows <- vector()
  applied <- apply(normdata, 1, pctpresentrow, threshold)
  failedrows <- which(applied == FALSE)
  return(failedrows)
}


makeFilterPrefix <- function(filterPercentPresent,filterFlaggedSpots,filterRGR,filterControls,percentPresentThreshold,RGRthreshold){
  #FUNCTION: makeFilterprefix
  # returns an appropriate prefix to add to filtering outFileNames to reflect how the data was filtered.
  #ARGUMENTS:
  #    filterPercentPresent,filterFlaggedSpots,filterRGR,filterControls (logical) : filtering treatments.
  outputPrefix <- "filtered-"
  if (filterPercentPresent) outputPrefix <- paste(c(outputPrefix,100*percentPresentThreshold,"p"),sep="",collapse="")
  if (filterFlaggedSpots) outputFilename <- paste(c("-Flag",  outputPrefix), sep = "", collapse = "")
  if (filterRGR) outputFilename <- paste(c("-RGR",RGRthreshold[1],"_",RGRthreshold[2],  outputPrefix), sep = "", collapse = "")         
  if (filterControls) outputFilename <- paste(c("Control",  outputPrefix), sep = "", collapse = "")               
  return(paste(c(outputPrefix,"_"),sep="",collapse=""))
}


subd0 <- function(files,ZeroCols=c(3,4)){
  # ARGUMENTS: 
  #      files (string vector) : the files to be processed.
  #       ZeroCols (numerical vector) : the columns in the input files that correspond to day 0 samples. 
  
  #error checking
  for (i in 1:length(files)) if (!file.exists(files[i])) stop (c("file does not exist : ", files[i]))
  #data : list of data.frames of input data
  data <- lapply(files,function(f){read.delim(f)})
  outFileNames <- vector()
  for (i in 1:length(data)){
    #d0s: matrix containing expression data for all day-0 samples
    d0s <- matrix(nrow = nrow(data[[i]]), ncol = length(ZeroCols))
    for (c in 1:length(ZeroCols))  d0s[,c] <- data[[i]][,ZeroCols[c]]
    
    #avgs: numerical vector of the average expression data for all day-0 samples
    avgs <- sapply(1:nrow(d0s),function(row,d0s){mean(na.omit(d0s[row,]))},d0s)
    
    #subtract the average of day-0 expression from all arrays.
    for (j in 3:ncol(data[[i]])) data[[i]][,j] <- data[[i]][,j] - avgs
    
    outFileNames <- c(outFileNames,paste(c("d0subbed_",files[i]),sep="",collapse=""))
    write.table(data[[i]],outFileNames[i],row.names=FALSE,quote=FALSE,sep="\t")
  }
  return(outFileNames)
}
