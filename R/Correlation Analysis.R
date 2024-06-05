library(pacman)
p_load("ggplot2", "dplyr", "clinfun", "reshape2", "MASS", "ggrepel", "seqinr", "metap", 'see', 'gghalves', 'ggcorrplot', 'extrafont')
font_import()
loadfonts(device = "win")

# Loads the original dataset (i.e. the pNpS one, and then converts the hit and length data into a mutsPerBase statistic)
pNpSdata <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/Combined NS-MA-TajD Raw Data v3.csv")
pNpSdata <- data.frame(gene = pNpSdata$geneName, mutsPerBase = pNpSdata$numberOfHits/pNpSdata$length, numberOfHits = pNpSdata$numberOfHits, length = pNpSdata$length, ratio = pNpSdata$ratio, dN = pNpSdata$dN, dS = pNpSdata$dS, tajD = pNpSdata$tajD, piNorm = pNpSdata$piNorm, mutationType = pNpSdata$mutationType, chromosome = pNpSdata$chromosome, ma = pNpSdata$MAGene)

# Loads the dNdS dataset, before removing all values that *don't* have a dNdS ratio (dRatio) 
dNdSdata <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/Combined Selection Measures v3.csv") # Note that the reported mutation rate is specific to C. elegans, while the dN/dS ratio is between C. elegans and C. briggsae
dNdSdata <- dNdSdata[which(!is.na(dNdSdata$dRatio)),]
dNdSdata$McKt <- dNdSdata$dRatio / dNdSdata$pRatio # Calculates the McDonald-Kreitman (McKt) test statistic 

mapNpSdata <- pNpSdata[which(pNpSdata$ma == 'TRUE'),] # Takes the pNpS data, and filters out all the genes that are only present in the CaeNDR dataset
mapNpSdata$numberOfHits <- mapNpSdata$numberOfHits - 1 # Remove the extra +1 given to all the data
mapNpSdata$mutsPerBase <- mapNpSdata$numberOfHits / mapNpSdata$length # Recalculate the mutation rate following this removal

madNdSdata <- dNdSdata[which(dNdSdata$ma == 'TRUE'),] # Takes the pNpS data, and filters out all the genes that are only present in the CaeNDR dataset
madNdSdata$numberOfHits <- madNdSdata$numberOfHits - 1 # Remove the extra +1 given to all the data
madNdSdata$mutsPerBase <- madNdSdata$numberOfHits / madNdSdata$length # Recalculate the mutation rate following this removal

#'* Performs correlation testing (both Spearman's and Pearson's) on the data provided (where x and y are the two sets to be correlated). If logged = TRUE, the Pearson's values are log10'd *
corTesting <- function(x, y, log = T) {
  corTest <- cor.test(x, y, method = "spearman")
  print(paste0("Spearman’s - rho = ", signif(corTest$estimate, 5), " (p-value of ", signif(corTest$p.value, 5), ")"))
  if (log == T) {
    corTest <- cor.test(log10(x), log10(y), method = "pearson")
    print(paste0("Pearson’s (log10) - corr = ", signif(corTest$estimate, 5), " (p-value of ", signif(corTest$p.value, 5), ")"))
  } else {
    corTest <- cor.test(x, y, method = "pearson")
    print(paste0("Pearson’s - corr = ", signif(corTest$estimate, 5), " (p-value of ", signif(corTest$p.value, 5), ")"))
  }
}

#'* Performs T- and Wilcoxon testing on the data provided (where x and y are two sets to be correlated). If pair = TRUE, the paired versions of the tests are used. The difference in means (for the T-test) is calculated as the mean of x minus the mean of y *
tWilTesting <- function(x, y, pair = F) {
  wilTest <- wilcox.test(x, y, paired = pair)
  print(paste0("Wilcoxon - V = ", signif(wilTest$statistic, 5), " (p-value of ", signif(wilTest$p.value, 5), ")"))
  tTest <- t.test(x, y, paired = pair)
  print(paste0("T-Test - T = ", signif(tTest$statistic, 5), " (p-value of ", signif(tTest$p.value, 5), ", Δμ = ", signif(mean(x) - mean(y),5), ")"))
}

#'* Performs a Jonckheere test on the data provided (where x is a binned dataframe, variable is the factor to test the trend of, and desiredBins are a list of bin names/numbers to be considered) *
jonckTesting <- function(x, variable, desiredBins) {
  jonckData <- x[which(x$binNum %in% desiredBins),]
  groups <- c(which(jonckData$binNum == desiredBins[1]))
  for (i in 2:length(desiredBins)) {
    groups <- append(groups, which(jonckData$binNum == desiredBins[i]))
  }
  jonckTest <- jonckheere.test(jonckData[,variable], g = groups)
  print(paste0("Jonckheere - JT = ", round(jonckTest$statistic, 5), " (p-value of ", signif(jonckTest$p.value, 5), ")"))
}

#'* Bins by the pN/pS ratio into specifically assigned intervals (as suggested by Peter), adding the binNum to a new column of the dataframe *
binBypNpS <- function(dataFrame) {
  dataFrame <- dataFrame[order(dataFrame[,"ratio"], decreasing=T),]
  for (i in 1:nrow(dataFrame)) { #'[ For each gene, figure out what bin it fits into, and assign it the correct selection type (binNum, a new column of the dataframe) ]
    if (dataFrame$ratio[i] < 0.05) { # Extremely purifying selection (pN/pS <= 0.05)
      dataFrame$binNum[i] <- "Extremely Purifying"
    } else if (0.05 <= dataFrame$ratio[i] & dataFrame$ratio[i] < 0.5) { # Strongly purifying selection (0.05 <= pN/pS < 0.5)
      dataFrame$binNum[i] <- "Strongly Purifying" 
    } else if (0.5 <= dataFrame$ratio[i] & dataFrame$ratio[i] < 0.8) { # Weakly purifying selection (0.5 <= pN/pS < 0.8)
      dataFrame$binNum[i] <- "Weakly Purifying"
    } else if (0.8 <= dataFrame$ratio[i] & dataFrame$ratio[i] < 1.2) { # ~Neutral (0.8 <= pN/pS < 1.2)
      dataFrame$binNum[i] <- "Neutral"
    } else if (1.2 <= dataFrame$ratio[i] & dataFrame$ratio[i] < 2.5) { # Diversifying (1.2 <= pN/pS < 2.5)
      dataFrame$binNum[i] <- "Diversifying"
    } else  { # Extremely diversifying selection (2.5 <= pN/pS)
      dataFrame$binNum[i] <- "Extremely Diversifying"
    }
  }
  dataFrame$binNum <- factor(dataFrame$binNum, levels = c("Extremely Purifying", "Strongly Purifying", "Weakly Purifying", "Neutral", "Diversifying", "Extremely Diversifying"), ordered = T)
  return(dataFrame)
}
desiredpNpSRatios = c("Strongly Purifying", "Weakly Purifying", "Neutral", "Diversifying")

#'* Bins by the dN/dS ratio into specifically assigned intervals, adding the binNum to a new column of the dataframe *
binBydNdS <- function(dataFrame) {
  dataFrame <- dataFrame[order(dataFrame[,"dRatio"], decreasing=T),]
  dataFrame <- dataFrame[complete.cases(dataFrame$dRatio),]
  for (i in 1:nrow(dataFrame)) { #'[ For each gene, figure out what bin it fits into, and assign it the correct selection type (binNum, a new column of the dataframe) ]
    if (dataFrame$dRatio[i] < 0.03) { # Extremely purifying selection (dN/dS <= 0.03)
      dataFrame$binNum[i] <- "Extremely Purifying"
    } else if (0.03 <= dataFrame$dRatio[i] & dataFrame$dRatio[i] < 0.08) { # Strongly purifying selection (0.03 <= dN/dS < 0.08)
      dataFrame$binNum[i] <- "Strongly Purifying" 
    } else if (0.08 <= dataFrame$dRatio[i] & dataFrame$dRatio[i] < 0.15) { # Purifying selection (0.08 <= dN/dS < 0.15)
      dataFrame$binNum[i] <- "Purifying"
    } else if (0.15 <= dataFrame$dRatio[i] & dataFrame$dRatio[i] < 0.25) { # Weakly purifying selection (0.15 <= dN/dS < 0.25)
      dataFrame$binNum[i] <- "Weakly Purifying"
    } else  { # Diversifying selection (0.25 <= dN/dS)
      dataFrame$binNum[i] <- "Diversifying"
    }
  }
  dataFrame$binNum <- factor(dataFrame$binNum, levels = c("Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying", "Diversifying"), ordered = T)
  return(dataFrame)
}
desireddNdSRatios = c( "Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying")

#'* Bins by the McDonald-Kreitman (McKt) test statistic into specifically assigned intervals, adding the binNum to a new column of the dataframe *
binByMcKt <- function(dataFrame) {
  dataFrame <- dataFrame[order(dataFrame[,"McKt"], decreasing=T),]
  for (i in 1:nrow(dataFrame)) { #'[ For each gene, figure out what bin it fits into, and assign it the correct selection type (binNum, a new column of the dataframe) ]
    if (dataFrame$McKt[i] < 0.02) { # Extremely purifying selection (McKt <= 0.05)
      dataFrame$binNum[i] <- "Extremely Purifying"
    } else if (0.02 <= dataFrame$McKt[i] & dataFrame$McKt[i] < 0.15) { # Strongly purifying selection (0.02 <= McKt < 0.15)
      dataFrame$binNum[i] <- "Strongly Purifying" 
    } else if (0.15 <= dataFrame$McKt[i] & dataFrame$McKt[i] < 0.4) { # Purifying selection (0.15 <= McKt < 0.4)
      dataFrame$binNum[i] <- "Purifying"
    } else { # Weakly purifying selection (0.4 <= McKt)
      dataFrame$binNum[i] <- "Weakly Purifying"
    }
  }
  dataFrame$binNum <- factor(dataFrame$binNum, levels = c("Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying"), ordered = T)
  return(dataFrame)
}

#'* Bins by the Tajima's D value into specifically assigned intervals *
binByTajD <- function(dataFrame, zCut = 0) {
  dataFrame <- dataFrame[complete.cases(dataFrame$tajD),] # Remove any genes that don't have D values (in the CaeNDR data, 1,099 don't)
  dataFrame <- dataFrame[order(dataFrame[,"tajD"], decreasing=T),] # Order the dataframe (for neatness' sake)
  
  if (zCut > 0) { # If you want to have a Z-value cut-off, this number is set to greater than 0 at the function call (2 is usually a good value)
    tajDMean = mean(dataFrame$tajD) # Calculate the mean of the D values
    tajDSD = sd(dataFrame$tajD) # Calculate the SD of the D values
    dataFrame <- dataFrame[abs((dataFrame$tajD - tajDMean) / tajDSD) <= zCut,] # Remove any genes that have D values that fall outside of the specified range
  }
  
  for (i in 1:nrow(dataFrame)) { #'[ For each gene, figure out what bin it fits into, and assign it the correct selection type (binNum, a new column of the dataframe) ]
    if (dataFrame$tajD[i] < -2.5) {
      dataFrame$binNum[i] <- "Extremely Purifying"
    } else if (-2.5 <= dataFrame$tajD[i] & dataFrame$tajD[i] < -1.8) {
      dataFrame$binNum[i] <- "Strongly Purifying"
    } else if (-1.8 <= dataFrame$tajD[i] & dataFrame$tajD[i] < -1.0) {
      dataFrame$binNum[i] <- "Purifying"
    } else if (-1.0 <= dataFrame$tajD[i] & dataFrame$tajD[i] < -0.3) {
      dataFrame$binNum[i] <- "Weakly Purifying"
    } else if (-0.3 <= dataFrame$tajD[i] & dataFrame$tajD[i] < 0.3) {
      dataFrame$binNum[i] <- "Neutral"
    } else if (0.3 <= dataFrame$tajD[i]) {
      dataFrame$binNum[i] <- "Balancing"
    }
  }
  dataFrame$binNum <- factor(dataFrame$binNum, levels = c("Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying", "Neutral", "Balancing"))
  return(dataFrame) # Return the dataframe, now with bin values
}
desiredTajDs = c("Strongly Purifying", "Purifying", "Weakly Purifying", "Neutral", "Balancing")

#'* Bins by the pi value into specifically assigned intervals. Note that it also LOGS THE PI VALUES *
binByPi <- function(dataFrame, zCut = 0) {
  dataFrame <- dataFrame[which(complete.cases(dataFrame$piNorm) & dataFrame$piNorm != 0),] # Remove any genes that don't have pi values, or where the value is 0 (as this won't be loggable) 
  dataFrame$piLog <- log10(dataFrame$piNorm)
  dataFrame <- dataFrame[order(dataFrame[,"piNorm"], decreasing=T),]
  # dataFrame$piNorm <- log10(dataFrame$piNorm)
  # piMeanLog <- mean(dataFrame$piNorm, na.rm = TRUE)
  # piSDLog <- sd(dataFrame$piNorm, na.rm = TRUE)
  # dataFrame <- dataFrame[abs((dataFrame$piNorm - piMeanLog) / piSDLog) <= 2,]
  
  if (zCut > 0) { # If you want to have a Z-value cut-off, this number is set to greater than 0 at the function call (2 is usually a good value)
    piLogMean = mean(dataFrame$piLog) # Calculate the mean of the pi values
    piLogSD = sd(dataFrame$piLog) # Calculate the SD of the pi values
    dataFrame <- dataFrame[abs((dataFrame$piLog - piLogMean) / piLogSD) <= zCut,] # Remove any genes that have pi values that fall outside of the specified range
  }
  
  for (i in 1:nrow(dataFrame)) { #'[ For each gene, figure out what bin it fits into, and assign it the correct selection type (binNum, a new column of the dataframe) ]
    if (dataFrame$piLog[i] < -3.75) {
      dataFrame$binNum[i] <- "Extremely Purifying"
    } else if (-3.75 <= dataFrame$piLog[i] & dataFrame$piLog[i] < -3.25) {
      dataFrame$binNum[i] <- "Strongly Purifying"
    } else if (-3.25 <= dataFrame$piLog[i] & dataFrame$piLog[i] < -2.75) {
      dataFrame$binNum[i] <- "Purifying"
    } else if (-2.75 <= dataFrame$piLog[i]) {
      dataFrame$binNum[i] <- "Weakly Purifying"
    }
  }
  dataFrame$binNum <- factor(dataFrame$binNum, levels = c("Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying"))
  return(dataFrame)
}

#'* Plots a particular (specified) measure of (purifying) selection * 
plotBoxViolinPlot <- function(data, selectionMeasure, datasetName, xLabels) {
  if (length(xLabels) == 4) { # Gets the appropriate number of colours based off the number of x-axis labels provided
    colours <- c('#c92e62', '#54850f', '#7228a0', '#216392')
  } else if (length(xLabels) == 5) {
    colours <- c('#c92e62', '#c9832f', '#54850f', '#7228a0', '#216392')
  } else if (length(xLabels) == 6) {
    colours <- c('#c92e62', '#c9832f', '#54850f', '#288879', '#7228a0', '#216392')
  }
  
  # Plots the data as a half-boxplot, half-violinplot
  ggplot(data = data, aes(x = binNum, y = mutsPerBase, color = as.factor(binNum)))  + 
    geom_half_boxplot(outlier.shape = NA, position = position_nudge(x = 0, y = 0), width = 0.5, lwd = 1, fatten = 0.8, fill = '#AFAFAF') + 
    geom_violinhalf(position = position_nudge(x = 0, y = 0), lwd = 0.8, fill = '#AFAFAF') +
    labs(title = paste0("Mutation Rate of ", datasetName, " Genes for the Different Types of Selection (Binned by ", selectionMeasure, ")"), x = paste0(selectionMeasure, " Bin (Range)"), y = "Mutation Rate (per Gene)", "color" = "Selection Type", caption = paste0("n = ", nrow(data))) + 
    scale_x_discrete(labels = xLabels) +
    scale_colour_manual(values = colours) +
    theme(text=element_text(size=14,  family="Arial Nova"),
      axis.text = element_text(colour = '#000000'),
      panel.grid.major = element_line(linewidth = 1, colour = '#CFCFCF'), panel.grid.minor = element_line(linewidth = 0.2, colour = '#CFCFCF'),  
      panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
      legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
      plot.background = element_rect(fill = '#CFCFCF')) +
    scale_y_continuous(limits = c(0, 0.002), minor_breaks = seq(0, 0.003, 0.0001))
}

outputData <- function(caendrData, maData, selectionMeasure, xLabels) {
  if (selectionMeasure == 'pN/pS') {
    corTesting(caendrData$ratio, caendrData$mutsPerBase, F) # Calculate the correlations between the selection measure and the gene-specific mutation rate for the CaeNDR data
    caendrData <- binBypNpS(caendrData) # Bin the CaeNDR data by the selection measure
    jonckTesting(caendrData, 'mutsPerBase', desiredpNpSRatios) # Calculate the Jonckheere statistic for the CaeNDR data

    print('---')
    
    corTesting(maData$ratio, maData$mutsPerBase, F) # Calculate the correlations between the selection measure and the gene-specific mutation rate for the MA data
    maData <- binBypNpS(maData) # Bin the MA data by the selection measure
    jonckTesting(maData, 'mutsPerBase', desiredpNpSRatios) # Calculate the Jonckheere statistic for the MA data
    
    caendrData <- caendrData[caendrData$binNum %in% desiredpNpSRatios,] # Abbreviate the CaeNDR data to only include the bins of interest (for the purposes of plotting)
    maData <- maData[maData$binNum %in% desiredpNpSRatios,] # Abbreviate the MA data to only include the bins of interest (for the purposes of plotting)
  } else if (selectionMeasure == "Tajima's D") {
    corTesting(caendrData$tajD, caendrData$mutsPerBase, F) # Calculate the correlations between the selection measure and the gene-specific mutation rate for the CaeNDR data
    caendrData <- binByTajD(caendrData) # Bin the CaeNDR data by the selection measure
    jonckTesting(caendrData, 'mutsPerBase', desiredTajDs) # Calculate the Jonckheere statistic for the CaeNDR data
    
    print('---')
    
    corTesting(maData$tajD, maData$mutsPerBase, F) # Calculate the correlations between the selection measure and the gene-specific mutation rate for the MA data
    maData <- binByTajD(maData) # Bin the MA data by the selection measure
    jonckTesting(maData, 'mutsPerBase', desiredTajDs) # Calculate the Jonckheere statistic for the MA data
    
    caendrData <- caendrData[caendrData$binNum %in% desiredTajDs,] # Abbreviate the CaeNDR data to only include the bins of interest (for the purposes of plotting)
    maData <- maData[maData$binNum %in% desiredTajDs,] # Abbreviate the MA data to only include the bins of interest (for the purposes of plotting)
  } else if (selectionMeasure == "π") {
    corTesting(caendrData$piNorm, caendrData$mutsPerBase, F) # Calculate the correlations between the selection measure and the gene-specific mutation rate for the CaeNDR data
    caendrData <- binByPi(caendrData) # Bin the CaeNDR data by the selection measure
    jonckTesting(caendrData, 'mutsPerBase', c("Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying")) # Calculate the Jonckheere statistic for the CaeNDR data
    
    print('---')
    
    corTesting(maData$piNorm, maData$mutsPerBase, F) # Calculate the correlations between the selection measure and the gene-specific mutation rate for the MA data
    maData <- binByPi(maData) # Bin the MA data by the selection measure
    jonckTesting(maData, 'mutsPerBase', c("Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying")) # Calculate the Jonckheere statistic for the MA data
  } else if (selectionMeasure == "dN/dS") {
    corTesting(caendrData$dRatio, caendrData$mutsPerBase, F) # Calculate the correlations between the selection measure and the gene-specific mutation rate for the CaeNDR data
    caendrData <- binBydNdS(caendrData) # Bin the CaeNDR data by the selection measure
    jonckTesting(caendrData, 'mutsPerBase', desireddNdSRatios) # Calculate the Jonckheere statistic for the CaeNDR data
    
    print('---')
    
    corTesting(maData$dRatio, maData$mutsPerBase, F) # Calculate the correlations between the selection measure and the gene-specific mutation rate for the MA data
    maData <- binBydNdS(maData) # Bin the MA data by the selection measure
    jonckTesting(maData, 'mutsPerBase', desireddNdSRatios) # Calculate the Jonckheere statistic for the MA data
    
    caendrData <- caendrData[caendrData$binNum %in% desireddNdSRatios,] # Abbreviate the CaeNDR data to only include the bins of interest (for the purposes of plotting)
    maData <- maData[maData$binNum %in% desireddNdSRatios,] # Abbreviate the MA data to only include the bins of interest (for the purposes of plotting)
  } else if (selectionMeasure == "McDonald-Kreitman [McKt]") {
    corTesting(caendrData$McKt, caendrData$mutsPerBase, F) # Calculate the correlations between the selection measure and the gene-specific mutation rate for the CaeNDR data
    caendrData <- binByMcKt(caendrData) # Bin the CaeNDR data by the selection measure
    jonckTesting(caendrData, 'mutsPerBase', c("Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying")) # Calculate the Jonckheere statistic for the CaeNDR data
    
    print('---')
    
    corTesting(maData$McKt, maData$mutsPerBase, F) # Calculate the correlations between the selection measure and the gene-specific mutation rate for the MA data
    maData <- binByMcKt(maData) # Bin the MA data by the selection measure
    jonckTesting(maData, 'mutsPerBase', c("Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying")) # Calculate the Jonckheere statistic for the MA data
  }
  
  print(plotBoxViolinPlot(caendrData, selectionMeasure, "AG", xLabels))
  print(plotBoxViolinPlot(maData, selectionMeasure, "MA", xLabels))
}

outputData(pNpSdata, mapNpSdata, "pN/pS", c("0.05 <= pN/pS <0.5", "0.5 <= pN/pS < 0.8", "0.8 <= pN/pS < 1.2", "1.2 <= pN/pS < 2.5"))
outputData(pNpSdata, mapNpSdata, "Tajima's D", c("-2.5 <= tajD < -1.8", "-1.8 <= tajD < -1.0", "-1.0 <= tajD < -0.3", "-0.3 <= tajD < 0.3", "0.3 <= tajD"))
outputData(pNpSdata, mapNpSdata, "π", c("log10(π) < -3.75", "-3.75 <= log10(π) < -3.25", "-3.25 <= log10(π) < -2.75", "-2.75 <= log10(π)"))
outputData(dNdSdata, madNdSdata, "dN/dS", c("dN/dS < 0.03", "0.03 <= dN/dS < 0.08", "0.08 <= dN/dS < 0.15", "0.15 <= dN/dS < 0.25"))
outputData(dNdSdata, madNdSdata, "McDonald-Kreitman [McKt]", c("McKt < 0.02", "0.02 <= McKt < 0.15", "0.15 <= McKt < 0.4", "0.4 <= McKt"))

count <- 0
for (i in 1:nrow(mapNpSdata)) {
  for (j in 1:nrow(mapNpSdata)) {
    if (mapNpSdata[i, 'gene'] == mapNpSdata[j, 'gene']) {
      count <- count + 1
    }
  }
}

plotCorrelationPlot <- function(method = 'spearman') {
  correlations <- c() 
  pValues <- c()
  selectionTypes <- c('ratio', 'tajD', 'piNorm', 'dRatio', 'McKt')
  for (i in 1:3) {
    corTest <- cor.test(mapNpSdata[,selectionTypes[i]], mapNpSdata$mutsPerBase, method = method)
    correlations <- append(correlations, corTest$estimate[[1]])
    pValues <- append(pValues, corTest$p.value)
    corTest <- cor.test(pNpSdata[,selectionTypes[i]], pNpSdata$mutsPerBase, method = method)
    correlations <- append(correlations, corTest$estimate[[1]])
    pValues <- append(pValues, corTest$p.value)
  }
  for (i in 1:2) {
    corTest <- cor.test(madNdSdata[,selectionTypes[i + 3]], madNdSdata$mutsPerBase, method = method)
    correlations <- append(correlations, corTest$estimate[[1]])
    pValues <- append(pValues, corTest$p.value)
    corTest <- cor.test(dNdSdata[,selectionTypes[i + 3]], dNdSdata$mutsPerBase, method = method)
    correlations <- append(correlations, corTest$estimate[[1]])
    pValues <- append(pValues, corTest$p.value)
  }
  
  correlations <- matrix(correlations, nrow = 2, ncol = 5)
  rownames(correlations) <- c('MA', 'AG')
  colnames(correlations) <- c('pN/pS', 'tajD', 'π', 'dN/dS', 'McKt')
  correlations <- t(correlations)
  pValues <- matrix(pValues, nrow = 2, ncol = 5)
  rownames(pValues) <- c('MA', 'AG')
  colnames(pValues) <- c('pN/pS', 'tajD', 'π', 'dN/dS', 'McKt')
  pValues <- t(pValues)
  
  correlationLabels <- correlations
  for (i in 1:5) {
    if (pValues[i, 'AG'] <= 0.05) {
      correlationLabels[i, 'AG'] <- formatC(round(correlations[i,'AG'], 3), format = 'f', flag='0', digits = 3)
    } else {
      correlationLabels[i, 'AG'] <- NA
    }
    if (pValues[i, 'MA'] <= 0.05) {
      correlationLabels[i, 'MA'] <- formatC(round(correlations[i,'MA'], 3), format = 'f', flag='0', digits = 3)
    } else {
      correlationLabels[i, 'MA'] <- NA
    }
    
  }
  
  correlationPlot <- ggcorrplot(correlations, p.mat = pValues, insig = 'blank', digits = 3, outline.color = '#4F4F4F') + 
    scale_fill_gradient(high = '#c92e62', low = '#216392', limit=c(-0.025, 0.215)) + 
    theme(text=element_text(size=14,  family="Arial Nova"),
          axis.text = element_text(colour = '#000000'),
          panel.grid.major = element_line(colour = '#CFCFCF'), panel.grid.minor = element_line(colour = '#CFCFCF'),
          panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
          legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
          plot.background = element_rect(fill = '#CFCFCF')) +
    labs(x = 'Measure of Selection', y = 'Dataset')
  if (method == 'spearman') {
    correlationPlot <- correlationPlot + labs(title = "Spearman's Correlations Between Measures of Selection and Mutation Rate (per Gene)", fill = 'Correlation (rho)')
  } else if (method == 'pearson') {
    correlationPlot <- correlationPlot + labs(title = "Pearson's Correlations Between Measures of Selection and Mutation Rate (per Gene)", fill = 'Correlation (corr)')
  }
  for (i in 1:5) {
    correlationPlot <- correlationPlot + annotate('text', x = i, y = 2, label = correlationLabels[i,'AG'], col = '#CFCFCF', size = 6, family = "Arial Nova")
    correlationPlot <- correlationPlot + annotate('text', x = i, y = 1, label = correlationLabels[i,'MA'], col = '#CFCFCF', size = 6, family = "Arial Nova")
  }
  
  print(correlationPlot)
}

plotCorrelationPlot('spearman')
