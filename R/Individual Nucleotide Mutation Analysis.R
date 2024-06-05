# Takes an input .csv of genes that have mutation rates calculated for each of the four nucleotides, along with a dN/dS ratio, and performs statistical correlation testing on each rate

library(pacman)
p_load("ggplot2", "ggrepel", "dplyr", "reshape2", "clinfun", 'gghalves', 'see', 'ggcorrplot', 'extrafont')

data <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/Individual Nucleotide Mutation Rates (Genome, +1).csv")
# oldMaData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/OLD Individual MA Nucleotide Mutation Rates.csv")
# oldCaendrData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/OLD Individual CaeNDR Nucleotide Mutation Rates (+1).csv")
maData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/NEW Individual MA Nucleotide Mutation Rates.csv")
caendrData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/NEW Individual CaeNDR Nucleotide Mutation Rates (+1).csv")

#'* Bins by the dN/dS ratio into specifically assigned intervals (as suggested by Peter), adding the binNum to a new column of the dataframe *
binByOlddNdS <- function(dataFrame) {
  dataFrame <- dataFrame[order(dataFrame[,"ratio"], decreasing=T),]
  for (i in 1:nrow(dataFrame)) { #'[ For each gene, figure out what bin it fits into, and assign it the correct selection type (binNum, a new column of the dataframe) ]
    if (dataFrame$ratio[i] < 0.05) { # Extremely purifying selection (dN/dS <= 0.05)
      dataFrame$binNum[i] <- "Extremely Purifying"
    } else if (0.05 <= dataFrame$ratio[i] & dataFrame$ratio[i] < 0.5) { # Strongly purifying selection (0.05 <= dN/dS < 0.5)
      dataFrame$binNum[i] <- "Strongly Purifying" 
    } else if (0.5 <= dataFrame$ratio[i] & dataFrame$ratio[i] < 0.8) { # Weakly purifying selection (0.5 <= dN/dS < 0.8)
      dataFrame$binNum[i] <- "Weakly Purifying"
    } else if (0.8 <= dataFrame$ratio[i] & dataFrame$ratio[i] < 1.2) { # ~Neutral (0.8 <= dN/dS < 1.2)
      dataFrame$binNum[i] <- "Neutral"
    } else if (1.2 <= dataFrame$ratio[i] & dataFrame$ratio[i] < 2.5) { # Diversifying (1.2 <= dN/dS < 2.5)
      dataFrame$binNum[i] <- "Diversifying"
    } else  { # Extremely diversifying selection (2.5 <= dN/dS)
      dataFrame$binNum[i] <- "Extremely Diversifying"
    }
  }
  dataFrame$binNum <- factor(dataFrame$binNum, levels = c("Extremely Purifying", "Strongly Purifying", "Weakly Purifying", "Neutral", "Diversifying", "Extremely Diversifying"), ordered = T)
  return(dataFrame)
}

#'* Bins by the pN/pS ratio into specifically assigned intervals (as suggested by Peter), adding the binNum to a new column of the dataframe *
binBypNpS <- function(dataFrame) {
  dataFrame <- dataFrame[order(dataFrame[,"pRatio"], decreasing=T),]
  for (i in 1:nrow(dataFrame)) { #'[ For each gene, figure out what bin it fits into, and assign it the correct selection type (binNum, a new column of the dataframe) ]
    if (dataFrame$pRatio[i] < 0.05) { # Extremely purifying selection (pN/pS <= 0.05)
      dataFrame$binNum[i] <- "Extremely Purifying"
    } else if (0.05 <= dataFrame$pRatio[i] & dataFrame$pRatio[i] < 0.5) { # Strongly purifying selection (0.05 <= pN/pS < 0.5)
      dataFrame$binNum[i] <- "Strongly Purifying" 
    } else if (0.5 <= dataFrame$pRatio[i] & dataFrame$pRatio[i] < 0.8) { # Weakly purifying selection (0.5 <= pN/pS < 0.8)
      dataFrame$binNum[i] <- "Weakly Purifying"
    } else if (0.8 <= dataFrame$pRatio[i] & dataFrame$pRatio[i] < 1.2) { # ~Neutral (0.8 <= pN/pS < 1.2)
      dataFrame$binNum[i] <- "Neutral"
    } else if (1.2 <= dataFrame$pRatio[i] & dataFrame$pRatio[i] < 2.5) { # Diversifying (1.2 <= pN/pS < 2.5)
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

#'* Bins by the pi value into specifically assigned intervals. Note that it also LOGS THE PI VALUES *
binByPi <- function(dataFrame, zCut = 0) {
  dataFrame <- dataFrame[which(complete.cases(dataFrame$pi) & dataFrame$pi != 0),] # Remove any genes that don't have pi values, or where the value is 0 (as this won't be loggable) 
  dataFrame$piLog <- log10(dataFrame$pi)
  dataFrame <- dataFrame[order(dataFrame[,"pi"], decreasing=T),]
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


calculateNucleotideCorrelations <- function(data, selectionMeasure) {
  results <- data.frame('nucleotide' = c('A', 'C', 'G', 'T'), 'rho' = c(NA, NA, NA, NA), 'spearmanP-value' = c(NA, NA, NA, NA), 'corr' = c(NA, NA, NA, NA), 'pearsonP-value' = c(NA, NA, NA, NA))
  index <- 1
  
  if (selectionMeasure == 'pN/pS') {
    selectionMeasure = 'pRatio'
  } else if (selectionMeasure == 'π') {
    selectionMeasure = 'pi'
  } else if (selectionMeasure == 'dN/dS') {
    selectionMeasure = 'dRatio'
  }
  
  for (nuc in c('A', 'C', 'G', 'T')) {
    oneNucData <- data[which(!is.na(data[nuc])),c('gene', nuc, selectionMeasure)]
    oneNucData <- oneNucData[which(!is.na(oneNucData[selectionMeasure])),]
    names(oneNucData) <- c('gene', 'mutsPerBase', selectionMeasure)
    corTest <- cor.test(oneNucData[,selectionMeasure], oneNucData$mutsPerBase, method = "spearman")
    results$rho[index] <- signif(corTest$estimate, 5)
    results$spearmanP.value[index] <- signif(corTest$p.value, 5)
    corTest <- cor.test(oneNucData[,selectionMeasure], oneNucData$mutsPerBase, method = "pearson")
    results$corr[index] <- signif(corTest$estimate, 5)
    results$pearsonP.value[index] <- signif(corTest$p.value, 5)
    index <- index + 1
  }
  print(results)
}

plotNucleotideJoncks <- function(data, dataset, selectionMeasure, ymax) {
  if (selectionMeasure == 'pN/pS') {
    binnedData <- binBypNpS(data[complete.cases(data$pRatio),]) # Bin the data by pN/pS ratio, removing any values that don't have a pN/pS ratio (likely owing to WBID mismatches)
    bins <- c("Strongly Purifying", "Weakly Purifying", "Neutral", "Diversifying")
    xLabels <- c("0.05 <= pN/pS < 0.5", "0.5 <= pN/pS < 0.8", "0.8 <= pN/pS < 1.2", "1.2 <= pN/pS < 2.5")
  } else if (selectionMeasure == 'π') {
    binnedData <- binByPi(data) # Bin the data by pi (removal is done by the function)
    bins <- c("Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying")
    xLabels <- c("log10(π) < -3.75", "-3.75 <= log10(π) < -3.25", "-3.25 <= log10(π) < -2.75", "-2.75 <= log10(π)")
  } else if (selectionMeasure == 'dN/dS') {
    binnedData <- binBydNdS(data[complete.cases(data$dRatio),]) # Bin the data by dN/dS, removing any values that don't have one
    bins <- c("Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying")
    xLabels <- c("dN/dS < 0.03", "0.03 <= dN/dS < 0.08", "0.08 <= dN/dS < 0.15", "0.15 <= dN/dS < 0.25")
  }
  
  if (ymax <= 0.025) {
    minorLines <- 0.001
  } else if (0.025 < ymax & ymax <= 0.1) {
    minorLines <- 0.0025
  } else {
    minorLines <- 0.01
  }
  
  for (nuc in c('A', 'C', 'G', 'T')) { # For each nucleotide, calculate the Jonckheere value, print it, and plot the bins
    jonckData <- binnedData[which(!is.na(binnedData[,nuc]) & binnedData$binNum %in% bins),]
    jonckTesting(jonckData, nuc, bins)
    
    # Plot the data as a half-boxplot, half-violinplot
    boxPlot <- ggplot(data = jonckData, aes(x = binNum, y = .data[[nuc]], color = as.factor(binNum)))  + 
      geom_half_boxplot(outlier.shape = NA, position = position_nudge(x = 0, y = 0), width = 0.5, lwd = 1, fatten = 0.8, fill = '#AFAFAF') + 
      geom_violinhalf(position = position_nudge(x = 0, y = 0), lwd = 0.8, fill = '#AFAFAF') +
      labs(title = paste0("Individual ", nuc, " Mutation Rate of ", dataset, " Genes (Binned by ", selectionMeasure, ")"), 
           x = paste0(selectionMeasure, " Bin (Range)"), 
           y = "Mutation Rate per Nucleotide (per Gene)", "color" = "Selection Type", caption = paste0("n = ", nrow(jonckData))) +
      scale_x_discrete(labels = xLabels) +
      scale_colour_manual(values = c('#c92e62', '#54850f', '#7228a0', '#216392')) +
      theme(text=element_text(size=14,  family="Arial Nova"),
            axis.text = element_text(colour = '#000000'),
            panel.grid.major = element_line(linewidth = 1, colour = '#CFCFCF'), panel.grid.minor = element_line(linewidth = 0.2, colour = '#CFCFCF'),  
            panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
            legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
            plot.background = element_rect(fill = '#CFCFCF')) +
      scale_y_continuous(limits = c(0, ymax), minor_breaks = seq(0, 1, minorLines))
  
    print(boxPlot)
  }
}

plotCorrelationPlot <- function(data, dataset, method = 'spearman') {
  correlations <- c() 
  pValues <- c()
  for (measure in c('pRatio', 'pi', 'dRatio')) {
    for (nuc in c('A', 'C', 'G', 'T')) {
      corTest <- cor.test(data[,measure], data[,nuc], method = method)
      correlations <- append(correlations, corTest$estimate[[1]])
      pValues <- append(pValues, corTest$p.value)
    }
  }
  
  correlations <- matrix(correlations, nrow = 4, ncol = 3)
  rownames(correlations) <- c('A', 'C', 'G', 'T')
  colnames(correlations) <- c('pN/pS', 'π', 'dN/dS')
  pValues <- matrix(pValues, nrow = 4, ncol = 3)
  rownames(pValues) <- c('A', 'C', 'G', 'T')
  colnames(pValues) <- c('pN/pS', 'π', 'dN/dS')
  
  correlationLabels <- correlations
  for (i in 1:3) {
    for (nuc in c('A', 'C', 'G', 'T')) {
      if (pValues[nuc, i] <= 0.05) {
        correlationLabels[nuc, i] <- formatC(round(correlations[nuc, i], 3), format = 'f', flag='0', digits = 3)
      } else {
        correlationLabels[nuc, i] <- NA
      }
    }
  }
  
  correlationPlot <- ggcorrplot(correlations, p.mat = pValues, insig = 'blank', digits = 3, outline.color = '#4F4F4F', tl.srt = 0) + 
    # scale_fill_gradient(high = '#c92e62', low = '#216392', limits = c(-0.07, 0.57)) + 
    theme(text=element_text(size=14,  family="Arial Nova"),
          axis.text = element_text(colour = '#000000'), plot.title = element_text(size=13),
          panel.grid.major = element_line(colour = '#CFCFCF'), panel.grid.minor = element_line(colour = '#CFCFCF'),
          panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
          legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
          plot.background = element_rect(fill = '#CFCFCF')) +
    labs(x = 'Measure of Selection', y = 'Nucleotide')
  
  if (method == 'spearman') { # Set the title/legend heading based on the correlation method used
    correlationPlot <- correlationPlot + labs(title = paste0("Spearman's Correlations Between Selection Measures and Nucleotide-specific Mutation Rate (per Gene, ", dataset, " dataset)"), fill = 'Correlation (rho)')
  } else if (method == 'pearson') {
    correlationPlot <- correlationPlot + labs(title = paste0("Pearson's Correlations Between Selection Measures and Nucleotide-specific Mutation Rate (per Gene, ", dataset, " dataset)"), fill = 'Correlation (corr)')
  }
  
  if (dataset == 'AG') { # Set the range of the colours based on the dataset used (as they have drastically different ranges). This allows for comparison of correlation methods within the two datasets without MA looking like a homogeneous mess of colour
    correlationPlot <- correlationPlot + scale_fill_gradient(high = '#c92e62', low = '#216392', limits = c(0.08, 0.23))
  } else if (dataset == 'MA') {
    correlationPlot <- correlationPlot + scale_fill_gradient(high = '#c92e62', low = '#216392', limits = c(-0.06, 0.11)) # , limits = c(-0.07, 0.56)
  }
  
  for (i in 1:3) {
    for (j in 1:4) {
      correlationPlot <- correlationPlot + annotate('text', x = j, y = i, label = correlationLabels[j, i], col = '#CFCFCF', size = 6, family = "Arial Nova")
    }
  }
  
  print(correlationPlot)
}

plotCorrelationPlot(caendrData, 'AG', 'spearman')
plotCorrelationPlot(caendrData, 'AG', 'pearson')
plotCorrelationPlot(maData, 'MA', 'spearman')
plotCorrelationPlot(maData, 'MA', 'pearson')

calculateNucleotideCorrelations(caendrData, 'pN/pS')
calculateNucleotideCorrelations(maData, 'pN/pS')
calculateNucleotideCorrelations(caendrData, 'π')
calculateNucleotideCorrelations(maData, 'π')
calculateNucleotideCorrelations(caendrData, 'dN/dS')
calculateNucleotideCorrelations(maData, 'dN/dS')

plotNucleotideJoncks(caendrData, 'AG', 'pN/pS', 0.01)
plotNucleotideJoncks(caendrData, 'AG', 'π', 0.01)
plotNucleotideJoncks(caendrData, 'AG', 'dN/dS', 0.01)
plotNucleotideJoncks(maData, 'MA', 'pN/pS', 0.002)
plotNucleotideJoncks(maData, 'MA', 'π', 0.0002)
plotNucleotideJoncks(maData, 'MA', 'dN/dS', 0.0002)

binnedData <- binBydNdS(caendrData[!is.na(caendrData$pRatio),])
#- binnedData[(is.na(binnedData))] <- 0

# Calculates the Wilcoxon V value for the four different nucleotide-specific mutations 
wilcoxTesting <- function(data, lowBin, highBin) {
  results <- matrix(NA, nrow = 4, ncol = 4)
  rownames(results) <- c('A', 'C', 'G', 'T')
  colnames(results) <- c('wilcoxonV', 'pValue', 'lowMedian', 'highMedian')
  
  for (nuc in c('A', 'C', 'G', 'T')) {
    results[nuc, 'lowMedian'] <- round(median(data[which(!is.na(data[, nuc]) & data$binNum == lowBin), nuc]), 5)
    results[nuc, 'highMedian'] <- round(median(data[which(!is.na(data[, nuc]) & data$binNum == highBin), nuc]), 5)
    wilTest <- wilcox.test(data[which(!is.na(data[, nuc]) & data$binNum == lowBin), nuc], data[which(!is.na(data[, nuc]) & data$binNum == highBin), nuc], paired = F)
    results[nuc, 'wilcoxonV'] <- round(wilTest$statistic, 5)
    results[nuc, 'pValue'] <- signif(wilTest$p.value, 5)
  }

  print(results)
}

# Calculates the T-test T value for the four different nucleotide-specific mutations
tTesting <- function(data, lowBin, highBin) {
  results <- matrix(NA, nrow = 4, ncol = 4)
  rownames(results) <- c('A', 'C', 'G', 'T')
  colnames(results) <- c('tTestStat', 'pValue', 'lowMean', 'highMean')
  
  for (nuc in c('A', 'C', 'G', 'T')) {
    results[nuc, 'lowMean'] <- round(mean(data[which(!is.na(data[, nuc]) & data$binNum == lowBin), nuc]), 5)
    results[nuc, 'highMean'] <- round(mean(data[which(!is.na(data[, nuc]) & data$binNum == highBin), nuc]), 5)
    # print(round(sd(data[which(!is.na(data[, nuc]) & data$binNum == lowBin), nuc]), 5))
    # print(round(sd(data[which(!is.na(data[, nuc]) & data$binNum == highBin), nuc]), 5))
    
    tTest <- t.test(data[which(!is.na(data[, nuc]) & data$binNum == lowBin), nuc], data[which(!is.na(data[, nuc]) & data$binNum == highBin), nuc], paired = F, var.equal = F)
    results[nuc, 'tTestStat'] <- round(tTest$statistic, 5)
    results[nuc, 'pValue'] <- signif(tTest$p.value, 5)
  }
  
  print(results)
}

wilcoxTesting(binBypNpS(caendrData[!is.na(caendrData$pRatio),]), 'Strongly Purifying', 'Diversifying')
wilcoxTesting(binBypNpS(maData[!is.na(maData$pRatio),]), 'Strongly Purifying', 'Diversifying')
wilcoxTesting(binByPi(caendrData[!is.na(caendrData$pi),]), 'Extremely Purifying', 'Weakly Purifying')
wilcoxTesting(binByPi(maData[!is.na(maData$pi),]), 'Extremely Purifying', 'Weakly Purifying')
wilcoxTesting(binBydNdS(caendrData[!is.na(caendrData$dRatio),]), 'Extremely Purifying', 'Weakly Purifying')
wilcoxTesting(binBydNdS(maData[!is.na(maData$dRatio),]), 'Extremely Purifying', 'Weakly Purifying')

tTesting(binBypNpS(caendrData[!is.na(caendrData$pRatio),]), 'Strongly Purifying', 'Diversifying')
tTesting(binBypNpS(maData[!is.na(maData$pRatio),]), 'Strongly Purifying', 'Diversifying')
tTesting(binByPi(caendrData[!is.na(caendrData$pi),]), 'Extremely Purifying', 'Weakly Purifying')
tTesting(binByPi(maData[!is.na(maData$pi),]), 'Extremely Purifying', 'Weakly Purifying')
tTesting(binBydNdS(caendrData[!is.na(caendrData$dRatio),]), 'Extremely Purifying', 'Weakly Purifying')
tTesting(binBydNdS(maData[!is.na(maData$dRatio),]), 'Extremely Purifying', 'Weakly Purifying')












# loessLine <- predict(loess(log10(A) ~ log10(ratio), data = data, span = 0.25)) # Using the loess (local regression) function, create a smoothed line with a relatively low degree of smoothing (0.25)
# loessFrame <- data.frame(distance = c(seq(-OFFSET_LIMIT + 50, OFFSET_LIMIT - 50, 100)), frequency = loessLine)
# 
# ggplot(data = data[!is.na(data$ratio),], aes(x = ratio)) + 
#   geom_point(aes(y = A), col = '#FF8C00', alpha = 0.10) +
#   geom_line(aes(y = predict(loess(A ~ ratio, data = data[!is.na(data$ratio),], span = 0.25))), col = '#FF8C00') +
#   geom_point(aes(y = C), col = '#7A9500', alpha = 0.10) +
#   geom_line(aes(y = predict(loess(C ~ ratio, data = data[!is.na(data$ratio),], span = 0.25))), col = '#7A9500') +
#   geom_point(aes(y = G), col = '#008A98', alpha = 0.10) +
#   geom_line(aes(y = predict(loess(G ~ ratio, data = data[!is.na(data$ratio),], span = 0.25))), col = '#008A98') +
#   geom_point(aes(y = T), col = '#716EC3', alpha = 0.10) +
#   geom_line(aes(y = predict(loess(T ~ ratio, data = data[!is.na(data$ratio),], span = 0.25))), col = '#716EC3') +
#   ylim(0, 0.025) + xlim(0, 2.5)

binnedData <- binBydNdS(data[complete.cases(data$ratio),])
#- binnedData[(is.na(binnedData))] <- 0
signif(median(binnedData[which(!is.na(binnedData$T) & binnedData$binNum == 'Strongly Purifying'),'T']), 5)
signif(median(binnedData[which(!is.na(binnedData$T) & binnedData$binNum == 'Diversifying'),'T']), 5)

wilTest <- wilcox.test(binnedData[which(!is.na(binnedData$A) & binnedData$binNum == 'Strongly Purifying'),'A'], binnedData[which(!is.na(binnedData$A) & binnedData$binNum == 'Diversifying'),'A'])
paste0("Wilcoxon - V = ", round(wilTest$statistic, 5), " (p-value of ", round(wilTest$p.value, 5), ")")
wilTest <- wilcox.test(binnedData[which(!is.na(binnedData$C) & binnedData$binNum == 'Strongly Purifying'),'C'], binnedData[which(!is.na(binnedData$C) & binnedData$binNum == 'Diversifying'),'C'])
paste0("Wilcoxon - V = ", round(wilTest$statistic, 5), " (p-value of ", round(wilTest$p.value, 5), ")")
wilTest <- wilcox.test(binnedData[which(!is.na(binnedData$G) & binnedData$binNum == 'Strongly Purifying'),'G'], binnedData[which(!is.na(binnedData$G) & binnedData$binNum == 'Diversifying'),'G'])
paste0("Wilcoxon - V = ", round(wilTest$statistic, 5), " (p-value of ", round(wilTest$p.value, 5), ")")
wilTest <- wilcox.test(binnedData[which(!is.na(binnedData$T) & binnedData$binNum == 'Strongly Purifying'),'T'], binnedData[which(!is.na(binnedData$T) & binnedData$binNum == 'Diversifying'),'T'])
paste0("Wilcoxon - V = ", round(wilTest$statistic, 5), " (p-value of ", round(wilTest$p.value, 5), ")")
ggplot(data = binnedData[binnedData$binNum %in% c("Strongly Purifying", "Diversifying"),], aes(x = binNum, y = T, color = as.factor(binNum))) + 
  geom_boxplot(outlier.shape = NA) + ylim(0, 0.01) + 
  labs(title = "Individual Guanine Mutation Rate of MA Genes for the Different Types of Selection (Binned by dN/dS)", x = "dN/dS Ratio Bin (Range)", y = "Mutation Rate per Nucleotide (per base per gene)", "color" = "Selection Type") + 
  scale_x_discrete(labels = c("0.05 <= dN/dS < 0.5", "1.2 <= dN/dS < 2.5")) + 
  scale_colour_manual(values = c('#FF8C00', '#008A98'))


results <- data.frame('nucleotide' = c('A', 'C', 'G', 'T'), 'stronglyPurifying' = c(NA, NA, NA, NA), 'diversifying' = c(NA, NA, NA, NA))
results$stronglyPurifying[1] <- sum(binnedData[which(!is.na(binnedData$A) & binnedData$binNum == 'Strongly Purifying'), 'A'])
results$stronglyPurifying[2] <- sum(binnedData[which(!is.na(binnedData$C) & binnedData$binNum == 'Strongly Purifying'), 'C'])
results$stronglyPurifying[3] <- sum(binnedData[which(!is.na(binnedData$G) & binnedData$binNum == 'Strongly Purifying'), 'G'])
results$stronglyPurifying[4] <- sum(binnedData[which(!is.na(binnedData$T) & binnedData$binNum == 'Strongly Purifying'), 'T'])
results$diversifying[1] <- sum(binnedData[which(!is.na(binnedData$A) & binnedData$binNum == 'Diversifying'), 'A'])
results$diversifying[2] <- sum(binnedData[which(!is.na(binnedData$C) & binnedData$binNum == 'Diversifying'), 'C'])
results$diversifying[3] <- sum(binnedData[which(!is.na(binnedData$G) & binnedData$binNum == 'Diversifying'), 'G'])
results$diversifying[4] <- sum(binnedData[which(!is.na(binnedData$T) & binnedData$binNum == 'Diversifying'), 'T'])

pValues <- c()
for (i in 1:4) { #'[ For each substitution, calculate the pairwise p-values (chi-squared) between the CaeNDR and MA datasets for the strongly purifying bin ]
  chiresults <- chisq.test(rbind(c(results$stronglyPurifying[[i]], results$diversifying[[i]]), c(sum(results$stronglyPurifying) - results$stronglyPurifying[[i]], sum(results$diversifying) - results$diversifying[[i]])))
  pValues <- append(pValues, chiresults$p.value)
}
results$pValue <- c(unlist(p.adjust(pValues, method = 'bonferroni'))) # Add the p-values to the dataframe, adjusting them via the Bonferroni method
ratios <- c()
for (i in 1:4) { #'[ For each substitution, calculate the difference in ratios ((Diversifying / Diversifying_Total) - (Strongly_Purifying / Strongly_Purifying_Total)) between the CaeNDR and MA data, and then append it to the list of ratios ]
  ratios <- append(ratios, (results$diversifying[i] / sum(results$diversifying)) - (results$stronglyPurifying[i] / sum(results$stronglyPurifying)))
}
results$difference <- c(unlist(ratios)) # Add the ratios to the dataframe

results$enriched <- ""
results$enriched <- ifelse(results$difference <= 0, 'Strongly Purifying', 'Diversifying') # Determine which bin the dinucleotide is more strongly enriched in (i.e. in which its count is higher). This uses the ratios, as these are scaled relative to the total number of elements in both bins
# results$enriched <- ifelse(results$pValue <= 0.05, results$enriched, 'None') # Remove the enrichment if the p-value is not significant
results$displayLabel <- ifelse(results$pValue <= 0.05, results$nucleotide, NA) # Give the data point a label only if the p-value is significant. If not, the label is NA (it won't be shown on the plot)
#'[ Plot the data ]
ggplot(results, aes(x = difference, y = -log10(pValue), col = enriched, label = nucleotide)) + 
  geom_point() +
  geom_text_repel() + 
  labs(title = "Mononucleotide Differences Between the 'Diversifying' Bin Compared to the 'Strongly Purifying' Bin (Binned by dN/dS, MA Dataset)", x = "Difference Between Relative Proportions (Diversifying - Strongly Purifying)", y = "P-value (-log10, Bonferroni-Adjusted)", col = "Enrichment") + 
  scale_color_manual(values = c('Diversifying' = '#008A98', 'None' = '#C7C7C7', 'Strongly Purifying' = '#FF8C00'))
# geom_hline(yintercept = -log10(0.05), col = '#474747', linetype = 'dotted')

# # Comparing the sums of each bin by chi-squared shows no significant difference in the mutation rate
# spMutNum <- nrow(binnedData[which(!is.na(binnedData$A) & binnedData$binNum == 'Strongly Purifying'),]) +
#   nrow(binnedData[which(!is.na(binnedData$C) & binnedData$binNum == 'Strongly Purifying'),]) +
#   nrow(binnedData[which(!is.na(binnedData$G) & binnedData$binNum == 'Strongly Purifying'),]) +
#   nrow(binnedData[which(!is.na(binnedData$T) & binnedData$binNum == 'Strongly Purifying'),])
#   
# dMutNum <- nrow(binnedData[which(!is.na(binnedData$A) & binnedData$binNum == 'Diversifying'),]) +
#   nrow(binnedData[which(!is.na(binnedData$C) & binnedData$binNum == 'Diversifying'),]) +
#   nrow(binnedData[which(!is.na(binnedData$G) & binnedData$binNum == 'Diversifying'),]) +
#   nrow(binnedData[which(!is.na(binnedData$T) & binnedData$binNum == 'Diversifying'),])

# meltedResults <- melt(results, id = c('nucleotide', 'pValue', 'difference'))
# meltedResults$mutNum <- c(spMutNum, spMutNum, spMutNum, spMutNum, dMutNum, dMutNum, dMutNum, dMutNum)
# ggplot(meltedResults, aes(x = variable, y = value/mutNum, fill = nucleotide)) + geom_bar(position="dodge", stat = "identity") +
#   labs(title = 'Relative (Summed) Nucleotide-specific Mutation Rates in the Strongly Purifying and Diversifying Genes', x = 'dN/dS Ratio Bin (Range)', y = 'Sum of Nucleotide-specific Mutation Rate / Number of Mutations', fill = 'Nucleotide') +
#   scale_x_discrete(labels = c("0.05 <= dN/dS < 0.5", "1.2 <= dN/dS < 2.5")) + 
#   scale_fill_manual(values = c('#FF8C00', '#7A9500', '#00833C', '#008A98'))





