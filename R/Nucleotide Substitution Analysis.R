library(pacman)
p_load("ggplot2", "reshape2", "ggrepel")

NSNucData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/Summarised NS SNPs v3.csv")
MANucData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/Summarised MA SNPs v3.csv")

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

#'* Bins by the pi value into specifically assigned intervals. Note that it also LOGS THE PI VALUES *
binByPi <- function(dataFrame, zCut = 0) {
  dataFrame$piNorm <- as.numeric(dataFrame$piNorm) # Coerce the pi value to a numeric (rather than chr), because sometimes it just isn't... 
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

#'* Bins by the dN/dS ratio into specifically assigned intervals, adding the binNum to a new column of the dataframe *
binBydNdS <- function(dataFrame) {
  dataFrame <- dataFrame[order(dataFrame[,"dRatio"], decreasing=T),]
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
desireddNdSRatios = c("Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying")


#'* Calculates the relative enrichment of all six possible nucleotide substitutions between the two specified bins (of a specified selection measure), and plots the output *
# NOTE THAT SPACES IN THE lowBin AND highBin MUST BE REPLACED WITH PERIODS (.)
plotSubstitutionEnrichment <- function(nucData, dataset, selectionMeasure, lowBin, highBin) {
  # Bin by the correct selection measure
  if (selectionMeasure == 'pN/pS') {
    nucData <- binBypNpS(nucData)
    bins = desiredpNpSRatios
    n = nrow(nucData)
  } else if (selectionMeasure == 'π') {
    nucData <- binByPi(nucData)
    bins = c("Extremely Purifying", "Strongly Purifying", "Purifying", "Weakly Purifying")
    n = nrow(nucData[which(complete.cases(nucData$piNorm) & nucData$piNorm != 0),])
  } else if (selectionMeasure == 'dN/dS') {
    nucData <- binBydNdS(nucData)
    bins = desireddNdSRatios
    n = nrow(nucData[complete.cases(nucData$dRatio),])
  } else {
    print('ERR: Unknown selection measure')
  }
  
  substitutions <- c("A.G", "G.C", "C.T", "C.A", "T.A", "A.C") # The six substitutions in the nucData
  collatedData <- matrix(NA, nrow = 6, ncol = length(bins)) # Create a matrix to store the sums
  for (i in 1:length(bins)) { # For each of the bins in the nucData
    for (j in 1:6) { # For each of the six substitutions
      collatedData[j, i] <- sum(nucData[which(nucData$binNum == bins[i]),substitutions[j]], na.rm = T) # Sum the total, removing any NA values (that gene doesn't have a substitution of that nucleotide)
    }
  }
  collatedData <- rbind(collatedData, c(sum(collatedData[,1]), sum(collatedData[,2]), sum(collatedData[,3]), sum(collatedData[,4]))) # Sum the totals for each bin
  colnames(collatedData) <- bins # Label the columns
  rownames(collatedData) <- c("AG", "GC", "CT", "CA", "TA", "AC", "Total") # Label the rows
  collatedData <- data.frame(collatedData) # Convert the matrix to a dataframe
  
  pValues <- c()
  for (i in 1:6) { # For each substitution, calculate the pairwise p-values (chi-squared) between the low and high bins
    chiTest <- chisq.test(rbind(c(collatedData[i, lowBin], collatedData[i, highBin]), c(collatedData[7, lowBin] - collatedData[i, lowBin], collatedData[7, highBin] - collatedData[i, highBin])))
    pValues <- append(pValues, chiTest$p.value)
  }
  collatedData$pValue <- c(unlist(p.adjust(pValues, method = 'bonferroni')), NA) # Add the p-values to the dataframe, adjusting them via the Bonferroni method
  ratios <- c()
  for (i in 1:6) { # For each substitution, calculate the difference in ratios ((highBin / highBinTotal) - (lowBin / lowBinTotal)), and then append it to the list of ratios
    ratios <- append(ratios, (collatedData[i, highBin] / collatedData[7, highBin]) - (collatedData[i, lowBin] / collatedData[7, lowBin]))
  }
  collatedData$difference <- c(unlist(ratios), NA) # Add the ratios to the dataframe
  
  print(collatedData)
  
  # Determine the labelling
  collatedData$enriched <- ""
  collatedData$enriched <- ifelse(collatedData$difference <= 0, lowBin, highBin) # Determine which dataset the substitution is more strongly enriched in (i.e. in which its count is higher). This uses the ratios, as these are scaled relative to the total number of elements in both datasets
  collatedData$enriched <- ifelse(collatedData$pValue <= 0.05, collatedData$enriched, 'None') # Remove the enrichment if the p-value is not significant
  for (i in 1:6) {
    collatedData$displayLabel[i] <- ifelse(collatedData$pValue[i] <= 0.05, rownames(collatedData)[i], NA) # Give the data point a label only if the p-value is significant. If not, the label is NA (it won't be shown on the plot)
    if (!is.na(collatedData$displayLabel[i])) {
      if (collatedData$displayLabel[i] == 'AG') {
        collatedData$displayLabel[i] <- 'A:T -> G:C'
      } else if (collatedData$displayLabel[i] == 'GC') {
        collatedData$displayLabel[i] <- 'G:C -> C:G'
      } else if (collatedData$displayLabel[i] == 'CT') {
        collatedData$displayLabel[i] <- 'C:G -> T:A'
      } else if (collatedData$displayLabel[i] == 'CA') {
        collatedData$displayLabel[i] <- 'C:G -> A:T'
      } else if (collatedData$displayLabel[i] == 'TA') {
        collatedData$displayLabel[i] <- 'T:A -> A:T'
      } else if (collatedData$displayLabel[i] == 'AC') {
        collatedData$displayLabel[i] <- 'A:T -> C:G'
      }
    }
  }
  
  colours <- setNames(c('#c92e62', '#4F4F4F', '#216392'), c(gsub(".", " ", lowBin, fixed = T), 'None', gsub(".", " ", highBin, fixed = T)))
  
  # Plot the data 
  enrichmentPlot <- ggplot(collatedData[1:6,], aes(x = difference, y = -log10(pValue), col = gsub(".", " ", enriched, fixed = T), label = displayLabel)) + 
    geom_point() +
    geom_text_repel() +
    labs(title = paste0("Nucleotide Substitution Differences Between the '", gsub(".", " ", highBin, fixed = T), "' and '", gsub(".", " ", lowBin, fixed = T), "' ", selectionMeasure, " Bins (", dataset, " dataset)"), 
      x = paste0("Difference Between Relative Proportions (Proportion of ", gsub(".", " ", highBin, fixed = T), " - Proportion of ", gsub(".", " ", lowBin, fixed = T), ")"), 
      y = "-log10(P-value) (Bonferroni-Adjusted)", col = "Enrichment", caption = paste0("n (genes) = ", n)) +
    scale_color_manual(values = colours) +
    geom_hline(yintercept = -log10(0.05), col = '#2F2F2F', linetype = 'dotted') +
    theme(text=element_text(size=14,  family="Arial Nova"),
      axis.text = element_text(colour = '#000000'),
      panel.grid.major = element_line(linewidth = 1, colour = '#CFCFCF'), panel.grid.minor = element_line(linewidth = 0.2, colour = '#CFCFCF'),  
      panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
      legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
      plot.background = element_rect(fill = '#CFCFCF')) +
    scale_y_continuous(limits = c(-1, 200), minor_breaks = seq(0, 300, 25)) +
    scale_x_continuous(limits = c(-0.055, 0.086))
  
  print(enrichmentPlot)
}

plotSubstitutionEnrichment(NSNucData, 'AG', 'pN/pS', 'Strongly.Purifying', 'Diversifying')
plotSubstitutionEnrichment(MANucData, 'MA', 'pN/pS', 'Strongly.Purifying', 'Diversifying')
plotSubstitutionEnrichment(NSNucData, 'AG', 'π', 'Extremely.Purifying', 'Weakly.Purifying')
plotSubstitutionEnrichment(MANucData, 'MA', 'π', 'Extremely.Purifying', 'Weakly.Purifying')
plotSubstitutionEnrichment(NSNucData, 'AG', 'dN/dS', 'Extremely.Purifying', 'Weakly.Purifying')
plotSubstitutionEnrichment(MANucData, 'MA', 'dN/dS', 'Extremely.Purifying', 'Weakly.Purifying')


#' # Uses either MANucData or NSNucData (i.e. where all six substitution types are separated out)
#' #'* Takes binned nucleotide substitution data and creates a dataframe structure that counts the six different types of substitution for each bin (i.e. it summarises [collates] every individual substitution into single totals) *
#' collateSubs <- function(nucData) {
#'   nucData <- binBypNpS(nucData) # Bin the data by pN/pS
#'   nucData <- data.frame(binNum = nucData$binNum, AG = nucData$A.G, GC = nucData$G.C, CT = nucData$C.T, CA = nucData$C.A, TA = nucData$T.A, AC = nucData$A.C) # Form a new dataframe of just the nucleotide substitution counts
#'   ns = c(sum(nucData[nucData$binNum == "Diversifying", c("AG", "GC", "CT", "CA", "TA", "AC")], na.rm = TRUE), sum(nucData[nucData$binNum == "Neutral", c("AG", "GC", "CT", "CA", "TA", "AC")], na.rm = TRUE), sum(nucData[nucData$binNum == "Weakly Purifying", c("AG", "GC", "CT", "CA", "TA", "AC")], na.rm = TRUE), sum(nucData[nucData$binNum == "Strongly Purifying", c("AG", "GC", "CT", "CA", "TA", "AC")], na.rm = TRUE)) # Get the total counts for each of the four bins (this is to allow for the calculation of ratio'd totals if needed)
#'   meltNucData <- melt(nucData, id=c("binNum")) # Melt the dataframe in order to make it more easily manipulatible
#'   meltNucData$binNum <- factor(meltNucData$binNum, levels = c("Diversifying", "Neutral", "Weakly Purifying", "Strongly Purifying", "Extremely Purifying", "Extremely Diversifying"), ordered = T) # Order the melted dataframe
#'   
#'   subTotals <- c( #Get the total counts for every substitution, ignoring NA values 
#'     sum(meltNucData[which(meltNucData$variable == "AG" & meltNucData$binNum == "Diversifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "GC" & meltNucData$binNum == "Diversifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "CT" & meltNucData$binNum == "Diversifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "CA" & meltNucData$binNum == "Diversifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "TA" & meltNucData$binNum == "Diversifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "AC" & meltNucData$binNum == "Diversifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "AG" & meltNucData$binNum == "Neutral"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "GC" & meltNucData$binNum == "Neutral"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "CT" & meltNucData$binNum == "Neutral"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "CA" & meltNucData$binNum == "Neutral"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "TA" & meltNucData$binNum == "Neutral"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "AC" & meltNucData$binNum == "Neutral"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "AG" & meltNucData$binNum == "Weakly Purifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "GC" & meltNucData$binNum == "Weakly Purifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "CT" & meltNucData$binNum == "Weakly Purifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "CA" & meltNucData$binNum == "Weakly Purifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "TA" & meltNucData$binNum == "Weakly Purifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "AC" & meltNucData$binNum == "Weakly Purifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "AG" & meltNucData$binNum == "Strongly Purifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "GC" & meltNucData$binNum == "Strongly Purifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "CT" & meltNucData$binNum == "Strongly Purifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "CA" & meltNucData$binNum == "Strongly Purifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "TA" & meltNucData$binNum == "Strongly Purifying"),"value"], na.rm = TRUE),
#'     sum(meltNucData[which(meltNucData$variable == "AC" & meltNucData$binNum == "Strongly Purifying"),"value"], na.rm = TRUE)
#'   )
#'   
#'   # Create a dataframe of the data, with the different bins, the substitutions for each bin, and the totals (and ratio'd totals)
#'   collatedSubData <- data.frame(
#'     selection = c(rep("Diversifying", 6), rep("Neutral", 6), rep("Weakly Purifying", 6), rep("Strongly Purifying", 6)),
#'     subs = rep(c("AG", "GC", "CT", "CA", "TA", "AC"), 4),
#'     totals = subTotals,
#'     ratioTotals = subTotals / c(rep(ns[1], 6), rep(ns[2], 6), rep(ns[3], 6), rep(ns[4], 6))
#'   )
#'   collatedSubData$selection <- factor(collatedSubData$selection, levels = c("Strongly Purifying", "Weakly Purifying", "Neutral", "Diversifying")) # Factor the bins (for ordering purposes)
#'   return(collatedSubData)
#' }
#' 
#' collatedCaeNDRSubData <- collateSubs(NSNucData)
#' ggplot(data = collatedCaeNDRSubData, aes(x = selection, y = ratioTotals, fill = subs)) + geom_bar(position="dodge", stat="identity", color = "black") + geom_text(aes(label = sprintf('%0.3f', ratioTotals)), position = position_dodge(width = 0.9), vjust = -0.5) + labs(title = "Relative Proportion of (CaeNDR) Nucleotide Substitutions, Binned by pN/pS Ratio", x = "pN/pS Ratio Bin (Range)", y = "Proportion of Substitutions", fill = "Substitution") + scale_x_discrete(labels = c("0.05 <= pN/pS < 0.5", "0.5 <= pN/pS < 0.8",  "0.8 <= pN/pS < 1.2", "1.2 <= pN/pS < 2.5")) + scale_fill_manual(values = c('#FF8C00', '#7A9500', '#00833C', '#008A98', '#0081C4', '#716EC3'), labels=c("A:T -> C:G","A:T -> G:C", "C:G -> A:T", "C:G -> T:A", "G:C -> C:G", "T:A -> A:T"))
#' collatedMASubData <- collateSubs(MANucData)
#' ggplot(data = collatedMASubData, aes(x = selection, y = ratioTotals, fill = subs)) + geom_bar(position="dodge", stat="identity", color = "black") + geom_text(aes(label = sprintf('%0.3f', ratioTotals)), position = position_dodge(width = 0.9), vjust = -0.5) + labs(title = "Relative Proportion of (MA) Nucleotide Substitutions, Binned by pN/pS Ratio", x = "pN/pS Ratio Bin (Range)", y = "Proportion of Substitutions", fill = "Substitution") + scale_x_discrete(labels = c("0.05 <= pN/pS < 0.5", "0.5 <= pN/pS < 0.8",  "0.8 <= pN/pS < 1.2", "1.2 <= pN/pS < 2.5")) + scale_fill_manual(values = c('#FF8C00', '#7A9500', '#00833C', '#008A98', '#0081C4', '#716EC3'), labels=c("A:T -> C:G","A:T -> G:C", "C:G -> A:T", "C:G -> T:A", "G:C -> C:G", "T:A -> A:T"))
#' 
#' CaeNDRSubData <- data.frame(sub = c('AG', 'GC', 'CT', 'CA', 'TA', 'AC'), stronglyPurifying = collatedCaeNDRSubData[collatedCaeNDRSubData$selection == 'Strongly Purifying','totals'], diversifying = collatedCaeNDRSubData[collatedCaeNDRSubData$selection == 'Diversifying','totals']) # Extract the collated data into a new dataframe for easier manipulation 
#' MASubData <- data.frame(sub = c('AG', 'GC', 'CT', 'CA', 'TA', 'AC'), stronglyPurifying = collatedMASubData[collatedMASubData$selection == 'Strongly Purifying','totals'], diversifying = collatedMASubData[collatedMASubData$selection == 'Diversifying','totals']) # Extract the collated data into a new dataframe for easier manipulation
#' CaeNDRSubData <- rbind(CaeNDRSubData, c(NA, sum(CaeNDRSubData[, 'stronglyPurifying']), sum(CaeNDRSubData[, 'diversifying']))) # Add the totals for the two bins (i.e. all six substitutions)
#' CaeNDRSubData$sub[[7]] <- 'Total' # Add this total to the dataframe
#' MASubData <- rbind(MASubData, c(NA, sum(MASubData[, 'stronglyPurifying']), sum(MASubData[, 'diversifying']))) # Add the totals for the two bins (i.e. all six substitutions)
#' MASubData$sub[[7]] <- 'Total' # Add this total to the dataframe
#' 
#' stronglyPurifyingSubData <- data.frame(sub = c('AG', 'GC', 'CT', 'CA', 'TA', 'AC', 'Total'), caendr = CaeNDRSubData$stronglyPurifying[1:7], ma = MASubData$stronglyPurifying[1:7]) # Take the CaeNDR and MA dataframes and extract the strongly purifying data (into a new dataframe)
#' diversifyingSubData <- data.frame(sub = c('AG', 'GC', 'CT', 'CA', 'TA', 'AC', 'Total'), caendr = CaeNDRSubData$diversifying[1:7], ma = MASubData$diversifying[1:7]) # Take the CaeNDR and MA dataframes and extract the diversifying data (into a new dataframe)
#' 
#' pValues <- c()
#' for (i in 1:6) { # For each substitution, calculate the pairwise p-values (chi-squared) between the CaeNDR and MA datasets for the strongly purifying bin
#'   chiTest <- chisq.test(rbind(c(stronglyPurifyingSubData$caendr[[i]], stronglyPurifyingSubData$ma[[i]]), c(stronglyPurifyingSubData$caendr[[7]] - stronglyPurifyingSubData$caendr[[i]], stronglyPurifyingSubData$ma[[7]] - stronglyPurifyingSubData$ma[[i]])))
#'   pValues <- append(pValues, chiTest$p.value)
#' }
#' stronglyPurifyingSubData$pValue <- c(unlist(p.adjust(pValues, method = 'bonferroni')), NA) # Add the p-values to the dataframe, adjusting them via the Bonferroni method
#' ratios <- c()
#' for (i in 1:6) { # For each substitution, calculate the difference in ratios ((CaeNDR / CaeNDR_Total) - (MA / MA_Total)) between the CaeNDR and MA data, and then append it to the list of ratios
#'   ratios <- append(ratios, (stronglyPurifyingSubData$caendr[i] / stronglyPurifyingSubData$caendr[7]) - (stronglyPurifyingSubData$ma[i] / stronglyPurifyingSubData$ma[7]))
#' }
#' stronglyPurifyingSubData$difference <- c(unlist(ratios), NA) # Add the ratios to the dataframe
#' 
#' pValues <- c()
#' for (i in 1:6) { # For each substitution, calculate the pairwise p-values (chi-squared) between the CaeNDR and MA datasets for the strongly purifying bin
#'   chiTest <- chisq.test(rbind(c(diversifyingSubData$caendr[[i]], diversifyingSubData$ma[[i]]), c(diversifyingSubData$caendr[[7]] - diversifyingSubData$caendr[[i]], diversifyingSubData$ma[[7]] - diversifyingSubData$ma[[i]])))
#'   pValues <- append(pValues, chiTest$p.value)
#' }
#' diversifyingSubData$pValue <- c(unlist(p.adjust(pValues, method = 'bonferroni')), NA) # Add the p-values to the dataframe, adjusting them via the Bonferroni method
#' ratios <- c()
#' for (i in 1:6) { # For each substitution, calculate the difference in ratios ((CaeNDR / CaeNDR_Total) - (MA / MA_Total)) between the CaeNDR and MA data, and then append it to the list of ratios
#'   ratios <- append(ratios, (diversifyingSubData$caendr[i] / diversifyingSubData$caendr[7]) - (diversifyingSubData$ma[i] / diversifyingSubData$ma[7]))
#' }
#' diversifyingSubData$difference <- c(unlist(ratios), NA) # Add the ratios to the dataframe
#' 
#' stronglyPurifyingSubData$enriched <- ""
#' stronglyPurifyingSubData$enriched <- ifelse(stronglyPurifyingSubData$difference <= 0, 'MA', 'CaeNDR') # Determine which dataset the substitution is more strongly enriched in (i.e. in which its count is higher). This uses the ratios, as these are scaled relative to the total number of elements in both datasets
#' stronglyPurifyingSubData$enriched <- ifelse(stronglyPurifyingSubData$pValue <= 0.05, stronglyPurifyingSubData$enriched, 'None') # Remove the enrichment if the p-value is not significant
#' stronglyPurifyingSubData$displayLabel <- ifelse(stronglyPurifyingSubData$pValue <= 0.05, stronglyPurifyingSubData$sub, NA) # Give the data point a label only if the p-value is significant. If not, the label is NA (it won't be shown on the plot)
#' # Plot the data
#' ggplot(stronglyPurifyingSubData, aes(x = difference, y = -log10(pValue), col = enriched, label = displayLabel)) + 
#'   geom_point() +
#'   geom_text_repel() +
#'   labs(title = "Nucleotide Substitution Differences Between the Polymorphism (CaeNDR) and Mutation Accumulation (MA) Datasets (in the 'Strongly Purifying' pN/pS Bin)", x = "Difference Between Relative Proportions (Proportion of CaeNDR - Proportion of MA)", y = "P-value (-log10, Bonferroni-Adjusted)", col = "Enrichment") + 
#'   scale_color_manual(values = c('CaeNDR' = '#008A98', 'None' = '#C7C7C7', 'MA' = '#FF8C00')) +
#'   geom_hline(yintercept = -log10(0.05), col = '#474747', linetype = 'dotted')
#' 
#' diversifyingSubData$enriched <- ""
#' diversifyingSubData$enriched <- ifelse(diversifyingSubData$difference <= 0, 'MA', 'CaeNDR') # Determine which dataset the substitution is more strongly enriched in (i.e. in which its count is higher). This uses the ratios, as these are scaled relative to the total number of elements in both datasets
#' diversifyingSubData$enriched <- ifelse(diversifyingSubData$pValue <= 0.05, diversifyingSubData$enriched, 'None') # Remove the enrichment if the p-value is not significant
#' diversifyingSubData$displayLabel <- ifelse(diversifyingSubData$pValue <= 0.05, diversifyingSubData$sub, NA) # Give the data point a label only if the p-value is significant. If not, the label is NA (it won't be shown on the plot)
#' # Plot the data
#' ggplot(diversifyingSubData, aes(x = difference, y = -log10(pValue), col = enriched, label = displayLabel)) + 
#'   geom_point() +
#'   geom_text_repel() +
#'   labs(title = "Nucleotide Substitution Differences Between the Polymorphism (CaeNDR) and Mutation Accumulation (MA) Datasets (in the 'Diversifying' pN/pS Bin)", x = "Difference Between Relative Proportions (Proportion of CaeNDR - Proportion of MA)", y = "P-value (-log10, Bonferroni-Adjusted)", col = "Enrichment") + 
#'   scale_color_manual(values = c('CaeNDR' = '#008A98', 'None' = '#C7C7C7', 'MA' = '#FF8C00')) +
#'   geom_hline(yintercept = -log10(0.05), col = '#474747', linetype = 'dotted')
#' 
#' 
#' pValues <- c()
#' for (i in 1:6) { # For each substitution, calculate the pairwise p-values (chi-squared) between the strongly purifying and diversifying bins for the CaeNDR dataset
#'   chiTest <- chisq.test(rbind(c(CaeNDRSubData$stronglyPurifying[[i]], CaeNDRSubData$diversifying[[i]]), c(CaeNDRSubData$stronglyPurifying[[7]] - CaeNDRSubData$stronglyPurifying[[i]], CaeNDRSubData$diversifying[[7]] - CaeNDRSubData$diversifying[[i]])))
#'   pValues <- append(pValues, chiTest$p.value)
#' }
#' CaeNDRSubData$pValue <- c(unlist(p.adjust(pValues, method = 'bonferroni')), NA) # Add the p-values to the dataframe, adjusting them via the Bonferroni method
#' ratios <- c()
#' for (i in 1:6) { # For each substitution, calculate the difference in ratios ((d / d_Total) - (sp / sp_Total)) between the strongly purifying (sp) and diversifying (d) bins in the CaeNDR data, and then append it to the list of ratios
#'   ratios <- append(ratios, (CaeNDRSubData$diversifying[i] / CaeNDRSubData$diversifying[7]) - (CaeNDRSubData$stronglyPurifying[i] / CaeNDRSubData$stronglyPurifying[7]))
#' }
#' CaeNDRSubData$difference <- c(unlist(ratios), NA) # Add the ratios to the dataframe
#' 
#' pValues <- c()
#' for (i in 1:6) { # For each substitution, calculate the pairwise p-values (chi-squared) between the strongly purifying and diversifying bins for the MA dataset
#'   chiTest <- chisq.test(rbind(c(MASubData$stronglyPurifying[[i]], MASubData$diversifying[[i]]), c(MASubData$stronglyPurifying[[7]] - MASubData$stronglyPurifying[[i]], MASubData$diversifying[[7]] - MASubData$diversifying[[i]])))
#'   pValues <- append(pValues, chiTest$p.value)
#' }
#' MASubData$pValue <- c(unlist(p.adjust(pValues, method = 'bonferroni')), NA) # Add the p-values to the dataframe, adjusting them via the Bonferroni method
#' ratios <- c()
#' for (i in 1:6) { # For each substitution, calculate the difference in ratios ((d / d_Total) - (sp / sp_Total)) between the strongly purifying (sp) and diversifying (d) bins in the MA data, and then append it to the list of ratios
#'   ratios <- append(ratios, (MASubData$diversifying[i] / MASubData$diversifying[7]) - (MASubData$stronglyPurifying[i] / MASubData$stronglyPurifying[7]))
#' }
#' MASubData$difference <- c(unlist(ratios), NA) # Add the ratios to the dataframe
#' 
#' CaeNDRSubData$enriched <- ""
#' CaeNDRSubData$enriched <- ifelse(CaeNDRSubData$difference <= 0, 'Strongly Purifying', 'Diversifying') # Determine which dataset the substitution is more strongly enriched in (i.e. in which its count is higher). This uses the ratios, as these are scaled relative to the total number of elements in both datasets
#' CaeNDRSubData$enriched <- ifelse(CaeNDRSubData$pValue <= 0.05, CaeNDRSubData$enriched, 'None') # Remove the enrichment if the p-value is not significant
#' CaeNDRSubData$displayLabel <- ifelse(CaeNDRSubData$pValue <= 0.05, CaeNDRSubData$sub, NA) # Give the data point a label only if the p-value is significant. If not, the label is NA (it won't be shown on the plot)
#' for (i in 1:6) {
#'   if (CaeNDRSubData$displayLabel[i] == 'AG') {
#'     CaeNDRSubData$displayLabel[i] <- 'A:T -> G:C'
#'   } else if (CaeNDRSubData$displayLabel[i] == 'GC') {
#'     CaeNDRSubData$displayLabel[i] <- 'G:C -> C:G'
#'   } else if (CaeNDRSubData$displayLabel[i] == 'CT') {
#'     CaeNDRSubData$displayLabel[i] <- 'C:G -> T:A'
#'   } else if (CaeNDRSubData$displayLabel[i] == 'CA') {
#'     CaeNDRSubData$displayLabel[i] <- 'C:G -> A:T'
#'   } else if (CaeNDRSubData$displayLabel[i] == 'TA') {
#'     CaeNDRSubData$displayLabel[i] <- 'T:A -> A:T'
#'   } else if (CaeNDRSubData$displayLabel[i] == 'AC') {
#'     CaeNDRSubData$displayLabel[i] <- 'A:T -> C:G'
#'   } 
#' }
#' # Plot the data 
#' ggplot(CaeNDRSubData, aes(x = difference, y = -log10(pValue), col = enriched, label = displayLabel)) + 
#'   geom_point() +
#'   geom_text_repel() +
#'   labs(title = "Nucleotide Substitution Differences Between the 'Diversifying' (1.2 <= pN/pS < 2.5) and 'Strongly Purifying' (0.05 <= pN/pS < 0.5) Bins of the CaeNDR Dataset", x = "Difference Between Relative Proportions (Proportion of Diversifying - Proportion of Strongly Purifying)", y = "P-value (-log10, Bonferroni-Adjusted)", col = "Enrichment") + 
#'   scale_color_manual(values = c('Diversifying' = '#008A98', 'None' = '#C7C7C7', 'Strongly Purifying' = '#FF8C00')) +
#'   geom_hline(yintercept = -log10(0.05), col = '#474747', linetype = 'dotted') +
#'   geom_vline(xintercept = 0, col = '#474747', linetype='dotted')
  # coord_cartesian(ylim = c(0, max(-log10(CaeNDRSubData$pValue) + (-log10(CaeNDRSubData$pValue) / 10))), clip="off") +
  # annotate('text', x = max(CaeNDRSubData$difference), y = -log10(0.05), label = 'Significance Line (p < 0.05)')

# pValues <- c()
# for (i in 1:6) {
#   chiTest <- chisq.test(rbind(c(MASubData$diversifying[[i]], MASubData$stronglyPurifying[[i]]), c(MASubData$diversifying[[7]] - MASubData$diversifying[[i]], MASubData$stronglyPurifying[[7]] - MASubData$stronglyPurifying[[i]])))
#   pValues <- append(pValues, chiTest$p.value)
# }
# MASubData$pValue <- c(unlist(p.adjust(pValues, method = 'bonferroni')),NA)
# 
# allPValues <- c()
# for (i in 1:1000) {
#   collatedCaeNDRSubData <- collateSubs(NSNucData[sample(nrow(NSNucData), 3355),])
#   CaeNDRSubData <- data.frame(sub = c('AG', 'GC', 'CT', 'CA', 'TA', 'AC'), stronglyPurifying = collatedCaeNDRSubData[collatedCaeNDRSubData$selection == 'Strongly Purifying','totals'], diversifying = collatedCaeNDRSubData[collatedCaeNDRSubData$selection == 'Diversifying','totals']) 
#   CaeNDRSubData <- rbind(CaeNDRSubData, c(NA, sum(CaeNDRSubData[, 'stronglyPurifying']), sum(CaeNDRSubData[, 'diversifying'])))
#   CaeNDRSubData$sub[[7]] <- 'Total'
#   pValues <- c()
#   for (i in 1:6) {
#     chiTest <- chisq.test(rbind(c(CaeNDRSubData$diversifying[[i]], CaeNDRSubData$stronglyPurifying[[i]]), c(CaeNDRSubData$diversifying[[7]] - CaeNDRSubData$diversifying[[i]], CaeNDRSubData$stronglyPurifying[[7]] - CaeNDRSubData$stronglyPurifying[[i]])))
#     pValues <- append(pValues, chiTest$p.value)
#   }
#   allPValues <- append(allPValues, pValues)
# }
# CaeNDRPValues <- data.frame(ag = allPValues[seq(1, 6000, 6)], gc = allPValues[seq(2, 6000, 6)],
#                             ct = allPValues[seq(3, 6000, 6)], ca = allPValues[seq(4, 6000, 6)], 
#                             ta = allPValues[seq(5, 6000, 6)], ac = allPValues[seq(6, 6000, 6)])
# ggplot(CaeNDRPValues, col = '#4F4F4F') + geom_histogram(aes(x = ag), binwidth = 0.5, fill = '#FF8C00', alpha = 0.4) + 
#   geom_histogram(aes(x = gc), binwidth = 0.5, fill = '#7A9500', alpha = 0.4) +
#   geom_histogram(aes(x = ct), binwidth = 0.5, fill = '#00833C', alpha = 0.4) +
#   geom_histogram(aes(x = ca), binwidth = 0.5, fill = '#008A98', alpha = 0.4) +
#   geom_histogram(aes(x = ta), binwidth = 0.5, fill = '#0081C4', alpha = 0.4) +
#   geom_histogram(aes(x = ac), binwidth = 0.5, fill = '#716EC3', alpha = 0.4) +
#   scale_x_continuous(trans='log10')
# 
# CaeNDRPValues <- data.frame()
# subNames <- c('AG', 'GC', 'CT', 'CA', 'TA', 'AC')
# for (i in 1:6) {
#   currentPValueCounts <- plyr::count(allPValues[seq(i, 60, 6)])
#   CaeNDRPValues$temp <- currentPValueCounts$x
#   CaeNDRPValues$tempCount <- currentPValueCounts$freq
#   names(CaeNDRPValues)[names(CaeNDRPValues) == 'temp'] <- subNames[[i]]
# }
# ggplot(CaeNDRPValues, aes(x = ag, y = )) + geom_point()