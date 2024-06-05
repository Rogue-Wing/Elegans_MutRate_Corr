#'* This script is designed to take two FASTA files of CDS transcripts, which should come from different bins of the same data (in this case, the strongly purifying and strongly diversifying bins) *
#'* mono, di and trinucleotide counts are then taken and plotted. *

library(pacman) # Load pacman for easier library loading
p_load("ggplot2", "ggrepel", "seqinr", "metap", "reshape2") # Load libraries. ggplot is for graphing, ggrepel is for labelling of points, seqinr allows for fasta loading, and metap allows for sumlog() of p-values. reshape2 allows for the melting of dataframes for visualisation 

# Load the fasta files
# lowCaeNDRBinFasta <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 CaeNDR Strongly Purifying Transcripts.fa")
# highCaeNDRBinFasta <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 CaeNDR Diversifying Transcripts.fa") 
# lowMABinFasta <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 MA Strongly Purifying Transcripts.fa")
# highMABinFasta <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 MA Diversifying Transcripts.fa") 
# sampledHighBinFasta <- sample(highBinFasta, 1744) # Sample the high bin to produce a list of elements of the same size as the low bin

# Load the FASTA files
lowpNpSCaeNDR <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 Strongly Purifying (pNpS) Transcripts (CaeNDR).fa")
highpNpSCaeNDR <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 Diversifying (pNpS) Transcripts (CaeNDR).fa")
lowpNpSMA <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 Strongly Purifying (pNpS) Transcripts (MA).fa")
highpNpSMA <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 Diversifying (pNpS) Transcripts (MA).fa")

lowpiCaeNDR <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 Extremely Purifying (pi) Transcripts (CaeNDR).fa")
highpiCaeNDR <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 Weakly Purifying (pi) Transcripts (CaeNDR).fa")
lowpiMA <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 Extremely Purifying (pi) Transcripts (MA).fa")
highpiMA <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 Weakly Purifying (pi) Transcripts (MA).fa")

lowdNdSCaeNDR <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 Extremely Purifying (dNdS) Transcripts (CaeNDR).fa")
highdNdSCaeNDR <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 Weakly Purifying (dNdS) Transcripts (CaeNDR).fa")
lowdNdSMA <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 Extremely Purifying (dNdS) Transcripts (MA).fa")
highdNdSMA <- read.fasta("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/WS235 Weakly Purifying (dNdS) Transcripts (MA).fa")




# allPValues <- c()
# monoTotals <- data.frame(bin = c('Strongly Purifying', 'Diversifying'),
#                             a = c(0, 0),
#                             c = c(0, 0),
#                             t = c(0, 0),
#                             g = c(0, 0),
#                             total = c(0, 0))
# for (i in 1:10) {
#   sampledLowBinFasta <- sample(lowBinFasta, 1000)
#   sampledHighBinFasta <- sample(highBinFasta, 1000)
#   monoCounts <- data.frame(bin = c('Strongly Purifying', 'Diversifying'), 
#                            a = c(sum(t(sapply(sampledLowBinFasta, count, 1))[, 'a']), sum(t(sapply(sampledHighBinFasta, count, 1))[, 'a'])), 
#                            c = c(sum(t(sapply(sampledLowBinFasta, count, 1))[, 'c']), sum(t(sapply(sampledHighBinFasta, count, 1))[, 'c'])), 
#                            g = c(sum(t(sapply(sampledLowBinFasta, count, 1))[, 'g']), sum(t(sapply(sampledHighBinFasta, count, 1))[, 'g'])), 
#                            t = c(sum(t(sapply(sampledLowBinFasta, count, 1))[, 't']), sum(t(sapply(sampledHighBinFasta, count, 1))[, 't'])),
#                            total = c(sum(t(sapply(sampledLowBinFasta, count, 1))[,]), sum(t(sapply(sampledHighBinFasta, count, 1))[,])))
#   
#   monoTotals$a <- monoTotals$a + monoCounts$a
#   monoTotals$c <- monoTotals$c + monoCounts$c
#   monoTotals$g <- monoTotals$g + monoCounts$g
#   monoTotals$t <- monoTotals$t + monoCounts$t
#   monoTotals$total <- monoTotals$total + monoCounts$total
#   
#   pValues = c()
#   for (i in 1:4) {
#     chiTest <- chisq.test(rbind(c(monoCounts[1, i + 1], monoCounts[2, i + 1]), c(monoCounts$total[[1]] - monoCounts[1, i + 1], monoCounts$total[[2]] - monoCounts[2, i + 1])))
#     pValues <- append(pValues, chiTest$p.value)
#   }
#   print(pValues)
#   allPValues <- append(allPValues, pValues)
# }
# combinedPValues <- c()
# for (i in 1:4) {
#   combinedPValues <- append(combinedPValues, sumlog(p.adjust(allPValues[seq(i, 40, 4)], method = 'bonferroni'))[[3]]) 
# }
# 
# monoNucRatios <- c((monoTotals$a[[2]] / monoTotals$total[[2]]) - (monoTotals$a[[1]] / monoTotals$total[[1]]), (monoTotals$c[[2]] / monoTotals$total[[2]]) - (monoTotals$c[[1]] / monoTotals$total[[1]]), (monoTotals$g[[2]] / monoTotals$total[[2]]) - (monoTotals$g[[1]] / monoTotals$total[[1]]), (monoTotals$t[[2]] / monoTotals$total[[2]]) - (monoTotals$t[[1]] / monoTotals$total[[1]]))

getMonoCounts <- function(lowBinFasta, highBinFasta) {
  # Create a dataframe of the counts
  monoCounts <- data.frame(bin = c('Strongly Purifying', 'Diversifying'), 
                           a = c(sum(t(sapply(lowBinFasta, count, 1))[, 'a']), sum(t(sapply(highBinFasta, count, 1))[, 'a'])), 
                           c = c(sum(t(sapply(lowBinFasta, count, 1))[, 'c']), sum(t(sapply(highBinFasta, count, 1))[, 'c'])), 
                           g = c(sum(t(sapply(lowBinFasta, count, 1))[, 'g']), sum(t(sapply(highBinFasta, count, 1))[, 'g'])), 
                           t = c(sum(t(sapply(lowBinFasta, count, 1))[, 't']), sum(t(sapply(highBinFasta, count, 1))[, 't'])),
                           total = c(sum(t(sapply(lowBinFasta, count, 1))[,]), sum(t(sapply(highBinFasta, count, 1))[,])))
  pValues = c()
  # Get the pairwise p-values for the four different mononucleotides
  for (i in 1:4) {
    chiTest <- chisq.test(rbind(c(monoCounts[1, i + 1], monoCounts[2, i + 1]), c(monoCounts$total[[1]] - monoCounts[1, i + 1], monoCounts$total[[2]] - monoCounts[2, i + 1])))
    pValues <- append(pValues, chiTest$p.value)
  }
  
  monoNucRatios <- c()
  for (i in 1:4) { # For each mononucleotide, calculate  the difference in ratios ((D / D_Total) - (SP / SP_Total)) between the low and high bin, and then append it to the list of ratios
    monoNucRatios <- append(monoNucRatios, (monoCounts[[i + 1]][2] / monoCounts$total[[2]]) - (monoCounts[[i + 1]][1] / monoCounts$total[[1]]))
  }
  pValues <- p.adjust(pValues, method = "bonferroni") # Adjust the p-values (as this is formed from multiple pairwise comparisons) using the bonferroni method
  monoCounts <- rbind(monoCounts, c('P-Values', unlist(pValues), NA)) # Add the p-values to the dataframe
  monoCounts <- rbind(monoCounts, c('Differences', t(unlist(monoNucRatios)), NA)) # Append the ratio differences to the dataframe
  
  return(monoCounts)
} #'* Gets the counts of each mononucleotide (i.e. a, c, g and t) from the FASTA files, before calculating the pairwise p-values between the bins *

plotMonoGraph <- function(lowBinFasta, highBinFasta, dataset, selectionMeasure, lowBinName, highBinName, yMax, xMin, xMax) {
  monoCounts <- getMonoCounts(lowBinFasta, highBinFasta)
  monoCountsFrame <- data.frame(stronglyPurifying = as.integer(t(monoCounts[1, 2:5])), diversifying = as.integer(t(monoCounts[2, 2:5])), pValues = as.numeric(t(monoCounts[3, 2:5])), difference = as.numeric(t(monoCounts[4, 2:5]))) # Create a dataframe of the count data with the first row (the one with the headers) removed
  names(monoCountsFrame) <- c("stronglyPurifying", "diversifying", "pValue", "difference") # Name the columns
  monoCountsFrame$pValue <- ifelse(monoCountsFrame$pValue < 10^(-yMax), 10^(-yMax), monoCountsFrame$pValue)  # Rounds any p-values smaller than 1e-300 to 1e-300 in order to make them plottable (there are a few p-values that are 0)
  monoCountsFrame$enriched <- ""
  monoCountsFrame$label <- names(monoCounts)[2:5] # Label each row with the appropriate dinucleotide
  monoCountsFrame$enriched <- ifelse(monoCountsFrame$difference <= 0, lowBinName, highBinName) # Determine which bin the dinucleotide is more strongly enriched in (i.e. in which its count is higher). This uses the ratios, as these are scaled relative to the total number of elements in both bins
  monoCountsFrame$enriched <- ifelse(monoCountsFrame$pValue <= 0.05, monoCountsFrame$enriched, 'None') # Remove the enrichment if the p-value is not significant
  monoCountsFrame$displayLabel <- ifelse(monoCountsFrame$pValue <= 0.05, monoCountsFrame$label, NA) # Give the data point a label only if the p-value is significant. If not, the label is NA (it won't be shown on the plot)
  
  colours <- setNames(c('#c92e62', '#216392'), c(lowBinName, highBinName))
  
  # Plot the data
  ggplot(monoCountsFrame, aes(x = difference, y = -log10(pValue), col = enriched, label = displayLabel)) +
    theme(text=element_text(size=14,  family="Arial Nova"),
          axis.text = element_text(colour = '#000000'),
          panel.grid.major = element_line(linewidth = 1, colour = '#CFCFCF'), panel.grid.minor = element_line(linewidth = 0.2, colour = '#CFCFCF'),  
          panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
          legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
          plot.background = element_rect(fill = '#CFCFCF')) +
    geom_point() +
    geom_text_repel() + 
    labs(title = paste0("Mononucleotide Differences Between the '", highBinName, "' and '", lowBinName, "' ", selectionMeasure, " Bins (", dataset, ")"), 
      x = paste0("Difference Between Relative Proportions (Proportion of ", highBinName, " - Proportion of ", lowBinName, ")"), 
      y = "-log10(P-value) (Bonferroni-Adjusted)", col = "Enrichment") + 
    scale_color_manual(values = colours) +
    geom_hline(yintercept = -log10(0.05), col = '#2F2F2F', linetype = 'dotted') +
    scale_x_continuous(limits = c(xMin, xMax)) +
    scale_y_continuous(limits = c(-1, yMax + 5), minor_breaks = seq(0, yMax, 25))
}

allColours <- c('#c92e62', '#c9832f', '#54850f', '#288879', '#7228a0', '#216392')

plotMultiMonoGraph <- function(lowBinFastaList, highBinFastaList, dataset, selectionMeasure, lowBinName, highBinName, yMax, xMin, xMax) {
  multiMonoCounts <- matrix(NA, nrow = 0, ncol = 7) # Create a matrix to store all the counts
  colnames(multiMonoCounts) <- c("stronglyPurifying", "diversifying", "pValue", "difference", "enriched", "label", "displayLabel")
  
  # For each low-high pair, calculate their data as in the non-multi version of the function above. Then add that data to the matrix
  for (i in 1:length(lowBinFastaList)) { # For each low-high pair
    monoCounts <- getMonoCounts(lowBinFastaList[[i]], highBinFastaList[[i]])
    monoCountsFrame <- data.frame(stronglyPurifying = as.integer(t(monoCounts[1, 2:5])), diversifying = as.integer(t(monoCounts[2, 2:5])), pValues = as.numeric(t(monoCounts[3, 2:5])), difference = as.numeric(t(monoCounts[4, 2:5]))) # Create a dataframe of the count data with the first row (the one with the headers) removed
    names(monoCountsFrame) <- c("stronglyPurifying", "diversifying", "pValue", "difference") # Name the columns
    monoCountsFrame$pValue <- ifelse(monoCountsFrame$pValue < 10^(-yMax), 10^(-yMax), monoCountsFrame$pValue)  # Rounds any p-values smaller than 1e-300 to 1e-300 in order to make them plottable (there are a few p-values that are 0)
    monoCountsFrame$enriched <- ""
    monoCountsFrame$label <- names(monoCounts)[2:5] # Label each row with the appropriate dinucleotide
    monoCountsFrame$enriched <- ifelse(monoCountsFrame$difference <= 0, paste0(selectionMeasure[i], " ", lowBinName[i]), paste0(selectionMeasure[i], " ", highBinName[i])) # Determine which bin the dinucleotide is more strongly enriched in (i.e. in which its count is higher). This uses the ratios, as these are scaled relative to the total number of elements in both bins
    monoCountsFrame$enriched <- ifelse(monoCountsFrame$pValue <= 0.05, monoCountsFrame$enriched, 'None') # Remove the enrichment if the p-value is not significant
    monoCountsFrame$displayLabel <- ifelse(monoCountsFrame$pValue <= 0.05, monoCountsFrame$label, NA) # Give the data point a label only if the p-value is significant. If not, the label is NA (it won't be shown on the plot)
      
    multiMonoCounts <- rbind(multiMonoCounts, monoCountsFrame)
    
    print(multiMonoCounts)
  }
  
  multiMonoCounts <- data.frame(multiMonoCounts) # Convert to a dataframe for ggplotting
  
  if (length(lowBinFastaList) == 1) {
    legendLabels <- c(paste0(selectionMeasure[1], " ", lowBinName[1]), 'None', paste0(selectionMeasure[1], " ", highBinName[1]))
    colours <- setNames(c('#c92e62', '#4F4F4F', '#216392'), legendLabels)
  } else if (length(lowBinFastaList) == 2) {
    legendLabels <- c(paste0(selectionMeasure[1], " ", lowBinName[1]), paste0(selectionMeasure[2], " ", lowBinName[2]), 'None', paste0(selectionMeasure[1], " ", highBinName[1]), paste0(selectionMeasure[2], " ", highBinName[2]))
    colours <- setNames(c('#c92e62', '#54850f', '#4F4F4F', '#216392', '#7228a0'), legendLabels)
  } else if (length(lowBinFastaList) == 3) {
    legendLabels <- c(paste0(selectionMeasure[1], " ", lowBinName[1]), paste0(selectionMeasure[2], " ", lowBinName[2]), paste0(selectionMeasure[3], " ", lowBinName[3]), 'None', paste0(selectionMeasure[1], " ", highBinName[1]), paste0(selectionMeasure[2], " ", highBinName[2]), paste0(selectionMeasure[3], " ", highBinName[3]))
    colours <- setNames(c('#c92e62', '#54850f', '#c9832f', '#4F4F4F', '#216392', '#7228a0', '#288879'), legendLabels)
  } else {
    print('ERR: Only a maximum of three pairs of bins can be plotted')
  }

  # Plot the data
  ggplot(multiMonoCounts, aes(x = difference, y = -log10(pValue), col = enriched, label = displayLabel)) +
    geom_point() +
    geom_text_repel() +
    theme(text=element_text(size=14,  family="Arial Nova"),
          axis.text = element_text(colour = '#000000'),
          panel.grid.major = element_line(linewidth = 1, colour = '#CFCFCF'), panel.grid.minor = element_line(linewidth = 0.2, colour = '#CFCFCF'),  
          panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
          legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
          plot.background = element_rect(fill = '#CFCFCF')) +
    labs(title = paste0("Combined Mononucleotide Differences Between High and Low Bins (Binned by pN/pS, π and dN/dS, ", dataset, " dataset)"), 
         x = paste0("Difference Between Relative Proportions (Proportion of High Bin - Proportion of Low Bin)"), 
         y = "-log10(P-value) (Bonferroni-Adjusted)", col = "Enrichment") + 
    scale_color_manual(values = colours, breaks = legendLabels) +
    geom_hline(yintercept = -log10(0.05), col = '#2F2F2F', linetype = 'dotted') +
    scale_x_continuous(limits = c(xMin, xMax)) +
    scale_y_continuous(limits = c(-1, yMax + 5), minor_breaks = seq(0, yMax, 25))
  
}

# monoCountsMelted <- melt(monoCounts, id = c('bin', 'total')) # Melt the dataframe so as to allow for the plotting of four bars for each bin onto the same graph
# monoCountsMelted$bin <- factor(monoCountsMelted$bin, levels = c("Strongly Purifying", "Diversifying")) # Set the bins as factors so that the order is fixed (i.e. the strongly purifying bin will always be the left-most set of bars)
# # Plot the data
# ggplot(monoCountsMelted[which(is.na(monoCountsMelted$bin) == F),], aes(x = bin, y = value / total, fill = variable)) + 
#   geom_bar(stat = 'identity', col = '#4F4F4F', position = position_dodge()) +
#   labs(title = "Proportions of CaeNDR Nucleotides in Gene Sequences, Binned by dN/dS Ratio", x = "dN/dS Ratio Bin (Range)", y = "Proportion of Nucleotides", fill = "Nucleotide") +
#   scale_fill_manual(values = c('#FF8C00', '#7A9500', '#00833C', '#008A98')) +
#   scale_x_discrete(labels = c("0.05 <= dN/dS < 0.5", "1.2 <= dN/dS < 2.5"))
# 
# getSampledDiCounts <- function(lowBinFasta, highBinFasta) {
#   allPValues <- c()
#   allDiNucRatios <- c()
#   for (i in 1:10) { # Iterate through 10 times, sampling the higher bin (as it's a lot bigger [5811 compared to 1744]) to get both the P-value and difference for each dinucleotide
#     lowDiCounts <- as.data.frame(t(sapply(lowBinFasta, count, 2)))
#     sampledHighBinFasta <- sample(highBinFasta, 1744)
#     highDiCounts <- as.data.frame(t(sapply(sampledHighBinFasta, count, 2)))
#     diCounts <- data.frame(bin = c('Strongly Purifying', 'Diversifying'))
#     for (i in 1:16) {
#       diCounts$temp <- c(sum(lowDiCounts[[i]]), sum(highDiCounts[[i]]))
#       names(diCounts)[names(diCounts) == 'temp'] <- colnames(lowDiCounts)[[i]]
#     }
#     diCounts$total <- c(sum(lowDiCounts[,]), sum(highDiCounts[,]))
#     
#     pValues = c()
#     for (i in 1:16) {
#       chiTest <- chisq.test(rbind(c(diCounts[1, i + 1], diCounts[2, i + 1]), c(diCounts$total[[1]] - diCounts[1, i + 1], diCounts$total[[2]] - diCounts[2, i + 1])))
#       pValues <- append(pValues, chiTest$p.value)
#     }
#     allPValues <- append(allPValues, pValues)
#     diNucRatios <- c()
#     for (i in 1:16) {
#       diNucRatios <- append(diNucRatios, (diCounts[[i + 1]][2] / diCounts$total[[2]]) - (diCounts[[i + 1]][1] / diCounts$total[[1]]))
#     }
#     allDiNucRatios <- append(allDiNucRatios, diNucRatios)
#   }
#   combinedPValues <- c()
#   for (i in 1:16) { # Combine the P-values together for each dinucleotide via sumlog
#     combinedPValues <- append(combinedPValues, sumlog(p.adjust(allPValues[seq(i, 40, 16)], method = 'bonferroni'))[[3]])
#     combinedPValues <- replace(combinedPValues, is.na(combinedPValues), 0) # If every P-value is 0, sumlog flags a warning, returning NA. Given this means that the overall P-value is 0, this catches the NA and sets it to 0
#   }
#   combinedDiNucRatios <- c()
#   for (i in 1:16) { # Combine the ratio differences between the high and low bins by averaging
#     combinedDiNucRatios <- append(combinedDiNucRatios, mean(allDiNucRatios[seq(i, 40, 16)])) 
#   }
#   pValues <- p.adjust(pValues, method = 'bonferroni')
#   diCounts <- rbind(diCounts, c('P-Values', t(unlist(combinedPValues)), NA))
#   diCounts <- rbind(diCounts, c('Differences', t(unlist(combinedDiNucRatios)), NA))
# }

getDiCounts <- function(lowBinFasta, highBinFasta) {
  lowDiCounts <- as.data.frame(t(sapply(lowBinFasta, count, 2))) # Count all the different unique dinucleotides in the lower bin
  highDiCounts <- as.data.frame(t(sapply(highBinFasta, count, 2))) # Count all the different unique dinucleotides in the upper bin
  diCounts <- data.frame(bin = c('Strongly Purifying', 'Diversifying')) # Create a dataframe  
  for (i in 1:16) { # Add all the dinucleotide counts to the dataframe, before setting their names (i.e. which dinucleotide it is) as the column name
    diCounts$temp <- c(sum(lowDiCounts[[i]]), sum(highDiCounts[[i]]))
    names(diCounts)[names(diCounts) == 'temp'] <- colnames(lowDiCounts)[[i]]
  }
  diCounts$total <- c(sum(lowDiCounts[,]), sum(highDiCounts[,])) # Sum the dincleotide counts (low and high bin) for each dinucleotide
  
  pValues = c()
  for (i in 1:16) { # For each dinucleotide, calculate  the pairwise p-value between the low and high bin, and then append it to the list of p-values
    chiTest <- chisq.test(rbind(c(diCounts[1, i + 1], diCounts[2, i + 1]), c(diCounts$total[[1]] - diCounts[1, i + 1], diCounts$total[[2]] - diCounts[2, i + 1])))
    pValues <- append(pValues, chiTest$p.value)
  }
  diNucRatios <- c()
  for (i in 1:16) { # For each dinucleotide, calculate  the difference in ratios ((D / D_Total) - (SP / SP_Total)) between the low and high bin, and then append it to the list of ratios
    diNucRatios <- append(diNucRatios, (diCounts[[i + 1]][2] / diCounts$total[[2]]) - (diCounts[[i + 1]][1] / diCounts$total[[1]]))
  }
  pValues <- p.adjust(pValues, method = 'bonferroni') # Adjust the p-values by the Bonferroni method
  diCounts <- rbind(diCounts, c('P-Values', t(unlist(pValues)), NA)) # Append the p-values to the dataframe
  diCounts <- rbind(diCounts, c('Differences', t(unlist(diNucRatios)), NA)) # Append the ratio differences to the dataframe
  
  return(diCounts)
} #'* Gets the counts of each dinucleotide (aa, ac, ag e.t.c.) from the FASTA files, before calculating the pairwise p-values between the bins, as well as the ratio difference *

plotDiGraph <- function(lowBinFasta, highBinFasta, dataset, selectionMeasure, lowBinName, highBinName, yMax, xMin, xMax) {
  diCounts <- getDiCounts(lowBinFasta, highBinFasta)
  diCountsFrame <- data.frame(stronglyPurifying = as.integer(t(diCounts[1, 2:17])), diversifying = as.integer(t(diCounts[2, 2:17])), pValues = as.numeric(t(diCounts[3, 2:17])), difference = as.numeric(t(diCounts[4, 2:17]))) # Create a dataframe of the count data with the first row (the one with the headers) removed
  names(diCountsFrame) <- c("stronglyPurifying", "diversifying", "pValue", "difference") # Name the columns
  diCountsFrame$pValue <- ifelse(diCountsFrame$pValue < 10^(-yMax), 10^(-yMax), diCountsFrame$pValue) # Rounds any p-values smaller than 1e-300 to 1e-300 in order to make them plottable (there are a few p-values that are 0)
  diCountsFrame$enriched <- ""
  diCountsFrame$label <- names(diCounts)[2:17] # Label each row with the appropriate dinucleotide
  diCountsFrame$enriched <- ifelse(diCountsFrame$difference <= 0, lowBinName, highBinName) # Determine which bin the dinucleotide is more strongly enriched in (i.e. in which its count is higher). This uses the ratios, as these are scaled relative to the total number of elements in both bins
  diCountsFrame$enriched <- ifelse(diCountsFrame$pValue <= 0.05, diCountsFrame$enriched, 'None') # Remove the enrichment if the p-value is not significant
  diCountsFrame$displayLabel <- ifelse(diCountsFrame$pValue <= 0.05, diCountsFrame$label, NA) # Give the data point a label only if the p-value is significant. If not, the label is NA (it won't be shown on the plot)
  
  colours <- setNames(c('#c92e62', '#216392'), c(lowBinName, highBinName))
  
  # Plot the data
  ggplot(diCountsFrame, aes(x = difference, y = -log10(pValue), col = enriched, label = displayLabel)) +
    theme(text=element_text(size=14,  family="Arial Nova"),
      axis.text = element_text(colour = '#000000'),
      panel.grid.major = element_line(linewidth = 1, colour = '#CFCFCF'), panel.grid.minor = element_line(linewidth = 0.2, colour = '#CFCFCF'),  
      panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
      legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
      plot.background = element_rect(fill = '#CFCFCF')) +
    geom_point() +
    geom_text_repel() + 
    labs(title = paste0("Dinucleotide Differences Between the '", highBinName, "' and '", lowBinName, "' ", selectionMeasure, " Bins (", dataset, ")"), 
      x = paste0("Difference Between Relative Proportions (Proportion of ", highBinName, " - Proportion of ", lowBinName, ")"), 
      y = "-log10(P-value) (Bonferroni-Adjusted)", col = "Enrichment") + 
    scale_color_manual(values = colours) +
    geom_hline(yintercept = -log10(0.05), col = '#2F2F2F', linetype = 'dotted') +
    scale_x_continuous(limits = c(xMin, xMax)) +
    scale_y_continuous(limits = c(-1, yMax + 5), minor_breaks = seq(0, yMax, 25))
}

plotMultiDiGraph <- function(lowBinFastaList, highBinFastaList, dataset, selectionMeasure, lowBinName, highBinName, yMax, xMin, xMax) {
  multiDiCounts <- matrix(NA, nrow = 0, ncol = 7) # Create a matrix to store all the counts
  colnames(multiDiCounts) <- c("stronglyPurifying", "diversifying", "pValue", "difference", "enriched", "label", "displayLabel")
  
  # For each low-high pair, calculate their data as in the non-multi version of the function above. Then add that data to the matrix
  for (i in 1:length(lowBinFastaList)) { # For each low-high pair
    diCounts <- getDiCounts(lowBinFastaList[[i]], highBinFastaList[[i]])
    diCountsFrame <- data.frame(stronglyPurifying = as.integer(t(diCounts[1, 2:17])), diversifying = as.integer(t(diCounts[2, 2:17])), pValues = as.numeric(t(diCounts[3, 2:17])), difference = as.numeric(t(diCounts[4, 2:17]))) # Create a dataframe of the count data with the first row (the one with the headers) removed
    names(diCountsFrame) <- c("stronglyPurifying", "diversifying", "pValue", "difference") # Name the columns
    diCountsFrame$pValue <- ifelse(diCountsFrame$pValue < 10^(-yMax), 10^(-yMax), diCountsFrame$pValue) # Rounds any p-values smaller than 1e-300 to 1e-300 in order to make them plottable (there are a few p-values that are 0)
    diCountsFrame$enriched <- ""
    diCountsFrame$label <- names(diCounts)[2:17] # Label each row with the appropriate dinucleotide
    diCountsFrame$enriched <- ifelse(diCountsFrame$difference <= 0, paste0(selectionMeasure[i], " ", lowBinName[i]), paste0(selectionMeasure[i], " ", highBinName[i])) # Determine which bin the dinucleotide is more strongly enriched in (i.e. in which its count is higher). This uses the ratios, as these are scaled relative to the total number of elements in both bins
    diCountsFrame$enriched <- ifelse(diCountsFrame$pValue <= 0.05, diCountsFrame$enriched, 'None') # Remove the enrichment if the p-value is not significant
    diCountsFrame$displayLabel <- ifelse(diCountsFrame$pValue <= 0.05, diCountsFrame$label, NA) # Give the data point a label only if the p-value is significant. If not, the label is NA (it won't be shown on the plot)
    
    multiDiCounts <- rbind(multiDiCounts, diCountsFrame)
  }
  
  multiDiCounts <- data.frame(multiDiCounts) # Convert to a dataframe for ggplotting
  
  if (length(lowBinFastaList) == 1) {
    legendLabels <- c(paste0(selectionMeasure[1], " ", lowBinName[1]), 'None', paste0(selectionMeasure[1], " ", highBinName[1]))
    colours <- setNames(c('#c92e62', '#4F4F4F', '#216392'), legendLabels)
  } else if (length(lowBinFastaList) == 2) {
    legendLabels <- c(paste0(selectionMeasure[1], " ", lowBinName[1]), paste0(selectionMeasure[2], " ", lowBinName[2]), 'None', paste0(selectionMeasure[1], " ", highBinName[1]), paste0(selectionMeasure[2], " ", highBinName[2]))
    colours <- setNames(c('#c92e62', '#54850f', '#4F4F4F', '#216392', '#7228a0'), legendLabels)
  } else if (length(lowBinFastaList) == 3) {
    legendLabels <- c(paste0(selectionMeasure[1], " ", lowBinName[1]), paste0(selectionMeasure[2], " ", lowBinName[2]), paste0(selectionMeasure[3], " ", lowBinName[3]), 'None', paste0(selectionMeasure[1], " ", highBinName[1]), paste0(selectionMeasure[2], " ", highBinName[2]), paste0(selectionMeasure[3], " ", highBinName[3]))
    colours <- setNames(c('#c92e62', '#54850f', '#c9832f', '#4F4F4F', '#216392', '#7228a0', '#288879'), legendLabels)
  } else {
    print('ERR: Only a maximum of three pairs of bins can be plotted')
  }
  
  # Plot the data
  ggplot(multiDiCounts, aes(x = difference, y = -log10(pValue), col = enriched, label = displayLabel)) +
    geom_point() +
    geom_text_repel() +
    theme(text=element_text(size=14,  family="Arial Nova"),
          axis.text = element_text(colour = '#000000'),
          panel.grid.major = element_line(linewidth = 1, colour = '#CFCFCF'), panel.grid.minor = element_line(linewidth = 0.2, colour = '#CFCFCF'),  
          panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
          legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
          plot.background = element_rect(fill = '#CFCFCF')) +
    labs(title = paste0("Combined Dinucleotide Differences Between High and Low Bins (Binned by pN/pS, π and dN/dS, ", dataset, " dataset)"), 
         x = paste0("Difference Between Relative Proportions (Proportion of High Bin - Proportion of Low Bin)"), 
         y = "-log10(P-value) (Bonferroni-Adjusted)", col = "Enrichment") + 
    scale_color_manual(values = colours, breaks = legendLabels) +
    geom_hline(yintercept = -log10(0.05), col = '#2F2F2F', linetype = 'dotted') +
    scale_x_continuous(limits = c(xMin, xMax)) +
    scale_y_continuous(limits = c(-1, yMax + 5), minor_breaks = seq(0, yMax, 25))
  
}




# diCountsMelted <- melt(diCounts, id = c('bin', 'total')) # Melt the dataframe so as to allow for the plotting of 16 bars for each bin on the same graph
# diCountsMelted$bin <- factor(diCountsMelted$bin, levels = c("Strongly Purifying", "Diversifying")) # Set the bins as factors so that the order is fixed (i.e. the strongly purifying bin will always be the left-most set of bars)
# # Plot the data
# ggplot(diCountsMelted, aes(x = bin, y = value / total, fill = variable)) + 
#   geom_bar(stat = 'identity', col = '#4F4F4F', position = position_dodge()) +
#   labs(title = "Proportions of CaeNDR Dinucleotides in Gene Sequences, Binned by dN/dS Ratio", x = "dN/dS Ratio Bin (Range)", y = "Proportion of Dinucleotides", fill = "Dinucleotide") +
#   scale_fill_manual(values = c("#3F3F00", "#3F7F10", "#3FBF20", "#3FFF30", "#7F3F40", "#7F7F50", "#7FBF60", "#7FFF70", "#BF3F80", "#BF7F90", "#BFBFA0", "#BFFFB0", "#FF3FC0", "#FF7FD0", "#FFBFE0", "#FFFFF0")) +
#   scale_x_discrete(labels = c("0.05 <= dN/dS < 0.5", "1.2 <= dN/dS < 2.5"))

getTriCounts <- function(lowBinFasta, highBinFasta) {
  lowTriCounts <- as.data.frame(t(sapply(lowBinFasta, count, 3))) # Count the trinucleotides in the low bin
  highTriCounts <- as.data.frame(t(sapply(highBinFasta, count, 3))) # Count the trinucleotides in the high bin
  triCounts <- data.frame(bin = c('Strongly Purifying', 'Diversifying'))
  for (i in 1:64) { # For each trinucleotide, take the sum of it and put in into the dataframe, before labelling it with the identity of the trinucleotide
    triCounts$temp <- c(sum(lowTriCounts[[i]]), sum(highTriCounts[[i]]))
    names(triCounts)[names(triCounts) == "temp"] <- colnames(lowTriCounts)[[i]]
  }
  triCounts$total <- c(sum(lowTriCounts[,]), sum(highTriCounts[,])) # Sum the totals of all the trinucleotide counts, and add them to the dataframe
  
  pValues = c()
  for (i in 1:64) { # For each trinucleotide, calculate the pairwise p-value between the low and high bin, then add it to the list of p-values
    chiTest <- chisq.test(rbind(c(triCounts[1, i + 1], triCounts[2, i + 1]), c(triCounts$total[[1]] - triCounts[1, i + 1], triCounts$total[[2]] - triCounts[2, i + 1])))
    pValues <- append(pValues, chiTest$p.value)
  }
  pValues <- p.adjust(pValues, method = "bonferroni")
  triNucRatios <- c()
  for (i in 1:64) { # For each trinucleotide, calculate  the difference in ratios ((D / D_Total) - (SP / SP_Total)) between the low and high bin, and then append it to the list of ratios
    triNucRatios <- append(triNucRatios, (triCounts[[i + 1]][2] / triCounts$total[[2]]) - (triCounts[[i + 1]][1] / triCounts$total[[1]]))
  }
  triCounts <- rbind(triCounts, c('P-Values', t(unlist(pValues)), '-')) # Add the p-values to the dataframe
  triCounts <- rbind(triCounts, c('Differences', t(unlist(triNucRatios)), '-')) # Add the ratio differences to the dataframe
  
  return(triCounts)
} #'* Gets the counts of each trinucleotide (aaa, aac, aag e.t.c.) from the FASTA files, before calculating the pairwise p-values between the bins, as well as the ratio difference *

plotTriGraph <- function(lowBinFasta, highBinFasta, dataset, selectionMeasure, lowBinName, highBinName, yMax, xMin, xMax) {
  triCounts <- getTriCounts(lowBinFasta, highBinFasta)
  triCountsFrame <- data.frame(stronglyPurifying = as.integer(t(triCounts[1, 2:65])), diversifying = as.integer(t(triCounts[2, 2:65])), pValues = as.numeric(t(triCounts[3, 2:65])), difference = as.numeric(t(triCounts[4, 2:65]))) # Create a dataframe of the count data with the first row (the one with the headers) removed
  names(triCountsFrame) <- c("stronglyPurifying", "diversifying", "pValue", "difference") # Name the columns
  triCountsFrame$pValue <- ifelse(triCountsFrame$pValue < 10^(-yMax), 10^(-yMax), triCountsFrame$pValue) # Rounds any p-values smaller than 1e-300 to 1e-300 in order to make them plottable (there are a few p-values that are 0)
  triCountsFrame$enriched <- ""
  triCountsFrame$label <- names(triCounts)[2:65] # Label each row with the appropriate dinucleotide
  triCountsFrame$enriched <- ifelse(triCountsFrame$difference <= 0, lowBinName, highBinName) # Determine which bin the dinucleotide is more strongly enriched in (i.e. in which its count is higher). This uses the ratios, as these are scaled relative to the total number of elements in both bins
  triCountsFrame$enriched <- ifelse(triCountsFrame$pValue <= 0.05, triCountsFrame$enriched, 'None') # Remove the enrichment if the p-value is not significant
  triCountsFrame$displayLabel <- ifelse(triCountsFrame$pValue <= 0.05, triCountsFrame$label, NA) # Give the data point a label only if the p-value is significant. If not, the label is NA (it won't be shown on the plot)
  
  colours <- setNames(c('#c92e62', '#216392'), c(lowBinName, highBinName))
  
  # Plot the data
  ggplot(triCountsFrame, aes(x = difference, y = -log10(pValue), col = enriched, label = displayLabel)) + 
    theme(text=element_text(size=14,  family="Arial Nova"),
      axis.text = element_text(colour = '#000000'),
      panel.grid.major = element_line(linewidth = 1, colour = '#CFCFCF'), panel.grid.minor = element_line(linewidth = 0.2, colour = '#CFCFCF'),  
      panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
      legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
      plot.background = element_rect(fill = '#CFCFCF')) +
    geom_point() +
    geom_text_repel(min.segment.length = unit(0, 'lines')) + 
    labs(title = paste0("Trinucleotide Differences Between the '", highBinName, "' and '", lowBinName, "' ", selectionMeasure, " Bins (", dataset, ")"), 
      x = paste0("Difference Between Relative Proportions (Proportion of ", highBinName, " - Proportion of ", lowBinName, ")"), 
      y = "-log10(P-value) (Bonferroni-Adjusted)", col = "Enrichment") + 
    scale_color_manual(values = colours) +
    geom_hline(yintercept = -log10(0.05), col = '#2F2F2F', linetype = 'dotted') +
    scale_x_continuous(limits = c(xMin, xMax)) +
    scale_y_continuous(limits = c(-1, yMax + 5), minor_breaks = seq(0, yMax, 25))
}

plotMultiTriGraph <- function(lowBinFastaList, highBinFastaList, dataset, selectionMeasure, lowBinName, highBinName, yMax, xMin, xMax) {
  multiTriCounts <- matrix(NA, nrow = 0, ncol = 7) # Create a matrix to store all the counts
  colnames(multiTriCounts) <- c("stronglyPurifying", "diversifying", "pValue", "difference", "enriched", "label", "displayLabel")
  
  # For each low-high pair, calculate their data as in the non-multi version of the function above. Then add that data to the matrix
  for (i in 1:length(lowBinFastaList)) { # For each low-high pair
    triCounts <- getTriCounts(lowBinFastaList[[i]], highBinFastaList[[i]])
    triCountsFrame <- data.frame(stronglyPurifying = as.integer(t(triCounts[1, 2:65])), diversifying = as.integer(t(triCounts[2, 2:65])), pValues = as.numeric(t(triCounts[3, 2:65])), difference = as.numeric(t(triCounts[4, 2:65]))) # Create a dataframe of the count data with the first row (the one with the headers) removed
    names(triCountsFrame) <- c("stronglyPurifying", "diversifying", "pValue", "difference") # Name the columns
    triCountsFrame$pValue <- ifelse(triCountsFrame$pValue < 10^(-yMax), 10^(-yMax), triCountsFrame$pValue) # Rounds any p-values smaller than 1e-300 to 1e-300 in order to make them plottable (there are a few p-values that are 0)
    triCountsFrame$enriched <- ""
    triCountsFrame$label <- names(triCounts)[2:65] # Label each row with the appropriate dinucleotide
    triCountsFrame$enriched <- ifelse(triCountsFrame$difference <= 0, paste0(selectionMeasure[i], " ", lowBinName[i]), paste0(selectionMeasure[i], " ", highBinName[i])) # Determine which bin the dinucleotide is more strongly enriched in (i.e. in which its count is higher). This uses the ratios, as these are scaled relative to the total number of elements in both bins
    triCountsFrame$enriched <- ifelse(triCountsFrame$pValue <= 0.05, triCountsFrame$enriched, 'None') # Remove the enrichment if the p-value is not significant
    triCountsFrame$displayLabel <- ifelse(triCountsFrame$pValue <= 0.05, triCountsFrame$label, NA) # Give the data point a label only if the p-value is significant. If not, the label is NA (it won't be shown on the plot)
    
    multiTriCounts <- rbind(multiTriCounts, triCountsFrame)
  }
  
  multiTriCounts <- data.frame(multiTriCounts) # Convert to a dataframe for ggplotting
  
  if (length(lowBinFastaList) == 1) {
    legendLabels <- c(paste0(selectionMeasure[1], " ", lowBinName[1]), 'None', paste0(selectionMeasure[1], " ", highBinName[1]))
    colours <- setNames(c('#c92e62', '#4F4F4F', '#216392'), legendLabels)
  } else if (length(lowBinFastaList) == 2) {
    legendLabels <- c(paste0(selectionMeasure[1], " ", lowBinName[1]), paste0(selectionMeasure[2], " ", lowBinName[2]), 'None', paste0(selectionMeasure[1], " ", highBinName[1]), paste0(selectionMeasure[2], " ", highBinName[2]))
    colours <- setNames(c('#c92e62', '#54850f', '#4F4F4F', '#216392', '#7228a0'), legendLabels)
  } else if (length(lowBinFastaList) == 3) {
    legendLabels <- c(paste0(selectionMeasure[1], " ", lowBinName[1]), paste0(selectionMeasure[2], " ", lowBinName[2]), paste0(selectionMeasure[3], " ", lowBinName[3]), 'None', paste0(selectionMeasure[1], " ", highBinName[1]), paste0(selectionMeasure[2], " ", highBinName[2]), paste0(selectionMeasure[3], " ", highBinName[3]))
    colours <- setNames(c('#c92e62', '#54850f', '#c9832f', '#4F4F4F', '#216392', '#7228a0', '#288879'), legendLabels)
  } else {
    print('ERR: Only a maximum of three pairs of bins can be plotted')
  }
  
  # Plot the data
  ggplot(multiTriCounts, aes(x = difference, y = -log10(pValue), col = enriched, label = displayLabel)) +
    geom_point() +
    geom_text_repel() +
    theme(text=element_text(size=14,  family="Arial Nova"),
          axis.text = element_text(colour = '#000000'),
          panel.grid.major = element_line(linewidth = 1, colour = '#CFCFCF'), panel.grid.minor = element_line(linewidth = 0.2, colour = '#CFCFCF'),  
          panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
          legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
          plot.background = element_rect(fill = '#CFCFCF')) +
    labs(title = paste0("Combined Trinucleotide Differences Between High and Low Bins (Binned by pN/pS, π and dN/dS, ", dataset, " dataset)"), 
         x = paste0("Difference Between Relative Proportions (Proportion of High Bin - Proportion of Low Bin)"), 
         y = "-log10(P-value) (Bonferroni-Adjusted)", col = "Enrichment") + 
    scale_color_manual(values = colours, breaks = legendLabels) +
    geom_hline(yintercept = -log10(0.05), col = '#2F2F2F', linetype = 'dotted') +
    scale_x_continuous(limits = c(xMin, xMax)) +
    scale_y_continuous(limits = c(-1, yMax + 5), minor_breaks = seq(0, yMax, 25))
  
}

# labs(title = paste0("Combined Trinucleotide Differences Between High and Low Bins (Binned by pN/pS, π and dN/dS, ", dataset, " dataset)"), 
# x = paste0("Difference Between Relative Proportions (Proportion of High Bin - Proportion of Low Bin)"),

plotMonoGraph(lowpNpSCaeNDR, highpNpSCaeNDR, 'AG', 'pN/pS', 'Strongly Purifying', 'Diversifying', 200, -0.008, 0.008)
plotMonoGraph(lowpNpSMA, highpNpSMA, 'MA', 'pN/pS', 'Strongly Purifying', 'Diversifying', 200, -0.008, 0.008)
plotMonoGraph(lowpiCaeNDR, highpiCaeNDR, 'AG', 'π', 'Extremely Purifying', 'Weakly Purifying', 150, -0.0075, 0.0075)
plotMonoGraph(lowpiMA, highpiMA, 'MA', 'π', 'Extremely Purifying', 'Weakly Purifying', 150, -0.0075, 0.0075)
plotMonoGraph(lowdNdSCaeNDR, highdNdSCaeNDR, 'AG', 'dN/dS', 'Extremely Purifying', 'Weakly Purifying', 200, -0.015, 0.0075)
plotMonoGraph(lowdNdSMA, highdNdSMA, 'MA', 'dN/dS', 'Extremely Purifying', 'Weakly Purifying', 200, -0.015, 0.0075)

plotDiGraph(lowpNpSCaeNDR, highpNpSCaeNDR, 'AG', 'pN/pS', 'Strongly Purifying', 'Diversifying', 250, -0.015, 0.0075)
plotDiGraph(lowpNpSMA, highpNpSMA, 'MA', 'pN/pS', 'Strongly Purifying', 'Diversifying', 250, -0.015, 0.0075)
plotDiGraph(lowpiCaeNDR, highpiCaeNDR, 'AG', 'π', 'Extremely Purifying', 'Weakly Purifying', 100, -0.005, 0.005)
plotDiGraph(lowpiMA, highpiMA, 'MA', 'π', 'Extremely Purifying', 'Weakly Purifying', 100, -0.005, 0.005)
plotDiGraph(lowdNdSCaeNDR, highdNdSCaeNDR, 'AG', 'dN/dS', 'Extremely Purifying', 'Weakly Purifying', 100, -0.0125, 0.0075)
plotDiGraph(lowdNdSMA, highdNdSMA, 'MA', 'dN/dS', 'Extremely Purifying', 'Weakly Purifying', 100, -0.0125, 0.0075)

plotTriGraph(lowpNpSCaeNDR, highpNpSCaeNDR, 'AG', 'pN/pS', 'Strongly Purifying', 'Diversifying', 300, -0.015, 0.005)
plotTriGraph(lowpNpSMA, highpNpSMA, 'MA', 'pN/pS', 'Strongly Purifying', 'Diversifying', 300, -0.015, 0.005)
plotTriGraph(lowpiCaeNDR, highpiCaeNDR, 'AG', 'π', 'Extremely Purifying', 'Weakly Purifying', 200, -0.005, 0.005)
plotTriGraph(lowpiMA, highpiMA, 'MA', 'π', 'Extremely Purifying', 'Weakly Purifying', 200, -0.005, 0.005)
plotTriGraph(lowdNdSCaeNDR, highdNdSCaeNDR, 'AG', 'dN/dS', 'Extremely Purifying', 'Weakly Purifying', 100, -0.0075, 0.005)
plotTriGraph(lowdNdSMA, highdNdSMA, 'MA', 'dN/dS', 'Extremely Purifying', 'Weakly Purifying', 100, -0.0075, 0.005)

plotMultiMonoGraph(list(lowpNpSCaeNDR, lowpiCaeNDR, lowdNdSCaeNDR), list(highpNpSCaeNDR, highpiCaeNDR, highdNdSCaeNDR), 'AG', c('pN/pS', 'π', 'dN/dS'), c('Strongly Purifying', 'Extremely Purifying', 'Extremely Purifying'), c('Diversifying', 'Weakly Purifying', 'Weakly Purifying'), 200, -0.015, 0.008)
plotMultiMonoGraph(list(lowpNpSMA, lowpiMA, lowdNdSMA), list(highpNpSMA, highpiMA, highdNdSMA), 'MA', c('pN/pS', 'π', 'dN/dS'), c('Strongly Purifying', 'Extremely Purifying', 'Extremely Purifying'), c('Diversifying', 'Weakly Purifying', 'Weakly Purifying'), 200, -0.015, 0.008)
plotMultiDiGraph(list(lowpNpSCaeNDR, lowpiCaeNDR, lowdNdSCaeNDR), list(highpNpSCaeNDR, highpiCaeNDR, highdNdSCaeNDR), 'AG', c('pN/pS', 'π', 'dN/dS'), c('Strongly Purifying', 'Extremely Purifying', 'Extremely Purifying'), c('Diversifying', 'Weakly Purifying', 'Weakly Purifying'), 200, -0.015, 0.0075)
plotMultiDiGraph(list(lowpNpSMA, lowpiMA, lowdNdSMA), list(highpNpSMA, highpiMA, highdNdSMA), 'MA', c('pN/pS', 'π', 'dN/dS'), c('Strongly Purifying', 'Extremely Purifying', 'Extremely Purifying'), c('Diversifying', 'Weakly Purifying', 'Weakly Purifying'), 200, -0.015, 0.0075)
plotMultiTriGraph(list(lowpNpSCaeNDR, lowpiCaeNDR, lowdNdSCaeNDR), list(highpNpSCaeNDR, highpiCaeNDR, highdNdSCaeNDR), 'AG', c('pN/pS', 'π', 'dN/dS'), c('Strongly Purifying', 'Extremely Purifying', 'Extremely Purifying'), c('Diversifying', 'Weakly Purifying', 'Weakly Purifying'), 200, -0.015, 0.008)
plotMultiTriGraph(list(lowpNpSMA, lowpiMA, lowdNdSMA), list(highpNpSMA, highpiMA, highdNdSMA), 'MA', c('pN/pS', 'π', 'dN/dS'), c('Strongly Purifying', 'Extremely Purifying', 'Extremely Purifying'), c('Diversifying', 'Weakly Purifying', 'Weakly Purifying'), 200, -0.015, 0.005)

# triCountsMelted <- melt(triCounts, id = c('bin', 'total')) # Melt the dataframe so as to allow for the plotting of 64 bars for each bin on the same graph
# triCountsMelted$bin <- factor(triCountsMelted$bin, levels = c("Strongly Purifying", "Diversifying")) # Set the bins as factors so that the order is fixed (i.e. the strongly purifying bin will always be the left-most set of bars)
# # Plot the data
# ggplot(triCountsMelted, aes(x = bin, y = value / total, fill = variable)) + 
#   geom_bar(stat = 'identity', col = '#4F4F4F', position = position_dodge()) +
#   labs(title = "Proportions of CaeNDR Trinucleotides in Gene Sequences, Binned by dN/dS Ratio", x = "dN/dS Ratio Bin (Range)", y = "Proportion of Trinucleotides", fill = "Trinucleotide") +
#   scale_fill_manual(values = c("#3F3F3F", "#3F3F7F", "#3F3FBF", "#3F3FFF", "#3F7F3F", "#3F7F7F", "#3F7FBF", "#3F7FFF", "#3FBF3F", "#3FBF7F", "#3FBFBF", "#3FBFFF", "#3FFF3F", "#3FFF7F", "#3FFFBF", "#3FFFFF", "#7F3F3F", "#7F3F7F", "#7F3FBF", "#7F3FFF", "#7F7F3F", "#7F7F7F", "#7F7FBF", "#7F7FFF", "#7FBF3F", "#7FBF7F", "#7FBFBF", "#7FBFFF", "#7FFF3F", "#7FFF7F", "#7FFFBF", "#7FFFFF", "#BF3F3F", "#BF3F7F", "#BF3FBF", "#BF3FFF", "#BF7F3F", "#BF7F7F", "#BF7FBF", "#BF7FFF", "#BFBF3F", "#BFBF7F", "#BFBFBF", "#BFBFFF", "#BFFF3F", "#BFFF7F", "#BFFFBF", "#BFFFFF", "#FF3F3F", "#FF3F7F", "#FF3FBF", "#FF3FFF", "#FF7F3F", "#FF7F7F", "#FF7FBF", "#FF7FFF", "#FFBF3F", "#FFBF7F", "#FFBFBF", "#FFBFFF", "#FFFF3F", "#FFFF7F", "#FFFFBF", "#FFFFFF")) +
#   scale_x_discrete(labels = c("0.05 <= dN/dS < 0.5", "1.2 <= dN/dS < 2.5"))

diMACounts <- getTriCounts(lowMABinFasta, highMABinFasta)
diCaeNDRCounts <- getTriCounts(lowCaeNDRBinFasta, highCaeNDRBinFasta)

differenceCounts <- c()
differenceSeqs <- c()
for (i in 1:64) {
  if ((diCaeNDRCounts[4, i + 1] < 0) & (diMACounts[4, i + 1] > 0)) {
    differenceCounts[length(differenceCounts) + 1] <- as.numeric(diCaeNDRCounts[4, i + 1]) - as.numeric(diMACounts[4, i + 1])
    differenceSeqs[length(differenceSeqs) + 1] <- names(diCaeNDRCounts[i + 1])
    print(i + 1)
  } else if ((diCaeNDRCounts[4, i + 1] > 0) & (diMACounts[4, i + 1] < 0)) {
    differenceCounts[length(differenceCounts) + 1] <- as.numeric(diCaeNDRCounts[4, i + 1]) - as.numeric(diMACounts[4, i + 1])
    differenceSeqs[length(differenceSeqs) + 1] <- names(diCaeNDRCounts[i + 1])
    print(i + 1)
  }
}


