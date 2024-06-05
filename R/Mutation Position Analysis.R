# Takes an input .csv file of relative SNP positions from the gene's TSS and TTS, groups them, and plots them

library(pacman)
p_load("ggplot2", "dplyr", "extrafont")

OFFSET_LIMIT <- 1500 # The number of bp either side of the TSS/TTS that should be considered
RESOLUTION <- 15 # The size of the bins, in bp (with the exception of the two extremes, which will be half the size)

LOW_PROMOTERS <- 4617 # pNpS = 2303, pi = 3471, dN/dS = 1823, Combined = 4617, Unbinned = 1
HIGH_PROMOTERS <- 4975 # pNpS = 5466, pi = 1091, dN/dS = 984, Combined = 4975, Unbinned = 1

data <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Coding Promoter Offsets (Nearest).csv")
# lowData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Strongly Purifying Coding Promoter Offsets (Nearest, -5000, pNpS).csv")
lowpNpSData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Strongly Purifying Coding Promoter Offsets (Nearest, pNpS).csv")
highpNpSData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Diversifying Coding Promoter Offsets (Nearest, pNpS).csv")
lowpiData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Extremely Purifying Coding Promoter Offsets (Nearest, pi).csv")
highpiData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Weakly Purifying Coding Promoter Offsets (Nearest, pi).csv")
lowdNdSData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Extremely Purifying Coding Promoter Offsets (Nearest, dNdS).csv")
highdNdSData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Weakly Purifying Coding Promoter Offsets (Nearest, dNdS).csv")
lowCombinedData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Low Coding Promoter Offsets (Nearest, Combined).csv")
highCombinedData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative High Coding Promoter Offsets (Nearest, Combined).csv")

processData <- function(dataFrame, span = 0.25) {
  dataFrame$distance <- dataFrame$distance * dataFrame$strand # Reverse the distance/displacement values on the negative strand to make them all positive-strand-facing
  cpCounts <- data.frame(table(cut(dataFrame$distance, c(seq(-OFFSET_LIMIT, OFFSET_LIMIT, RESOLUTION)), labels = c(seq(-OFFSET_LIMIT + (RESOLUTION / 2), OFFSET_LIMIT - (RESOLUTION / 2), RESOLUTION))))) # Divide the data set into a number of bins
  cpCounts$Var1 <- (-OFFSET_LIMIT + (RESOLUTION / 2)) + ((as.numeric(cpCounts$Var1) - 1) * RESOLUTION) # Label the bins with their numerical midpoints 
  names(cpCounts) <- c('midpoint', 'frequency') 
  
  loessLine <- predict(loess((frequency) ~ midpoint, data = cpCounts, span = span), se = T) # Using the loess (local regression) function, create a smoothed line with the specified degree of smoothing
  loessLine$distance <- c(seq(-OFFSET_LIMIT + (RESOLUTION / 2), OFFSET_LIMIT - (RESOLUTION / 2), RESOLUTION)) # Take the line data and give each point the correct x-value (i.e. the distance from the coding promoter)

  return(list(counts = cpCounts, loess = loessLine))  
}

#'* Groups the data into bins of 100bp, and then calculates the mutation rate of every gene on a bin-by-bin basis (the number of mutations for a gene / 100) *
processGroupData <- function(dataFrame) {
  dataFrame$distance <- dataFrame$distance * dataFrame$strand # Reverse the distance/displacement values on the negative strand to make them all positive-strand-facing
  for (point in 1:nrow(dataFrame)) { # For every individual mutation/loci in the dataset
    dataFrame$binMid[point] <- floor(dataFrame$distance[point] / 100) * 100 + 50 # Calculates closest bin to the distance of that mutation (effectively the nearest ..50bp value)
  }
  for (bin in seq(-OFFSET_LIMIT + 50, (OFFSET_LIMIT - 50), 100)) { # For every bin midpoint across the range
    uniqueGeneIDs <- c() # Create a vector to store the all the unique IDs (i.e. any repeats are ignored)
    for (geneID in dataFrame[dataFrame$binMid == bin, 'geneID']) { # For every gene ID
      if (!(geneID %in% uniqueGeneIDs)) { # If the gene ID isn't already in the vector
        uniqueGeneIDs <- append(uniqueGeneIDs, geneID) # Add it to the vector
      }
    }
    for (geneID in uniqueGeneIDs) { # For every gene in the vector (remember, this is for ONE BIN)
      mutCount <- count(dataFrame[which(dataFrame$binMid == bin & dataFrame$geneID == geneID),]) # Get the number of mutations for that gene. One instance of that gene ID means there's only one mutation in that bin (count of 1)
      dataFrame[which(dataFrame$binMid == bin & dataFrame$geneID == geneID), 'mutRate'] <- mutCount / 100 # Store the value of the count, divided by the length of the bin (100bp) as the mutation rate
    }
  }
  cpCounts <- data.frame(midpoint = -OFFSET_LIMIT + 50, frequency = mean(dataFrame[dataFrame$binMid == (-OFFSET_LIMIT + 50), 'mutRate'])) # Create the counts dataframe and add the first midpoint (-2950) and mutation rate to it
  for (bin in seq(-OFFSET_LIMIT + 150, (OFFSET_LIMIT - 50), 100)) { # For all the other bins, add the midpoint and mutation rate
    cpCounts <- rbind(cpCounts, c(bin, mean(dataFrame[dataFrame$binMid == bin, 'mutRate'])))
  }
  
  loessLine <- predict(loess((frequency) ~ midpoint, data = cpCounts, span = 0.25)) # Using the loess (local regression) function, create a smoothed line with a relatively low degree of smoothing (0.25)
  loessFrame <- data.frame(distance = c(seq(-OFFSET_LIMIT + 50, OFFSET_LIMIT - 50, 100)), frequency = loessLine) # Take the line data and give each point the correct x-value (i.e. the distance from the coding promoter)
  
  return(list(counts = cpCounts, loess = loessFrame))
}
  
  
#'* Plots a counts and loess combo, with the smooth loess line overlaid over a partially transparent version of the raw counts *
plotCounts <- function(cpCounts, loessFrame, n = 0) {
  # Plot the data
  ggplot(cpCounts, aes(x = midpoint, y = frequency)) + 
    geom_line(col = '#3F3F3F', linewidth=1, alpha = 0.1) + geom_point(aes(col = midpoint), alpha = 0.15) +
    geom_vline(xintercept = 0, col = '#2F2F2F', linetype = 'dashed') +
    geom_ribbon(aes(xmin = -OFFSET_LIMIT, xmax = OFFSET_LIMIT, ymin = (loessFrame$fit - 0.95 * loessFrame$se), ymax = (loessFrame$fit + 0.95 * loessFrame$se)), fill = '#1F1F1F', alpha = 0.3) +
    geom_line(aes(x = loessFrame$distance, y = loessFrame$fit, col = midpoint), linewidth = 1) +
    scale_colour_gradient(high = '#216392', low = '#c92e62') +
    theme(text=element_text(size=14,  family="Arial Nova"),
      axis.text = element_text(colour = '#000000'),
      panel.grid.major = element_line(linewidth = 1, colour = '#CFCFCF'), panel.grid.minor = element_line(linewidth = 0.2, colour = '#CFCFCF'),  
      panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
      legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
      plot.background = element_rect(fill = '#CFCFCF')) +
    labs(title = paste0("MA Mutational Frequency Around the Coding Promoters (Within ", OFFSET_LIMIT, "bp)"), 
      x = paste0("Distance from Coding Promoter (bp, to the nearest ", RESOLUTION, ")"), 
      y = "Number of Mutations (Substitutions)", caption = paste0('n = ', n)) +
    #- annotate('text', x = 0, y = max(cpCounts$frequency) + 5, label = 'Coding Promoter Midpoint') +
    annotate('text', x = -(OFFSET_LIMIT / 4) * 3, y = 85, label = 'Upstream') +
    annotate('text', x = (OFFSET_LIMIT / 4) * 3, y = 85, label = 'Gene Body') +
    scale_y_continuous(limits = c(0, 100), minor_breaks = seq(0, 100, 10)) +
    scale_x_continuous(minor_breaks = seq(-OFFSET_LIMIT, OFFSET_LIMIT, 250)) +
    guides(color = 'none')
}

#'* Takes  the cpCounts, bins them into 100bp windows, then performs Wilcoxon and T-tests on the bins, starting either side of the midpoint and expanding the bins by 100bp each time. Only the statistically significant results are returned *
expandingBinAnalysis <- function(cpCounts, normaliser = 1) { 
  cpCounts$bpBin <- as.numeric(cut(cpCounts$midpoint, c(seq(-OFFSET_LIMIT, OFFSET_LIMIT, 100))))
  
  wilOutputs <- data.frame('statistic' = NA, 'p-value' = NA)
  tOutputs <- data.frame('statistic' = NA, 'p-value' = NA, 'meanDiff'= NA)
  cpCounts$frequency <- cpCounts$frequency / normaliser # Divides all the frequencies (per RESOLUTION bp) by the input value (usually the number of promoters)
  for (range in 0:((OFFSET_LIMIT / 100) - 1)) {
    x <- cpCounts[which(cpCounts$bpBin <= (OFFSET_LIMIT / 100) & cpCounts$bpBin >= (OFFSET_LIMIT / 100) - range), 'frequency']
    y <- cpCounts[which(cpCounts$bpBin >= ((OFFSET_LIMIT / 100) + 1) & cpCounts$bpBin <= ((OFFSET_LIMIT / 100) + 1) + range), 'frequency']
    wilTest <- wilcox.test(x, rev(y), paired = TRUE) # Perform the test. y is reversed so that the pairs match up correctly
    wilOutputs <- rbind(wilOutputs, c(wilTest$statistic[[1]], wilTest$p.value))
    tTest <- t.test(x, rev(y), paired = TRUE) # Perform the test. y is reversed so that the pairs match up correctly
    tOutputs <- rbind(tOutputs, c(tTest$statistic, tTest$p.value, mean(x) - mean(y)))
  }
  wilOutputs$p.value <- p.adjust(wilOutputs$p.value, method = 'bonferroni') # Correct all the p-values
  wilOutputs <- wilOutputs[which(wilOutputs$p.value <= 0.05 & complete.cases(wilOutputs$p.value)),]
  tOutputs$p.value <- p.adjust(tOutputs$p.value, method = 'bonferroni') # Correct all the p-values
  tOutputs <- tOutputs[which(tOutputs$p.value <= 0.05 & complete.cases(tOutputs$p.value)),]
  
  return(list(wilOut = wilOutputs, tOut = tOutputs))
}

#'* Takes  the cpCounts, bins them into 100bp windows, then performs Wilcoxon and T-tests on the bins, starting either side of the midpoint and moving 100bp either direction each time. Only the statistically significant results are returned *
slidingBinAnalysis <- function(cpCounts, normaliser = 1) {
  cpCounts$bpBin <- as.numeric(cut(cpCounts$midpoint, c(seq(-OFFSET_LIMIT, OFFSET_LIMIT, 100))))
  
  wilOutputs <- data.frame('statistic' = NA, 'p-value' = NA, 'lowBinNum' = NA, 'highBinNum' = NA)
  tOutputs <- data.frame('statistic' = NA, 'p-value' = NA, 'lowBinNum' = NA, 'meanDiff'= NA)
  cpCounts$frequency <- cpCounts$frequency / normaliser # Divides all the frequencies (per RESOLUTION bp) by the input value (usually the number of promoters)
  for (range in 0:((OFFSET_LIMIT / 100) - 1)) {
    x <- cpCounts[which(cpCounts$bpBin == (OFFSET_LIMIT / 100) - range), 'frequency']
    y <- cpCounts[which(cpCounts$bpBin == ((OFFSET_LIMIT / 100) + 1) + range), 'frequency']
    wilTest <- wilcox.test(x, rev(y), paired = TRUE) # Perform the test. y is reversed so that the pairs match up correctly
    wilOutputs <- rbind(wilOutputs, c(wilTest$statistic[[1]], wilTest$p.value, (OFFSET_LIMIT / 100) - range, ((OFFSET_LIMIT / 100) + 1) + range))
    tTest <- t.test(x, rev(y), paired = TRUE) # Perform the test. y is reversed so that the pairs match up correctly
    tOutputs <- rbind(tOutputs, c(tTest$statistic[[1]], tTest$p.value, (OFFSET_LIMIT / 100) - range, mean(x) - mean(y)))
  }
  wilOutputs$actualp.value <- p.adjust(wilOutputs$p.value, method = 'bonferroni') # Correct all the p-values
  wilOutputs <- wilOutputs[which(wilOutputs$p.value <= 0.05 & complete.cases(wilOutputs$p.value)),]
  tOutputs$actualp.value <- p.adjust(tOutputs$p.value, method = 'bonferroni') # Correct all the p-values
  tOutputs <- tOutputs[which(tOutputs$p.value <= 0.05 & complete.cases(tOutputs$p.value)),]
  
  return(list(wilOut = wilOutputs, tOut = tOutputs))
}

#'* Considers the 30 bins either side of the midpoint and compares them (so 60 bins in total), outputting the two statistics and p-values *
twoSideAnalysis <- function(cpCounts, normaliser = 1) {
  cpCounts$bpBin <- as.numeric(cut(cpCounts$midpoint, c(seq(-OFFSET_LIMIT, OFFSET_LIMIT, 100))))
  
  cpCounts$frequency <- cpCounts$frequency / normaliser # Divides all the frequencies (per RESOLUTION bp) by the input value (usually the number of promoters) 
  x <- c(sum(cpCounts[cpCounts$bpBin == 1, 'frequency']))
  y <- c(sum(cpCounts[cpCounts$bpBin == ((OFFSET_LIMIT / 100) + 1), 'frequency']))
  for (bin in 1:(OFFSET_LIMIT / 100)) {
    x[bin] <- sum(cpCounts[cpCounts$bpBin == bin, 'frequency'])
    y[bin] <- sum(cpCounts[cpCounts$bpBin == bin + (OFFSET_LIMIT / 100), 'frequency'])
  }
  wilTest <- wilcox.test(x, rev(y), paired = TRUE) # Perform the test. y is reversed so that the pairs match up correctly
  print(paste0("Wilcoxon - V = ", signif(wilTest$statistic, 5), " (p-value of ", signif(wilTest$p.value, 5), ")"))
  tTest <- t.test(x, rev(y), paired = TRUE) # Perform the test. y is reversed so that the pairs match up correctly
  print(paste0("T-Test - T = ", signif(tTest$statistic, 5), " (p-value of ", signif(tTest$p.value, 5), ")"))
}

#'* Takes two sets of data and, for each equivalent bin, compares between the two with wilcoxon and t-testing *
doubleSlidingAnalysis <- function(lowCpCounts, highCpCounts, lowNormaliser = LOW_PROMOTERS, highNormaliser = HIGH_PROMOTERS) {
  lowCpCounts$bpBin <- as.numeric(cut(lowCpCounts$midpoint, c(seq(-OFFSET_LIMIT, OFFSET_LIMIT, 100))))
  highCpCounts$bpBin <- as.numeric(cut(highCpCounts$midpoint, c(seq(-OFFSET_LIMIT, OFFSET_LIMIT, 100))))
  
  wilOutputs <- data.frame('statistic' = NA, 'p-value' = NA)
  tOutputs <- data.frame('statistic' = NA, 'p-value' = NA, difference = NA, binMin = NA)
  lowCpCounts$frequency <- lowCpCounts$frequency / lowNormaliser # Divides all the frequencies (per RESOLUTION bp) by the input value (usually the number of promoters)
  highCpCounts$frequency <- highCpCounts$frequency / highNormaliser
  for (bin in 1:((OFFSET_LIMIT / 100) * 2)) {
    x <- lowCpCounts[lowCpCounts$bpBin == bin, 'frequency']
    y <- highCpCounts[highCpCounts$bpBin == bin, 'frequency']
    wilTest <- wilcox.test(x, y, paired = TRUE) # Perform the test
    wilOutputs <- rbind(wilOutputs, c(wilTest$statistic[[1]], wilTest$p.value))
    tTest <- t.test(x, y, paired = TRUE) # Perform the test
    tOutputs <- rbind(tOutputs, c(tTest$statistic[[1]], tTest$p.value, mean(y) - mean(x), -OFFSET_LIMIT - 100 + (bin * 100)))
  }
  wilOutputs$p.value <- p.adjust(wilOutputs$p.value, method = 'bonferroni') # Correct all the p-values
  wilOutputs <- wilOutputs[which(wilOutputs$p.value <= 0.05 & complete.cases(wilOutputs$p.value)),] # Return only the p-values that are significant 
  tOutputs$p.value <- p.adjust(tOutputs$p.value, method = 'bonferroni') # Correct all the p-values
  tOutputs <- tOutputs[which(tOutputs$p.value <= 0.05 & complete.cases(tOutputs$p.value)),] # Return only the p-values that are significant 
  
  return(list(wilOut = wilOutputs, tOut = tOutputs))
}

#'* As in the doubleSlidingAnalysis above, but this takes data that have been grouped and summed by mutation rate instead *
doubleSlidingGroupAnalysis <- function(lowCpCounts, highCpCounts, lowNormaliser = LOW_PROMOTERS, highNormaliser = HIGH_PROMOTERS) {

  chiOutputs <- data.frame('statistic' = NA, 'p-value' = NA, difference = NA, binMin = NA)
  lowCpCounts$frequency <- lowCpCounts$frequency / lowNormaliser # Divides all the frequencies (per RESOLUTION bp) by the input value (usually the number of promoters)
  highCpCounts$frequency <- highCpCounts$frequency / highNormaliser
  for (bin in seq(-OFFSET_LIMIT + 50, (OFFSET_LIMIT - 50), 100)) {
    x <- lowCpCounts[lowCpCounts$midpoint == bin, 'frequency']
    y <- highCpCounts[highCpCounts$midpoint == bin, 'frequency']
    chiTest <- chisq.test(rbind(c(x, y), c(sum(lowCpCounts$frequency) - x, sum(highCpCounts$frequency) - y)))
    chiOutputs <- rbind(chiOutputs, c(chiTest$statistic[[1]], chiTest$p.value, y- x, bin - 50))
  }
  chiOutputs <- chiOutputs[which(chiOutputs$p.value <= 0.05 & complete.cases(chiOutputs$p.value)),]
  
  return(chiOutputs)
}

#'* Takes two sets of data and, for each, creates 60 grouped bins and compares them via wilcoxon and t-testing *
summarisedAnalysis <- function(lowCpCounts, highCpCounts, lowNormaliser = LOW_PROMOTERS, highNormaliser = HIGH_PROMOTERS) {
  lowCpCounts$bpBin <- as.numeric(cut(lowCpCounts$midpoint, c(seq(-OFFSET_LIMIT, OFFSET_LIMIT, 100))))
  highCpCounts$bpBin <- as.numeric(cut(highCpCounts$midpoint, c(seq(-OFFSET_LIMIT, OFFSET_LIMIT, 100))))
  
  lowCpCounts$frequency <- lowCpCounts$frequency / lowNormaliser # Divides all the frequencies (per RESOLUTION bp) by the input value (usually the number of promoters)
  highCpCounts$frequency <- highCpCounts$frequency / highNormaliser
  x <- c(sum(lowCpCounts[lowCpCounts$bpBin == 1, 'frequency']))
  y <- c(sum(highCpCounts[highCpCounts$bpBin == 1, 'frequency']))
  for (bin in 2:((OFFSET_LIMIT / 100) * 2)) {
    x[bin] <- sum(lowCpCounts[lowCpCounts$bpBin == bin, 'frequency'])
    y[bin] <- sum(highCpCounts[highCpCounts$bpBin == bin, 'frequency'])
  }
  wilTest <- wilcox.test(x, y, paired = F)
  print(paste0("Wilcoxon - V = ", signif(wilTest$statistic, 5), " (p-value of ", signif(wilTest$p.value, 5), ", low median = ", signif(median(x), 5), ", high median = ", signif(median(y), 5),")"))
  tTest <- t.test(x, y, paired = F, var.equal = F)
  print(paste0("T-Test - T = ", signif(tTest$statistic, 5), " (p-value of ", signif(tTest$p.value, 5), ", low μ = ", signif(median(x), 5), ", high μ = ", signif(median(y), 5), ", Δμ = ", signif(mean(x) - mean(y),5), ")"))
}

#'* As in the summarisedAnalysis, but takes data that have been grouped and summed by mutation rate instead *
summarisedGroupAnalysis <- function(lowCpCounts, highCpCounts, lowNormaliser = LOW_PROMOTERS, highNormaliser = HIGH_PROMOTERS) {
  
  lowCpCounts$frequency <- lowCpCounts$frequency / lowNormaliser # Divides all the frequencies (per RESOLUTION bp) by the input value (usually the number of promoters)
  highCpCounts$frequency <- highCpCounts$frequency / highNormaliser
  x <- lowCpCounts[, 'frequency']
  y <- highCpCounts[,'frequency']
  wilTest <- wilcox.test(x, y, paired = TRUE)
  print(paste0("Wilcoxon - V = ", signif(wilTest$statistic, 5), " (p-value of ", signif(wilTest$p.value, 5), ")"))
  tTest <- t.test(x, y, paired = TRUE)
  print(paste0("T-Test - T = ", signif(tTest$statistic, 5), " (p-value of ", signif(tTest$p.value, 5), ", Δμ = ", signif(mean(x) - mean(y),5), ")"))
}

# Takes a lowData and highData, and plots both onto the same plot
comparativePlot <- function(lowData, highData, lowPromoterCount, highPromoterCount, smoothinAmount, dataset, lowBinName, highBinName, selectionMeasure) {
  processedData <- processData(lowData, smoothinAmount)
  lowCpCounts <- processedData$counts
  lowLoessFrame <- processedData$loess
  
  processedData <- processData(highData, smoothinAmount)
  highCpCounts <- processedData$counts
  highLoessFrame <- processedData$loess
  
  summarisedAnalysis(lowCpCounts, highCpCounts, lowPromoterCount, highPromoterCount)
  stats <- doubleSlidingAnalysis(lowCpCounts, highCpCounts, lowPromoterCount, highPromoterCount)
  tOutputs <- stats$tOut
  
  colours <- setNames(c('#c92e62', '#216392'), c(lowBinName, highBinName))
  
  ggplot() + geom_ribbon(aes(x = lowLoessFrame$distance, ymin = (lowLoessFrame$fit - 0.95 * lowLoessFrame$se) / lowPromoterCount, ymax = (lowLoessFrame$fit + 0.95 * lowLoessFrame$se) / lowPromoterCount), fill = '#c92e62', alpha = 0.3) +
    geom_line(aes(x = lowLoessFrame$distance, y = lowLoessFrame$fit / lowPromoterCount, col = lowBinName), linewidth = 1) + 
    geom_ribbon(aes(x = lowLoessFrame$distance, ymin = (highLoessFrame$fit - 0.95 * highLoessFrame$se) / highPromoterCount, ymax = (highLoessFrame$fit + 0.95 * highLoessFrame$se)  / highPromoterCount), fill = '#216392', alpha = 0.3) +
    geom_line(aes(x = highLoessFrame$distance, y = highLoessFrame$fit / highPromoterCount, col = highBinName), linewidth = 1) +
    theme(text=element_text(size=14,  family="Arial Nova"),
      axis.text = element_text(colour = '#000000'),
      panel.grid.major = element_line(linewidth = 1, colour = '#CFCFCF'), panel.grid.minor = element_line(linewidth = 0.2, colour = '#CFCFCF'),  
      panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
      legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
      plot.background = element_rect(fill = '#CFCFCF')) +
    scale_color_manual(values = colours) +
    geom_rect(data = tOutputs, aes(xmin = binMin, xmax = binMin + 100, ymin = -Inf, ymax = Inf, fill = -log10(p.value)), alpha = 0.4) +
    scale_fill_gradient(high = '#7A9500', low = '#3F3F3F') +
    labs(title = paste0(dataset, ' Mutational Count around the ', lowBinName, ' and ', highBinName, ' Coding Promoters (Within ', OFFSET_LIMIT, 'bp, Binned by ', selectionMeasure,')'), 
      x = paste0("Distance from Coding Promoter (bp, to the nearest ", RESOLUTION, ")"), 
      y = "Number of Mutations / Number of Promoters", col = "Promoter Bin", fill = "-log10(p-value)", caption = paste0('n = ', nrow(lowData), " / n = ", nrow(highData))) +
    annotate('text', x = -(OFFSET_LIMIT / 4) * 3, y= 0.0044, label = 'Upstream') +
    annotate('text', x = (OFFSET_LIMIT / 4) * 3, y = 0.0044, label = 'Gene Body') +
    geom_vline(xintercept = 0, col = '#2F2F2F', linetype = 'dashed') +
    scale_y_continuous(limits = c(0.001, 0.00475), minor_breaks = seq(0, 0.01, 0.00025)) +
    scale_x_continuous(minor_breaks = seq(-OFFSET_LIMIT, OFFSET_LIMIT, 250))
}

# comparativePlot(lowpNpSData, highpNpSData, 3661, 10993, 0.5, 'MA', 'Strongly Purifying', 'Diversifying', 'pN/pS') # Used for TSS analysis instead

comparativePlot(lowpNpSData, highpNpSData, 2303, 5466, 0.5, 'MA', 'Strongly Purifying', 'Diversifying', 'pN/pS')
comparativePlot(lowpiData, highpiData, 3471, 1091, 0.5, 'MA', 'Extremely Purifying', 'Weakly Purifying', 'π')
comparativePlot(lowdNdSData, highdNdSData, 1823, 984, 0.5, 'MA', 'Extremely Purifying', 'Weakly Purifying', 'dN/dS')
comparativePlot(lowCombinedData, highCombinedData, 4617, 4975, 0.5, 'MA', 'Low', 'High', 'pN/pS, π and dN/dS')

# comparativePlot(lowCombinedData, highCombinedData, 6425, 9756, 0.5, 'MA', 'Low', 'High', 'pN/pS, π and dN/dS')  # Used for TSS analysis instead

processedData <- processData(data, 0.25)
cpCounts <- processedData$counts
loessFrame <- processedData$loess
plotCounts(cpCounts, loessFrame, nrow(data))

stats <- expandingBinAnalysis(cpCounts)
wilOutputs <- stats$wilOut
tOutputs <- stats$tOut
stats <- slidingBinAnalysis(cpCounts)
wilOutputs <- stats$wilOut
tOutputs <- stats$tOut
twoSideAnalysis(cpCounts)

processedData <- processData(lowData, 0.25)
lowCpCounts <- processedData$counts
lowLoessFrame <- processedData$loess
plotCounts(lowCpCounts, lowLoessFrame)

processedData <- processData(highData, 0.25)
highCpCounts <- processedData$counts
highLoessFrame <- processedData$loess
plotCounts(highCpCounts, highLoessFrame)

comparativePlot(lowData, highData, 2303, 5466, 0.5)


stats <- expandingBinAnalysis(lowCpCounts, LOW_PROMOTERS)
wilOutputs <- stats$wilOut
tOutputs <- stats$tOut
stats <- slidingBinAnalysis(lowCpCounts, LOW_PROMOTERS)
wilOutputs <- stats$wilOut
tOutputs <- stats$tOut
twoSideAnalysis(lowCpCounts, LOW_PROMOTERS)

summarisedAnalysis(lowCpCounts, highCpCounts)
stats <- doubleSlidingAnalysis(lowCpCounts, highCpCounts)
wilOutputs <- stats$wilOut
tOutputs <- stats$tOut

ggplot() + geom_ribbon(aes(x = lowLoessFrame$distance, ymin = (lowLoessFrame$fit - 0.95 * lowLoessFrame$se) / LOW_PROMOTERS, ymax = (lowLoessFrame$fit + 0.95 * lowLoessFrame$se) / LOW_PROMOTERS), fill = '#008A98', alpha = 0.3) +
  geom_line(aes(x = lowLoessFrame$distance, y = lowLoessFrame$fit / LOW_PROMOTERS, col = 'Low'), linewidth = 1) + 
  geom_ribbon(aes(x = lowLoessFrame$distance, ymin = (highLoessFrame$fit - 0.95 * highLoessFrame$se) / HIGH_PROMOTERS, ymax = (highLoessFrame$fit + 0.95 * highLoessFrame$se)  / HIGH_PROMOTERS), fill = '#FF8C00', alpha = 0.3) +
  geom_line(aes(x = highLoessFrame$distance, y = highLoessFrame$fit / HIGH_PROMOTERS, col = 'High'), linewidth = 1) +
  #- geom_line(aes(x = lowCpCounts$midpoint, y = lowCpCounts$frequency / LOW_PROMOTERS, col = '#474747'), linewidth=1, alpha = 0.1) + geom_point(aes(x = lowCpCounts$midpoint, y = lowCpCounts$frequency / LOW_PROMOTERS, col = 'Low'), alpha = 0.15) + 
  #- geom_line(col = '#474747', linewidth=1, alpha = 0.1) + geom_point(aes(col = midpoint), alpha = 0.15) + 
  scale_color_manual(values = c('High' = '#FF8C00', 'Low' = '#008A98')) +
  geom_rect(data = tOutputs, aes(xmin = binMin, xmax = binMin + 100, ymin = -Inf, ymax = Inf, fill = -log10(p.value)), alpha = 0.4) +
  scale_fill_gradient(high = '#7A9500', low = '#3F3F3F') +
  labs(title = 'Comparison of the MA Mutational Count around the High and Low Coding Promoters (Within 3,000bp, Binned by pN/pS, pi and dN/dS)', x = paste0("Distance from Coding Promoter (bp, to the nearest ", RESOLUTION, ")"), y = "Number of Mutations / Number of Promoters", col = "Promoter Bin", fill = "-log10(p-value)") +
  #- annotate('text', x = 0, y = max(lowLoessFrame$frequency  / LOW_PROMOTERS, highLoessFrame$frequency / HIGH_PROMOTERS) + 0.000003, label = 'Coding Promoter Midpoint') +
  annotate('text', x = -(OFFSET_LIMIT / 4) * 3, y = max(lowLoessFrame$fit  / LOW_PROMOTERS, highLoessFrame$fit / HIGH_PROMOTERS) - 0.00000, label = 'Upstream') +
  annotate('text', x = (OFFSET_LIMIT / 4) * 3, y = max(lowLoessFrame$fit  / LOW_PROMOTERS, highLoessFrame$fit / HIGH_PROMOTERS) - 0.00000, label = 'Gene Body') +
  geom_vline(xintercept = 0, col = '#000000', linetype = 'dashed')

