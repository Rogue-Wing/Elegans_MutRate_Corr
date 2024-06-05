# Takes all the .csv input files for the binned promoters and calculates their overlap

library(pacman)
p_load("ggplot2", "dplyr", "ggVennDiagram", "extrafont", "ggcorrplot", "scales")

# lowpNpS <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Strongly Purifying Coding Promoter Offsets (Nearest, pNpS).csv")
# highpNpS <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Diversifying Coding Promoter Offsets (Nearest, pNpS).csv")
# lowpi <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Extremely Purifying Coding Promoter Offsets (Nearest, pi).csv")
# highpi <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Weakly Purifying Coding Promoter Offsets (Nearest, pi).csv")
# lowdNdS <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Extremely Purifying Coding Promoter Offsets (Nearest, dNdS).csv")
# highdNdS <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/MA Relative Weakly Purifying Coding Promoter Offsets (Nearest, dNdS).csv")

getOverlap <- function(offsets1, offsets2) {
  overlap <- 0
  for (i in 1:length(offsets1)) {
    for (j in 1:length(offsets2)) {
      if (offsets1[i] == offsets2[j]) { # The similarity data is only the genes, not individual mutants. That means that only the gene name needs to be different
        #- print(paste0(offsets1$geneID[i], ' ', offsets2$geneID[j]))
        overlap <- overlap + 1
      }
    }
  }
  return(overlap)
}

# getUniqueMutantCounts <- function(offsets1, offsets2) {
#   overlap <- 0
#   for (i in 1:nrow(offsets1)) {
#     for (j in 1:nrow(offsets2)) {
#       if (offsets1$geneID[i] == offsets2$geneID[j]) { # If the mutant's gene ID and locus match up, then they're the same mutant, and thus are overlapping
#         if (offsets1$locus[i] == offsets2$locus[j]) {
#           #- print(paste0(offsets1$geneID[i], ' ', offsets2$geneID[j]))
#           overlap <- overlap + 1
#         }
#       }
#     }
#   }
#   return(overlap)
# }

promoterSimilarityData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/Binned Promoters Similarity Output.csv")
lowpNpS <- promoterSimilarityData[promoterSimilarityData$lowpNpS == 'True','gene']
highpNpS <- promoterSimilarityData[promoterSimilarityData$highpNpS == 'True','gene']  
lowpi <- promoterSimilarityData[promoterSimilarityData$lowpi == 'True','gene']
highpi <- promoterSimilarityData[promoterSimilarityData$highpi == 'True','gene']
lowdNdS <- promoterSimilarityData[promoterSimilarityData$lowdNdS == 'True','gene']
highdNdS <- promoterSimilarityData[promoterSimilarityData$highdNdS == 'True','gene']

maSimilarityData <- read.csv("C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources/Binned MA Genes Similarity Output.csv")
lowpNpS <- maSimilarityData[maSimilarityData$lowpNpS == 'True','gene']
highpNpS <- maSimilarityData[maSimilarityData$highpNpS == 'True','gene']  
lowpi <- maSimilarityData[maSimilarityData$lowpi == 'True','gene']
highpi <- maSimilarityData[maSimilarityData$highpi == 'True','gene']
lowdNdS <- maSimilarityData[maSimilarityData$lowdNdS == 'True','gene']
highdNdS <- maSimilarityData[maSimilarityData$highdNdS == 'True','gene']

# print(getUniqueGeneCounts(similarityData[similarityData$lowdNdS == 'True',], similarityData[similarityData$highpNpS == 'True',]))

allColours <- c('#c92e62', '#c9832f', '#54850f', '#288879', '#7228a0', '#216392')

plotVennDiagram <- function(binList, binNames, dataset, binTitle, binnedBy, twoTone = F, fullData = promoterSimilarityData) {
  nonGenes <- nrow(fullData[!(fullData$gene %in% unlist(binList)), ])
  
  if (length(binList) == 2) {
    if (twoTone == T) {
      colours <- c('#c92e62', '#216392')
    } else {
      colours <- c('#c92e62', '#a36b27')
    }
    label <- annotate('text', x = 4, y = -3, label = nonGenes)    
  } else if (length(binList) == 3) {
    colours <- c('#c92e62', '#a36b27', '#54850f')
    label <- annotate('text', x = 7.5, y = -5.5, label = nonGenes)
  } else if (length(binList) == 4) {
    colours <- c('#c92e62', '#a36b27', '#54850f', '#288879')
    label <- annotate('text', x = 0.9, y = 0.3, label = nonGenes)
  } else if (length(binList) == 6) {
    colours <- c('#c92e62', '#a36b27', '#54850f', '#288879', '#7228a0', '#216392')
    label <- annotate('text', x = 1000, y = 250, label = nonGenes)
  }

  venn <- ggVennDiagram(binList, color = '#474747',
              category.names = binNames, label = 'count', label_alpha = 0.4,
              set_color = colours) +
            scale_fill_gradient(low = '#9F9F9F', high = '#4F4F4F') +
            theme(text=element_text(size=14,  family="Arial Nova"),
              panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
              legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
              plot.background = element_rect(fill = '#CFCFCF')) +
            labs(title = paste0('Count of Overlaps of ', dataset, ' Genes Between ', binTitle, ' Bins (Binned by ', binnedBy, ')'), fill = 'Gene Count') +
            scale_x_continuous(expand = expansion(mult = .2))
  
  venn <- venn + label
  print(venn)
}

plotVennDiagram(list(highpNpS, highpi, highdNdS), c("pN/pS (High)", "π (High)", "dN/dS (High)"), 'Promoter', 'High', 'pN/pS, π and dN/dS')
plotVennDiagram(list(lowpNpS, lowpi, lowdNdS), c("pN/pS (Low)", "π (Low)", "dN/dS (Low)"), 'Promoter', 'Low', 'pN/pS, π and dN/dS')
plotVennDiagram(list(highpNpS, lowdNdS, highdNdS, lowpi), c("pN/pS (High)", "dN/dS (Low)", "dN/dS (High)", "π (Low)"), 'Promoter', 'Low-High', 'pN/pS, π and dN/dS')
plotVennDiagram(list(lowpNpS, lowdNdS, highdNdS, highpi), c("pN/pS (Low)", "dN/dS (Low)", "dN/dS (High)", "π (High)"), 'Promoter', 'Low-High', 'pN/pS, π and dN/dS')
plotVennDiagram(list(lowpNpS, lowdNdS, lowpi, highpNpS, highdNdS, highpi), c("pN/pS (Low)", "π (Low)", "dN/dS (Low)", "pN/pS (High)", "π (High)", "dN/dS (High)"), 'Promoter', 'All', 'pN/pS, π and dN/dS')


plotVennDiagram(list(highpNpS, highpi, highdNdS), c("pN/pS (High)", "π (High)", "dN/dS (High)"), 'MA', 'High', 'pN/pS, π and dN/dS', F, maSimilarityData)
plotVennDiagram(list(lowpNpS, lowpi, lowdNdS), c("pN/pS (Low)", "π (Low)", "dN/dS (Low)"), 'MA', 'Low', 'pN/pS, π and dN/dS', F, maSimilarityData)
plotVennDiagram(list(highpNpS, lowdNdS, highdNdS, lowpi), c("pN/pS (High)", "dN/dS (Low)", "dN/dS (High)", "π (Low)"), 'MA', 'Low-High', 'pN/pS, π and dN/dS', F, maSimilarityData)
plotVennDiagram(list(lowpNpS, lowdNdS, highdNdS, highpi), c("pN/pS (Low)", "dN/dS (Low)", "dN/dS (High)", "π (High)"), 'MA', 'Low-High', 'pN/pS, π and dN/dS', F, maSimilarityData)
plotVennDiagram(list(lowpNpS, lowdNdS, lowpi, highpNpS, highdNdS, highpi), c("pN/pS (Low)", "π (Low)", "dN/dS (Low)", "pN/pS (High)", "π (High)", "dN/dS (High)"), 'MA', 'All', 'pN/pS, π and dN/dS')



calculateCorrelation <- function(geneNames1, geneNames2, total, pValues, geneCat1 = "x", geneCat2 = "y", silent = T) {
  
  overlapCount <- getOverlap(geneNames1, geneNames2)
  nonGeneCount <- total - length(geneNames1) - length(geneNames2) + overlapCount
  
  v <- ggVennDiagram(list(geneNames1, geneNames2), color = '#474747', category.names = c(geneCat1, geneCat2), label = 'count', label_alpha = 0.4, set_color = c('#FF8C00', '#008A98')) +
    scale_fill_gradient(low = '#575757', high = '#C7C7C7') +
    annotate('text', x = 4.5, y = -2, label = nonGeneCount)
  
  if (silent == F) {
    print(v)
  }
  
  data <- data.frame(
    "x" = c(overlapCount, length(geneNames1) - overlapCount),
    "x'" = c(length(geneNames2) - overlapCount, nonGeneCount),
    row.names = c(paste0(geneCat2), paste0(geneCat2, "'")),
    stringsAsFactors = FALSE
  )
  colnames(data) <- c(paste0(geneCat1), paste0(geneCat1, "'"))
  if (silent == F) {
    print(data)
  }
  
  
  chiTest <- chisq.test(data)
  if (silent == F) {
    print(chiTest)
    print(c(chiTest$observed[1], chiTest$expected[1], chiTest$residuals[1]))
  }
  pValues <- rbind(pValues, c(geneCat1, geneCat2, overlapCount, chiTest$residuals[1], chiTest$p.value, NA))
  
  return(pValues)
}

pNpSCount <- length(lowpNpS) + length(highpNpS)
dNdSCount <- length(lowdNdS) + length(highdNdS)
piCount <- length(lowpi) + length(highpi)

combinations <- list(genes1 = list(lowpNpS, lowpNpS, lowpNpS, lowpNpS, lowpNpS, highpNpS, highpNpS, highpNpS, highpNpS, lowpi, lowpi, lowpi, highpi, highpi, lowdNdS), 
                     genes2 = list(lowdNdS, lowpi, highdNdS, highpi, highpNpS, lowdNdS, lowpi, highdNdS, highpi, lowdNdS, highdNdS, highpi, lowdNdS, highdNdS, highdNdS), 
                     total = c(pNpSCount + dNdSCount, pNpSCount + piCount, pNpSCount + dNdSCount, pNpSCount + piCount, pNpSCount + pNpSCount, 
                               pNpSCount + dNdSCount, pNpSCount + piCount, pNpSCount + dNdSCount, pNpSCount + piCount,
                               piCount + dNdSCount, piCount + dNdSCount, piCount + piCount, piCount + dNdSCount, piCount + dNdSCount, dNdSCount + dNdSCount), 
                     name1 = c("pN/pS (Low)", "pN/pS (Low)", "pN/pS (Low)", "pN/pS (Low)", "pN/pS (Low)", "pN/pS (High)", "pN/pS (High)", 
                               "pN/pS (High)", "pN/pS (High)", "pi (Low)", "pi (Low)", "pi (Low)", "pi (High)", "pi (High)", "dN/dS (Low)"), 
                     name2 = c("dN/dS (Low)", "pi (Low)", "dN/dS (High)", "pi (High)", "pN/pS (High)", "dN/dS (Low)", "pi (Low)", 
                               "dN/dS (High)", "pi (High)", "dN/dS (Low)", "dN/dS (High)", "pi (High)", "dN/dS (Low)", "dN/dS (High)", "dN/dS (High)"))

pValues <- data.frame(pair1 = NA, pair2 = NA, overlap = NA, residual = NA, pValue = NA, adjustedPValue = NA, displayPValue = NA)
for (i in 1:15) {
  pValues <- calculateCorrelation(combinations$genes1[[i]], combinations$genes2[[i]], nrow(promoterSimilarityData), pValues, combinations$name1[i], combinations$name2[i])
}
pValues <- pValues[complete.cases(pValues$pair1),]
pValues$adjustedPValue <- signif(p.adjust(pValues$pValue, 'bonferroni'), 5)
pValues$displayPValue <- pValues$adjustedPValue
pValues$displayPValue <- -log10(pValues$adjustedPValue)
for (i in 1:15) {
  if (pValues$adjustedPValue[i] > 0.05) {
    pValues$overlap[i] <- NA 
  }
  if (pValues$residual[i] < 0) {
    pValues$displayPValue[i] <- pValues$displayPValue[i] * -1
  }
}

ggplot(pValues, aes(x = pair1, y = pair2, fill = displayPValue)) + 
  geom_tile() +
  geom_text(aes(label = overlap)) +
  scale_fill_gradientn(colours=c("#c92e62","#54850f","#4F4F4F", "#7228a0","#216392"), values=rescale(c(-1, -0.05, 0, 0.05, 1)), limit = c(min(pValues$displayPValue), max(pValues$displayPValue)), rescale = ~ scales::rescale_mid(.x, mid = 0)) +
  theme(text=element_text(size=14,  family="Arial Nova"),
        axis.text = element_text(colour = '#000000'),
        panel.grid.major = element_line(colour = '#CFCFCF'), panel.grid.minor = element_line(colour = '#CFCFCF'),
        panel.background = element_rect(fill = '#9F9F9F'), panel.border = element_rect(fill = NA, colour = '#4F4F4F', linewidth = 2),
        legend.background = element_rect(fill = '#CFCFCF', colour = '#4F4F4F', linewidth = 0.6),
        plot.background = element_rect(fill = '#CFCFCF')) +
  labs(title = "Heatmap of the Chi-squared p-values from Pairwise Testing of Bin Overlaps (Promoters)", x = 'Bin', y = 'Bin', 
       fill = "-log10(p-value)\n(Bonferroni-adjusted)", caption = paste0("n = ", pNpSCount + piCount + dNdSCount))
