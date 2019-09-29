## ################### ##
## General functions   ##
## Narinder Singh      ##
## ################### ##

# standard error
plot.karyotype <- function(lines = NULL) {
  # get specific line name and info
  line <- lines$TA[i]; line.info <- lines$notes[i]
  # sample to plot
  sample <- normal.tags.per.bin[, colnames(normal.tags.per.bin) %in% line]
  ploidy.cols <- c('gold', 'tomato', 'lightgreen', 'forestgreen', 'powderblue', 'royalblue')
  col.bin <- rep('white', length(sample))
  col.bin[sample >= 0.25] = ploidy.cols[1]
  col.bin[sample >= 0.75] = ploidy.cols[2]
  col.bin[sample >= 1.25] = ploidy.cols[3]
  col.bin[sample >= 1.75] = ploidy.cols[4]
  col.bin[sample >= 2.25] = ploidy.cols[5]
  col.bin[sample >= 2.75] = ploidy.cols[6]
  # color for tags per bin
  col.bin <- alpha(col.bin, 1)
  
  # plot
  barplot(height = chrom.info$length / 10^6,
          col = 'white', ylim = c(chrom.y.lim.len, 0), axes = F,
          ylab = 'Length in megabase')
  axis(side = 3,
       at = chrom.info$barplot.value, 
       labels = chrom.info$chrom_simple, las = 2, line = 0.05)
  abline(h = c(seq(0, chrom.y.lim.len, 50)), col = 'gray')
  rect(normal.tags.per.bin$barplot.value - spread, normal.tags.per.bin$upper,
       normal.tags.per.bin$barplot.value + spread, normal.tags.per.bin$lower, 
       col = col.bin, border = NA)
  barplot(height = chrom.info$length / 10^6,
          col = NA, ylim = c(chrom.y.lim.len, 0), axes = F, bg = 'lightgray',
          ylab = 'Length in megabase', xlab = paste(line, line.info), add = T)
  axis(side = 2,
       at = seq(0, chrom.y.lim.len, 50),
       labels = seq(0, chrom.y.lim.len, 50), las = 2, line = -0.25)
  
  # place centromeres
  rect(chrom.info$barplot.value - spread, chrom.info$centromere - 1,
       chrom.info$barplot.value + spread, chrom.info$centromere + 1, col = 'black')
  points(chrom.info$barplot.value, chrom.info$centromere, pch = 18, cex = 1.5)
  legend('bottomright', legend = c('0x', '1x', '2x', '3x', '4x', '5x', '6x'), horiz = T, bg = 'white',
         fill = c('white', ploidy.cols), y.intersp = 0.5, x.intersp = 0.5, cex = 1.15)
}
