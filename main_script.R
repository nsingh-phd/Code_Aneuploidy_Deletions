########################################
## Chromosome aberration analysis     ##
## Aberration - chromosomal deletions ##
## Narinder Singh                     ##
########################################

# install/load pacman for package management
if(!require(pacman)) install.packages(pacman); require(pacman)

# install/load required packages
p_load(data.table, ggplot2, dplyr)

# import functions
source('plot_function.R')

####################
# read data files ##
####################

# tags by taxa distribution file
tagDist.path <- '../Manuscript_Aneuploidy_Deletions/Manuscript_Aneuploidy_Deletions_TagTaxaDist.txt'
# unique sam file
sam.path <- '../Manuscript_Aneuploidy_Deletions/unique.sam'
# pdf outfile basename
pdf.outfile <- 'WGRC_deletions'
# read chromosome size, centromere etc. info
chrom.info <- read.table('chrom_info.txt', header = T, as.is = T)
chrom.info <- cbind(chrom.info, 'barplot.value' = barplot(chrom.info$length/1000000, plot = F))

# need chromosome size correction? T or F (logical)
need.chrSizeCorrection <- F
# analyze per rep (T) or combined (F)
analyze.per.rep <- T

# read in tag distribution and sam file from TASSEL5 pipeline
tagDist <- fread(input = tagDist.path, header = T, check.names = F, data.table = F)
tagDist <- tagDist[, grep('BLANK', colnames(tagDist), invert = T)] # remove blank wells
sam <- fread(input = sam.path, sep = '\t', header = F, 
             stringsAsFactors = F, strip.white = T, fill = T, data.table = F)

# get tags info
sam <- sam[, c(1,3,4)] # retain only tag, chrom and position cols
colnames(sam)[1:3] <- c('Tag', 'chrom', 'pos') # rename cols
sam$Tag <- sub('tagSeq=', '', sam$Tag) # strip tags of extra letters
sam <- sam[! sam$chrom %in% c('chrUn', '*') & !duplicated(sam[, 2:3]), ] # remove tags if unanchored or ambiguous positions

cat("Total number of chromosomes:", length(unique(sam$chrom)), '\n', # count of unique chromosomes
    'Levels:', paste(sort(unique(sam$chrom)), collapse = ', ')) # print out list of chromosome identifiers

# checking duplicated tags in tag dist and sam file
sum(duplicated(tagDist$Tag)); sum(duplicated(sam$Tag))

# checking if all tags in tag dist match in tags distribution file
all(sum(sam$Tag %in% tagDist$Tag))

# tags position correction
if (need.chrSizeCorrection) {
  # fix positions
  for (i in 1:nrow(chrom.info)) {
    indices <- grep(chrom.info$part_name[i], sam$chrom, ignore.case = T)
    sam$pos[indices] <- sam$pos[indices] + chrom.info$addThis[i]
  }
  # fix chrom names
  sam$chrom <- sub('_part.', '', sam$chrom)
  sam$chrom <- sub('chr', '', sam$chrom) # remove chr from chrom names
  cat("Total number of chromosomes:", length(unique(sam$chrom))) # chrom number
} else {
  sam$chrom <- sub('chr', '', sam$chrom) # remove chr from chrom names
}

# print out list of chromosome identifiers
sort(unique(sam$chrom))

# merge tag distribution and position data.frame, and order by chrom and positions
tagDist <- merge(sam, tagDist, by = 'Tag', sort = F)
tagDist <- tagDist[order(tagDist$chrom, tagDist$pos), ]

table(tagDist$chrom) # get chrom count info

# check if size correction done properly (any tag position should be higher than the previous)
for (i in 1:nrow(chrom.info)) {
  all.true <- all(diff(tagDist$pos[tagDist$chrom == chrom.info$chrom_simple[i]]) > 0)
  if (all.true) cat('Chr', i, '=', chrom.info$chrom_simple[i], 'is all good! \n') else cat('Something wrong !!!')
}

# remove unnecessary objects
rm(sam)

#######################
## DELETION ANALYSIS ##
#######################

# create a data frame with desired bin intervals
# bin size in mb
bin.size.mb = 100
step.size = 50
# max chrom length
max.chrom.length = bin.size.mb * ceiling(max(chrom.info$length) / 10^6 / bin.size.mb)
# total number of bins
num.bins = max.chrom.length / step.size
# bins data.frame
chrom.bins <- data.frame('bin' = 1:num.bins, 'lower' = (1:num.bins - 1) * step.size * 10^6,
                            'upper' = (((1:num.bins - 1) * step.size) + bin.size.mb) * 10^6)

# create a data frame with chrom, position, and empty bin column
chr.pos.bin <- data.frame(tagDist[, 2:3], 'bin' = NA, 'chr.bin' = NA, stringsAsFactors = F)
# loop to fill up the empty bin column
for (i in 1:nrow(chrom.bins)) {
  # use the lower and upper bounds of a bin to find which rows 
  # from the above chrom, pos fall in the range and assign a bin
  positions <- which(chr.pos.bin$pos > chrom.bins$lower[i] & chr.pos.bin$pos <= chrom.bins$upper[i])
  chr.pos.bin$bin[positions] <- chrom.bins$bin[i]
}
# create unique chrom-bin combination
chr.pos.bin$chr.bin = paste0(chr.pos.bin$chrom, '_', chr.pos.bin$bin)

# get sum of tags per line per chromosome bin
tagDist <- as_tibble(tagDist)
tags.per.bin <- tagDist %>%
  group_by(chr.pos.bin$chr.bin) %>%
  summarise_at(.vars = colnames(tagDist)[-c(1:3)], .funs = sum)

# normalize number of tags per bin
# total number of tags per line
tot.tags.per.line <- colSums(tags.per.bin[, -1])
tot.tags.per.line <- data.frame('FullSampleName' = names(tot.tags.per.line), 
                                'tagCount' = tot.tags.per.line, stringsAsFactors = F)

# median of total tags per line
median.tag.num = median(tot.tags.per.line$tagCount)
# coefficients to adjust total tags per line
norm.coefficients <- median.tag.num / tot.tags.per.line$tagCount 
# multiply tags/line/bin with norm coefficient
normal.tags.per.bin <- cbind(tags.per.bin[, 1],
                                sweep(x = tags.per.bin[, -1],
                                      MARGIN = 2, STATS = norm.coefficients, FUN = '*'))

# check total tags again (should be all same)
range(colSums(normal.tags.per.bin[, -1]))

# get median normal bin tags across all lines
median.bin.across.lines <- apply(normal.tags.per.bin[, -1], 1, median, na.rm=T)
# normalize bin tags by median of bin tags across samples
normal.tags.per.bin[, -1] <- normal.tags.per.bin[, -1] / median.bin.across.lines

# median bin tags per line
median.tags.per.line <- apply(normal.tags.per.bin[, -1], 2, median)
# normalize bin tags by median of bin tags per line
normal.tags.per.bin[, -1] <- sweep(x = normal.tags.per.bin[, -1],
                                      MARGIN = 2, STATS = median.tags.per.line, FUN = '/')

# histogram of raw reads and unique tags per line
pdf(file = 'Fig.1_Histogram_counts.pdf', height = 4, width = 11)
par(mfrow = c(1,3), oma = c(0,0.5,0,0))
# raw read counts
read.counts <- read.table('readCounts.txt', header = T, as.is = T)
cat('Total number of raw reads =', sum(read.counts$goodBarcodedReads))
read.counts <- read.counts[grep('blank', read.counts$FullSampleName, ignore.case = T, invert = T), ]
hist(read.counts$goodBarcodedReads / 10^6, breaks = 50, cex.lab = 1.5, cex.axis = 1.25,
     xlab = expression(paste("Raw read count (10"^"6",')')), ylab = 'Number of samples', main = '', col = 'black')
abline(v = median(read.counts$goodBarcodedReads / 10^6), lty = 2, col = 'red')
legend('topright', legend = expression(bold('(A)')), pch = '', bg = 'gray', x.intersp = -0.5, cex = 1.25)
# unique tag counts
hist(tot.tags.per.line$tagCount / 10^3, breaks = 50, cex.lab = 1.5, cex.axis = 1.25,
     xlab = expression(paste("Unique tag count (10"^"3",')')), ylab = "Number if samples", main = '', col = 'black')
abline(v = median(tot.tags.per.line$tagCount / 10^3), lty = 2, col = 'red')
legend('topright', legend = expression(bold('(B)')), pch = '', bg = 'gray', x.intersp = -0.5, cex = 1.25)
# distribution of relative normalized tag counts
hist(as.numeric(unlist(normal.tags.per.bin[, -1])), breaks = 60, cex.lab = 1.5, cex.axis = 1.25,
     xlab = 'Normalized tag count per bin', ylab = 'Number of bins', main = '', col = 'black')
legend('topright', legend = expression(bold('(C)')), pch = '', bg = 'gray', x.intersp = -0.5, cex = 1.25)
dev.off()

####
## explore the reason for skewed distribution
###

sample.meta <- read.csv('sampleMetadata.csv', header = T, as.is = T, na.strings = '-')
sample.meta <- sample.meta[order(sample.meta$Flowcell_order), ]
sample.meta <- sample.meta[grep('blank', sample.meta$FullSampleName, ignore.case = T, invert = T), ]

sample.meta <- merge(tot.tags.per.line, sample.meta, by = 'FullSampleName')
sample.meta <- merge(read.counts[, -3], sample.meta, by = 'FullSampleName')

table(sample.meta$Flowcell, sample.meta$plate_id)

# plot various distributions to find out which variable is causing the skewness
pdf(file = 'Fig.2_Boxplots_flowcell_plates.pdf', height = 8.5, width = 6)
par(mfrow = c(2,2), oma = c(0.5,0.5,0.5,0.5))

boxplot(sample.meta$goodBarcodedReads / 10^6 ~ sample.meta$Flowcell_order, main = 'Flowcells', cex.main = 1.5,
        cex.axis = 1.2, las = 2, xaxt = 'n', yaxt = 'n', col = alpha(c('blue', 'orange', 'maroon'), 0.5))
axis(side = 4, at = 0:50, labels = F)
axis(side = 1, at = 1:length(unique(sample.meta$Flowcell)), labels = F)
axis(side = 1, at = 1:length(unique(sample.meta$Flowcell)), labels = unique(sample.meta$Flowcell), 
     tick = F, cex.axis = 1.25, las = 2, line = 0.5)
mtext(text = expression(paste("Raw read count (10"^"6",')')), side = 2, cex = 1.25, line = 1.25)
legend('top', legend = expression(bold('(A)')), pch = '', bg = 'gray', x.intersp = -0.5, cex = 1)

boxplot(sample.meta$goodBarcodedReads / 10^6 ~ sample.meta$plate_id, main = 'DNA Plates', cex.main = 1.5,
        cex.axis = 1.2, las = 2, yaxt = 'n', col = alpha(c('blue', 'blue', 'blue', 'orange', 'orange', 'orange', 'maroon'), 0.5))
axis(side = 2, at = 0:50, hadj = 4, cex.axis = 1.25, las = 1)
legend('top', legend = expression(bold('(B)')), pch = '', bg = 'gray', x.intersp = -0.5, cex = 1)

boxplot(sample.meta$tagCount / 10^3 ~ sample.meta$Flowcell_order,
        xaxt = 'n', yaxt = 'n', col = alpha(c('blue', 'orange', 'maroon'), 0.5))
axis(side = 3, at = 1:length(unique(sample.meta$Flowcell)), labels = F)
axis(side = 4, at = seq(0, 2000, 200), labels = F)
mtext(text = expression(paste("Unique tag count (10"^"3",')')), side = 2, cex = 1.25, line = 1.25)
legend('top', legend = expression(bold('(C)')), pch = '', bg = 'gray', x.intersp = -0.5, cex = 1)

boxplot(sample.meta$tagCount / 10^3 ~ sample.meta$plate_id, xaxt = 'n', yaxt = 'n',
        col = alpha(c('blue', 'blue', 'blue', 'orange', 'orange', 'orange', 'maroon'), 0.5))
axis(side = 3, at = 1:length(unique(sample.meta$plate_id)), labels = F, cex.axis = 1.2)
axis(side = 2, at = seq(0, 2000, 200), labels = F)
axis(side = 2, at = seq(0, 2000, 200), tick = F, hadj = 0.5, cex.axis = 1.25, las = 1, line = 2)
legend('top', legend = expression(bold('(D)')), pch = '', bg = 'gray', x.intersp = -0.5, cex = 1)

dev.off()

# merge data frames to get chrom position, bin num and other info for plotting
# keep only unique bins
chr.pos.bin.nodup <- chr.pos.bin[!duplicated(chr.pos.bin$chr.bin), ]
# merge to get chrom and bin cols
normal.tags.per.bin <- merge(chr.pos.bin.nodup, normal.tags.per.bin, by.x = 4, by.y = 1, all.y = T)
# merge to get barplot values
normal.tags.per.bin <- merge(chrom.info[, c(2,8)], normal.tags.per.bin, by.x = 1, by.y = 2)
# merge to get lower and upper bounds
normal.tags.per.bin <- merge(chrom.bins, normal.tags.per.bin, by.x = 1, by.y = 5)

# order based on chromosome and bins
normal.tags.per.bin <- normal.tags.per.bin[
  order(normal.tags.per.bin$chrom_simple, 
        normal.tags.per.bin$bin), ]

# replace bins' upper bounds with original length of the chromosome if greater than the length
max.bin.upper <- aggregate(upper ~ chrom_simple, normal.tags.per.bin, max)
max.bin.upper <- merge(max.bin.upper, chrom.info[, c(2,5)], by = 1)
for (i in 1:nrow(max.bin.upper)) {
  index <- which(normal.tags.per.bin$chrom_simple == max.bin.upper$chrom_simple[i] &
                   normal.tags.per.bin$upper > max.bin.upper$length[i])
  normal.tags.per.bin$upper[index] <- max.bin.upper$length[i]
}
# check if it worked
max.bin.upper <- cbind(max.bin.upper, 
                          'last.win.upper' = aggregate(upper ~ chrom_simple, normal.tags.per.bin, max)[, 2])

# convert bp to Mb
normal.tags.per.bin$lower <- normal.tags.per.bin$lower / 10^6
normal.tags.per.bin$upper <- normal.tags.per.bin$upper / 10^6


# lines info
lines <- data.frame('TA' = colnames(normal.tags.per.bin)[-c(1:7)], stringsAsFactors = F)
lines <- merge(lines, sample.meta, by = 1, sort = F)
lines <- lines[order(lines$notes), ]

# paint karyotypes all
pdf(paste0(pdf.outfile,'.pdf'), width = 11, height = 8.5)
# plot multiple samples per sheet
par(mfrow=c(3,3))
# max chrom length to plot graphically and some padding
chrom.y.lim.len <- (100 * ceiling(max(chrom.info$length / 10^6 / 100))) + 80
# spread for chrom tags
spread = 0.45
for (i in 1:nrow(lines)) {
  plot.karyotype(lines = lines)
}
dev.off()

# lines for reference visualization
ref.lines <- read.table('ref_lines.txt', header = T, as.is = T)

# paint karyotypes for ref lines
pdf(paste0('Fig.3_Karyotypes.pdf'), width = 11, height = 6)
# plot multiple samples per sheet
par(mfrow=c(2,3))
# max chrom length to plot graphically and some padding
chrom.y.lim.len <- (100 * ceiling(max(chrom.info$length / 10^6 / 100))) + 80
# spread for chrom tags
spread = 0.45
for (i in 1:nrow(ref.lines)) {
  plot.karyotype(lines = ref.lines)
}
dev.off()

#########
## END ##
#########

