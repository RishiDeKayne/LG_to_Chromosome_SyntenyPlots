#script by Rishi De-Kayne 2019 (Eawag/UniBern) for plotting Linkage map vs. chromosome synteny
#https://github.com/RishiDeKayne/LG_to_Chromosome_SyntenyPlots

#########################################################################
#three input files are needed: 
# 1. FalconChromosomeLenghts.txt - a list of chromosome lengths 
#       - in my case this included the lengths of the fasta headings (all < 50) followed by the # bases in the chromosome fasta line
# 1. e.g. 
#45
#111295376
#45
#110734514
#45
#108247109
#44
#92930148

# 2. LM_stats.txt - a summary statstics table of the linkage (taken from Table 1: http://www.g3journal.org/content/8/12/3745.figures-only)
# 2. e.g. 
#Calb01 253 75.96 0.30 91.07 0.36 63.67 0.25 Ssa01 W02 1.43
#Calb02 228 83.57 0.37 101.33 0.44 69.58 0.31 Ssa01 W03 1.46
#Calb03 220 78.51 0.36 84.40 0.38 87.95 0.40 Ssa21 W32 0.96

# 3. SA_linkagemap_FalconPhaseMapped_Filtered.csv - this is the result of mapping all linkage map rad loci to the reference and filtering for MAPQ > 30 using: 
#       - bsub -n10 -W 4:00 -R "rusage[mem=3000]" "module load java gdc bwa/0.7.17; bwa mem -t 10 module load java gdc bwa/0.7.17; bwa mem -t 10 chromosome_level_assembly.fasta SpeciesAveragedLinkageMap.fasta > SA_linkagemap_FalconPhaseMapped.sam"
#       - awk '{ if($5 > 30) {print}}' SA_linkagemap_FalconPhaseMapped.sam > SA_linkagemap_FalconPhaseMapped_Filtered.txt

#       - this was then checked over and the column with mapping location separated by scaffold number and bp location e.g. - check columns in test file carefully
#@PG			ID:bwa	PN:bwa		VN:0.7.17-r1188	CL:bwa	mem	#NAME?	10	/cluster/project/gdc/shared/p298/Rishi/PostPhaseScaffoldingGenomes/WF_Canu.chr.fasta	/cluster/project/gdc/shared/p298/Rishi/LinkageMap/SpeciesAveragedMap.fasta						
#consensus_42108	SA10	0	0	PGA	scaffold21	16735146	60	90M	*	0	0	CACTCTCTCTCCTTCTCCCTCCTTTCCTCTCTATGTCCCAGAAAGCCTTCAAAGATGACCCGCAGAGCAGAGCGGTGGTATGATCCTGCA	*	NM:i:0	MD:Z:90	AS:i:90	XS:i:21	
#consensus_27024	SA10	2.3	0	PGA	scaffold21	6158036	60	90M	*	0	0	TGCAGGGTTCAAATGCAACCACTTGGGCAGCCTGGTCTCATAGGCTTGACGTAACATAGAAAATGTAAATCTGGGACACTCAAATTAGTA	*	NM:i:0	MD:Z:90	AS:i:90	XS:i:23	
#consensus_7502	SA10	2.63	16	PGA	scaffold21	26075329	60	90M	*	0	0	GCAGCTGCTGCTTGTGGGCCTCCTTAATTGCTTGCAGCTGGCGTTGGAGCAGCAGCTCTGATCTCAGGCCATCCTCCAGGAAAGCCTGCA	*	NM:i:1	MD:Z:25C64	AS:i:85	XS:i:20	

#for extra information about plotting check out the circlize manual: https://jokergoo.github.io/circlize_book/book/ 

#########################################################################

#-----------------------------------------------------------------------
#the first part of this script will not change the scaffold order meaning the circos plot is hard to read
#it should be run to check everything is working - the synteny pattern should be obvious but will look better once ordering is carried out later

#load the circlize library
library(circlize)

#load in chromosome lengths and linkage map lengths
setwd("~/GitHub/LG_to_Chromosome_SyntenyPlots/")
#make df with chromosome lengths
lengthsdf <- read.csv(file = "FalconChromosomeLengths.txt", header = FALSE)
#remove header values so we are left in our case with 40 lines with the number of bp
lengthsdf <- subset(lengthsdf, lengthsdf$V1 > 50)
#add scaffold number to length df - these start from 0 since our scaffolds are numbered that way 
lengthsdf$Scaffold <- 0:39

#now read in linkage map statistics
lm <- read.csv(file = "LM_stats.csv", header = FALSE, sep = ",")

#load in samfile already separated columns in excel - DOUBLE CHECK THE FORMAT HERE!
samfile <- read.csv(file = "SA_linkagemap_FalconPhaseMapped_Filtered.csv")

#remove any mappings which map to contigs rather than chromosome level scaffolds
samfile <- subset(samfile, as.character(samfile$PN.bwa) == "PGA")

#make synteny map structure file with each linkage group/chromosome and start and end
#find equivalent chromosomes - in our case we do not know which chromosome goes with which linkage group
#this loop will go through each linkage group in our lm df and search the sam file for mappings corresponding to that lg
#then it will look through the chromosomes those markers are mapping to and identify the most abundant one
#columns are added to the lm data frame to show the proportion of markers mapping to this new 'most abundant' syntenic chromosome 
for (i in 1:40){
  linkagegroup <- levels(lm$V1)[i]
  new_linkagegroup <- gsub("Calb", "SA", linkagegroup)
  sam_subset <- subset(samfile, as.character(samfile$X) == as.character(new_linkagegroup))
  scaffoldvect <- as.vector(sam_subset$X.2)
  abundant <- names(sort(summary(as.factor(scaffoldvect)), decreasing=T)[1:1])
  lm$scaffold[i] <- abundant
  lm$tot_markernumber[i] <- length(sam_subset$X.PG)
  abundant_number <- subset(sam_subset, as.character(sam_subset$X.2) == as.character(abundant))
  lm$scaf_markernumber[i] <- length(abundant_number$X.PG)
  lm$proportion_abundant <- ((lm$scaf_markernumber/lm$tot_markernumber)*100)
}

#order the lm dataset
orderedbyscaffolds <- lm[order(lm$scaffold), ]

#for the circos plot two files are key - the first is to specify the structure of the plot
#this describes the number of segments to plot - in my case 40 LGs and 40 chromosomes and the start and end of each
#circos can work out how to plot these so long as each segment has a unique name and two values, the start and end point
# e.g. SA1 0
#      SA1 100
#      SA2 0
#      SA2 200
#will plot two segments with the first half the size of the second

#we prepare this specifications df separately for LGs and chromosomes
#for linkage groups (lg_dat)
lg_dat_1 <- as.data.frame(lm$V1)
lg_dat_1$cM <- lm$V3

lg_dat_2 <- as.data.frame(lm$V1) 
lg_dat_2$cM <- 0

lg_dat <- rbind(lg_dat_2, lg_dat_1)

lg_dat <- lg_dat[order(lg_dat$`lm$V1`),]


#and for the chromosomes - this is a bit more involved since some chromosomes may not have been the 'most abundant' for any of the linkage groups
#now get scaffold order for plotting:
chrom_order_df <- unique(lm$scaffold) 
chrom_order_df

###••••••••#####
#check lm$scaffold manually - some scaffolds will be missing! you need to add these manually below
#they are likely not that key since they will have few markers mapping to them but should be included
missing <- c("scaffold6", "scaffold7","scaffold34", "scaffold37", "scaffold38", "scaffold39")
#now paste these together to get the full chromosome list
fullchroms <- c(chrom_order_df, missing)
chromosomes <- as.data.frame(fullchroms)

#now we use the lengths dataframe from earlier to assign the lengths of each chromosome 
#(these are currently in bp - we need to convert them to cM later so the lengths of each segment are comparable)
for (i in 1:40){
  searchscaffold <- as.character(chromosomes$fullchroms[i])
  new_searchscaffold <- gsub("scaffold", "", searchscaffold)
  lengthsubset <- subset(lengthsdf, as.character(lengthsdf$Scaffold) == as.character(new_searchscaffold))
  chromosomes$bp[i] <- lengthsubset$V1
}

#chromosomes are renamed with a letter 'w' to make plotting easier
#circos will always plot in a clockwise manner and in an alphabetic way with segment names - in my case I want LG on the right and chromosome on the left
#lGs can therefore be named 'c' and will be plotted first followed by 'W' segments
chromosomes$newchrom <- 1:40
chromosomes$newchrom <- paste("W", as.character(chromosomes$newchrom), sep = "")
for (i in 1:length(chromosomes$newchrom)){
  if (nchar(as.character(chromosomes$newchrom))[i] == 2){
    chromosomes$newchrom[i] <- as.character(gsub("W", "W0", as.character(chromosomes$newchrom)[i]))
  }
}

chromosomes$cm <- chromosomes$bp
#make new df just with the new chromosome names
renamedchromosomes <- chromosomes[,3:4]

#to calclulate conversion work out the max length of the chromosomes
totalbp <- sum(as.numeric(chromosomes$bp))
#and the max length of linkage groups 
totalcm <- sum(lg_dat$cM)

#this gives us total cM and total bp 
#work out bp in cm - we make a conversion factor to make the distances on each segment cmoparable
conversion <- totalcm/totalbp

#and prepare the chromosome segment data as we did for the linakge groups
chrom_dat_1 <- renamedchromosomes
chrom_dat_1$cm <- as.numeric(chromosomes$bp)*conversion

chrom_dat_2 <- chrom_dat_1
chrom_dat_2$cm <- 0

chrom_dat <- rbind(chrom_dat_2, chrom_dat_1)

chrom_dat <- chrom_dat[order(chrom_dat$newchrom),]

#now merge these two files to get a complete segment info df
#first make column the same
colnames(lg_dat) <- c("segment", "cM")
lg_dat$segment <- gsub("Calb", "c", lg_dat$segment)

colnames(chrom_dat) <- c("segment", "cM")
full_dat <- rbind (lg_dat, chrom_dat)

#full_dat should now have twice as many rows as you want segments - one for the start value and one for the end value 

#now try circos plot, track outlines are given using 'complete' - dots can be added from this - then want to add links
library(circlize)
#this just makes sure the graphical parameters are set correctly
circos.clear()

#par sets the graphical paremeters
circos.par("track.height" = 0.05, start.degree=90, cell.padding = c(0.02, 0, 0.02, 0))
#initialize the track itself specifying the full_dat file we created
circos.initialize(factors = full_dat$segment, x = full_dat$cM)
#now each segment should be plotted - it will look messy but dont worry
circos.track(factors = full_dat$segment, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
#you may get this warning: Note: 1 point is out of plotting region in sector 'c01', track '1'.
#dont worry about it!

#now we need to rename LG and scaffolds in samfile to match the segments we have
#this way circos knows where it should plot everything
#take beginning of our custom sam file
linkfile <- samfile[,1:7]
#change all chromosome bp values to cM for plot
linkfile$VN.0.7.17.r1188 <- linkfile$VN.0.7.17.r1188*conversion

#find chromosomes conversion in df: chromosomes
#renaming chromosomes in linkfile based on chromosomes df
for (i in 1:length(linkfile$X.PG)){
  oldchrom <- as.character(linkfile$X.2[i])
  oldlg <- as.character(linkfile$X)[i]
  newchrom_name_df <- subset(chromosomes, as.character(chromosomes$fullchroms) == as.character(oldchrom))
  linkfile$newchrom[i] <- newchrom_name_df$newchrom
  linkfile$newLG[i] <- gsub("SA", "c", oldlg)
}

#clear the old plot and try again using new chromosomes names - this time we will plot using the linkfile
circos.clear()
#then set circos plot parameters
circos.par("track.height" = 0.05, start.degree=90, cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = full_dat$segment, x = full_dat$cM)
circos.track(factors = full_dat$segment, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
#loop through linkage data frame and extract info to make links 
#here we specify the LG and cM position on that linkage group a and b (where the link should start)
#and the chromosome name and converted 'cM' position the link should go to
for (row in 1:(nrow(linkfile))){
  a <- linkfile$newLG[row]
  b <- linkfile$X.1[row]
  e <- linkfile$newchrom[row]
  d <- linkfile$VN.0.7.17.r1188[row]
  circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.5, col = "black")
}
#here you should see the links are not totally random but generally going from one linkage group to one chromosomes
#since they have not been reordered the first segment i.e. at 1 o'clock (going clockwise) should correspond to the first segment from 6 o'clock (going clockwise)


#-----------------------------------------------------------------------
#this second part of this script will change the scaffold order so the circos plot is easier to read

#load in chromosome lengths and linkage map lengths
#make df with chromosome lengths
lengthsdf <- read.csv(file = "FalconChromosomeLengths.txt", header = FALSE)
#remove header values so we are left in our case with 40 lines with the number of bp
lengthsdf <- subset(lengthsdf, lengthsdf$V1 > 50)
#add scaffold number to length df - these start from 0 since our scaffolds are numbered that way 
lengthsdf$Scaffold <- 0:39

#now read in linkage map statistics
lm <- read.csv(file = "LM_stats.csv", header = FALSE, sep = ",")

#load in samfile already separated columns in excel - DOUBLE CHECK THE FORMAT HERE!
samfile <- read.csv(file = "SA_linkagemap_FalconPhaseMapped_Filtered.csv")
samfile <- subset(samfile, as.character(samfile$PN.bwa) == "PGA")

#make synteny map structure file with each linkage group/chromosome and start and end
#find equivalent chromosomes
for (i in 1:40){
  linkagegroup <- levels(lm$V1)[i]
  new_linkagegroup <- gsub("Calb", "SA", linkagegroup)
  sam_subset <- subset(samfile, as.character(samfile$X) == as.character(new_linkagegroup))
  scaffoldvect <- as.vector(sam_subset$X.2)
  abundant <- names(sort(summary(as.factor(scaffoldvect)), decreasing=T)[1:1])
  lm$scaffold[i] <- abundant
  lm$tot_markernumber[i] <- length(sam_subset$X.PG)
  abundant_number <- subset(sam_subset, as.character(sam_subset$X.2) == as.character(abundant))
  lm$scaf_markernumber[i] <- length(abundant_number$X.PG)
  lm$proportion_abundant <- ((lm$scaf_markernumber/lm$tot_markernumber)*100)
}

orderedbyscaffolds <- lm[order(lm$scaffold), ]

#now get lg data for plotting
lg_dat_1 <- as.data.frame(lm$V1)
lg_dat_1$cM <- lm$V3

lg_dat_2 <- as.data.frame(lm$V1) 
lg_dat_2$cM <- 0

lg_dat <- rbind(lg_dat_2, lg_dat_1)

lg_dat <- lg_dat[order(lg_dat$`lm$V1`),]

#now get scaffold order for plotting:
chrom_order_df <- unique(lm$scaffold) 
chrom_order_df
missing <- c("scaffold6", "scaffold7","scaffold34", "scaffold37", "scaffold38", "scaffold39")
fullchroms <- c(chrom_order_df, missing)
chromosomes <- as.data.frame(fullchroms)

for (i in 1:40){
  searchscaffold <- as.character(chromosomes$fullchroms[i])
  new_searchscaffold <- gsub("scaffold", "", searchscaffold)
  lengthsubset <- subset(lengthsdf, as.character(lengthsdf$Scaffold) == as.character(new_searchscaffold))
  chromosomes$bp[i] <- lengthsubset$V1
}

#newrename chromosomes
#this time Z is used to denote a re-ordered chromosome
#this order is now reversed here (previously was 1:40)
#therefore the most abundant chromosome for linkage group 1/40 should now be numbered 40/40 (rather thann 1/40)
#this way it will be plotted at the top of the circos plot not at the bottom i.e.
# 40/40 chromosome ---- 1/40 LG
# 39/40 chromosome ---- 2/40 LG
# 38/40 chromosome ---- 3/40 LG
# 37/40 chromosome ---- 4/40 LG

chromosomes$newchrom <- 40:1
chromosomes$newchrom <- paste("Z", as.character(chromosomes$newchrom), sep = "")
for (i in 1:length(chromosomes$newchrom)){
  if (nchar(as.character(chromosomes$newchrom))[i] == 2){
    chromosomes$newchrom[i] <- as.character(gsub("Z", "Z0", as.character(chromosomes$newchrom)[i]))
  }
}

chromosomes$cm <- chromosomes$bp
renamedchromosomes <- chromosomes[,3:4]

totalbp <- sum(as.numeric(chromosomes$bp))
totalcm <- sum(lg_dat$cM)


#work out bp in cm
conversion <- totalcm/totalbp

chrom_dat_1 <- renamedchromosomes
chrom_dat_1$cm <- as.numeric(chromosomes$bp)*conversion

chrom_dat_2 <- chrom_dat_1
chrom_dat_2$cm <- 0

chrom_dat <- rbind(chrom_dat_2, chrom_dat_1)

chrom_dat <- chrom_dat[order(chrom_dat$newchrom),]

#now merge these two files:
#first make column the same
colnames(lg_dat) <- c("segment", "cM")
lg_dat$segment <- gsub("Calb", "c", lg_dat$segment)


colnames(chrom_dat) <- c("segment", "cM")
#this time we want values to be negative so plotting is carried out correctly
chrom_dat$cM <- (chrom_dat$cM)-(2*(chrom_dat$cM))


full_dat <- rbind (lg_dat, chrom_dat)
#full dat again speficies the plotting parameters

#test plot again to make sure segments are plotted correctly
circos.clear()
#then set circos plot parameters
circos.par("track.height" = 0.05, start.degree=90, cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = full_dat$segment, x = full_dat$cM)
circos.track(factors = full_dat$segment, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })

#now rename LG and scaffolds in samfile
#check that all the information is included at this step
linkfile <- samfile[,1:8]
#change bp to cM for plot
linkfile$VN.0.7.17.r1188 <- linkfile$VN.0.7.17.r1188*conversion

#now make chromosome bp/cM positions negative in linkfile
linkfile$VN.0.7.17.r1188 <- linkfile$VN.0.7.17.r1188-(2*as.numeric(linkfile$VN.0.7.17.r1188))


#find chromosomes conversion in df: chromosomes
for (i in 1:length(linkfile$X.PG)){
  oldchrom <- as.character(linkfile$X.2[i])
  oldlg <- as.character(linkfile$X)[i]
  newchrom_name_df <- subset(chromosomes, as.character(chromosomes$fullchroms) == as.character(oldchrom))
  linkfile$newchrom[i] <- newchrom_name_df$newchrom
  linkfile$newLG[i] <- gsub("SA", "c", oldlg)
}

#and now plot the reordered data
circos.clear()
#then set circos plot parameters
circos.par("track.height" = 0.05, start.degree=90, cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = full_dat$segment, x = full_dat$cM)
circos.track(factors = full_dat$segment, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
#loop through linkage data frame and extract info to make links as before
for (row in 1:(nrow(linkfile))){
  a <- linkfile$newLG[row]
  b <- linkfile$X.1[row]
  e <- linkfile$newchrom[row]
  d <- linkfile$VN.0.7.17.r1188[row]
  #circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.9, col = "black")
  circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h=0.7, col = add_transparency("darkgrey", transparency = 0.35))
  #this time I also add lines in the segments to see where the links are actually being plotted to
  circos.lines(c(as.integer(b)), 1, as.character(a), col = "darkgrey", type = 'h')
  circos.lines(c(as.integer(d)), 1, as.character(e), col = "darkgrey", type = 'h')
}

#because of the plotting parameters it may still be hard to see synteny but once we colour the links it will be more obvious

#-----------------------------------------------------------------------
#now we plot the reordered data properly with colours to make the synteny stand out
library(circlize)
circos.clear()
#name of file/parameters
#if you want to save as a .tiff uncomment out this line and the dev.off() at the end of the plot - try first to make sure plotting is carried out correctly
#tiff("yourplot.tiff",height=8,width=8,units="in",res=300,compression="lzw")
# Customize plot parameters
par(fig=c(0,1,0,1),mar=c(1,1,1,1))
#the next gapping paramters determine how the segments are spaced - in my case I want a bigger gap at 12 and 6 o'clock between the LGs and chromosomes
#create gaps vector
gapsnorm <- c(1)
gapswide <- c(1)
wfgaps <- rep(gapsnorm, 39)
salmgaps <- rep(gapswide, 39)
startgap <- c(9)
endgap <- c(9)
gapsall <- c(wfgaps, startgap, salmgaps, endgap)

#create names vector - this we will use to name our segments -remember plotting is clockwise! hence the backwards naming of chromosomes
Wnamevector <- c("LG1", "LG02", "LG03", "LG04", "LG05", "LG06", "LG07", "LG08", "LG09", "LG10", 
                 "LG11", "LG12", "LG13", "LG14", "LG15", "LG16", "LG17", "LG18", "LG19", "LG20", 
                 "LG21", "LG22", "LG23", "LG24", "LG25", "LG26", "LG27", "LG28", "LG29", "LG30", 
                 "LG31", "LG32", "LG33", "LG34", "LG35", "LG36", "LG37", "LG38", "LG39", "LG40")
Snamevector <- c("Chr40", "Chr39", "Chr38", "Chr37", "Chr36", "Chr35", "Chr34", "Chr33", "Chr32",
                 "Chr31", "Chr30", "Chr29", "Chr28", "Chr27", "Chr26", "Chr25", "Chr24", "Chr23", "Chr22", "Chr21", "Chr20", 
                 "Chr19", "Chr18", "Chr17", "Chr16", "Chr15", "Chr14", "Chr13", "Chr12", "Chr11", "Chr10", 
                 "Chr09", "Chr08", "Chr07", "Chr06", "Chr05", "Chr04", "Chr03", "Chr02", "Chr01")

namevector <- c(Wnamevector, Snamevector)

#then set circos plot parameters including the rotation, gaps between segments and gaps between tracks - this will probably all need tweaking for each plot
circos.par("track.height" = 0.08 , start.degree=87, cell.padding = c(0.005, 0, 0.005, 0), gap.degree = gapsall, track.margin = c(0.0045, 0.0045))
#initialize
circos.initialize(factors = full_dat$segment, x = full_dat$cM)

Wlet <- "W"
#colour vector
vec <- c(add_transparency("red", transparency = 0.8), 
         add_transparency("orange", transparency = 0.8), 
         add_transparency("green", transparency = 0.8), 
         add_transparency("blue", transparency = 0.8), 
         add_transparency("purple", transparency = 0.8))
colvec <- rep(vec, 8)
#initialize first track - chromosomes/arms
circos.track(factors = full_dat$segment, ylim = c(0,1),
             panel.fun = function(x, y) {
               if(grepl(Wlet, namevector[CELL_META$sector.numeric.index])){
                 circos.text(CELL_META$xcenter, 0.9 + uy(2, "mm"), 
                             namevector[CELL_META$sector.numeric.index], cex = 0.3)
               } else {
                 circos.text(CELL_META$xcenter, 1.45 + uy(2, "mm"), 
                             namevector[CELL_META$sector.numeric.index], cex = 0.3)
               }
             })
#here it should already look much better with spacing between the segments etc.

#loop through linkage data frame and extract info to make links
for (row in 1:(nrow(linkfile))){
  a <- linkfile$newLG[row]
  b <- linkfile$X.1[row]
  e <- linkfile$newchrom[row]
  d <- linkfile$VN.0.7.17.r1188[row]
  fact <- as.factor(linkfile$newLG)[row]
  #circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.9, col = "black")
  #circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = add_transparency("darkgrey", transparency = 0.35))

  #plot the links - they should 
  circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = colvec[fact])
  
  
  #this time I wanted to colour the lines in the segments based on the mapping quality score in the linkfile column linkfile$CL.bwa
  if ((linkfile$CL.bwa)[row] > 59){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "royalblue1", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 60){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "lightblue", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 50){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "orchid1", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 40){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "red", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] > 59){
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "royalblue1", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 60){
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "lightblue", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 50){
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "orchid1", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 40){
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "red", type = 'h')
  }
  # if ((linkfile$CL.bwa)[row] < 60){
  #   circos.lines(c(as.integer(d)), 1, as.character(e), col = "red", type = 'h')
  # }
  # if ((linkfile$CL.bwa)[row] > 59){
  #   circos.lines(c(as.integer(d)), 1, as.character(e), col = "lightblue", type = 'h')
  # }
}

#and here I add a legend to explain the colouring of lines in segments - the number of links was counted separately by subsetting the df and counting rows
legend("topright", c("30-39 - n=110", "40-49 - n=177", "50-59 - n=188", "60 - n=3173"), col = c("red", "orchid1", "lightblue", "royalblue1"), lwd = 1.5, cex = 0.55, title = "MapQ of 3648 links")

#dev.off()

#here the last few chromosomes are those which we placed randomly so they may not fit ideally 
#the rest should be ordered by their syntney
