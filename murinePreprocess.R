###########################################################################
## Publication:
## Identification of Potential Biomarkers of Vaccine Reactogenicity
## Authors:
## Paul F. McKay, Deniz Cizmeci, Yoann Aldon, Jeroen Maertzdorf,
## January Weiner 3rd, Stefan H. E. Kaufmann, David J. M. Lewis,
## Robert A. van den Berg, Giuseppe Del Giudice and Robin J. Shattock

## Preprocessing of raw trascriptomics data
## Deniz Cizmeci, Imperial College London - 2017

## Raw data is available:
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120661
###########################################################################


# Load Libraries ----------------------------------------------------------
library(limma)
library("ggplot2")
library(ggthemes)


# Import raw data ---------------------------------------------------------
targets <- readTargets("targets_Imperial_mouse.txt")
dat.raw <- read.maimages(targets, source="agilent", green.only=T)
dat.bc <- backgroundCorrect(dat.raw, method="normexp")
dat <- normalizeBetweenArrays(dat.bc)


# Filter probes -----------------------------------------------------------
dat2 <- dat[ dat$genes$ControlType == 0, ] #remove control probes
dat2 <- avereps(dat2, ID=dat2$genes$ProbeName) # Condense a microarray 
# data object so that values for within-array replicate probes are replaced
# with their average.
mdat <- dat2[ ! grepl( "^linc", dat2$genes$Description ), ] #remove long 
# non-coding RNA probes


# Rename columns ----------------------------------------------------------
mdat$targets$code <- row.names(mdat$targets)
identical(row.names(mdat$targets), colnames(mdat$E))
mdat$targets$treat.tp <- paste( mdat$targets$treatment, sprintf("%02d", mdat$targets$timepoint), sep=".")
mdat$targets$id <- paste( "m", mdat$targets$tissue, mdat$targets$treatment, sprintf("%02d", mdat$targets$timepoint), mdat$targets$animal, sep=".")
colnames(mdat) <- mdat$targets$id
mice <- mdat

mice$targets$tp <- as.factor(mice$targets$timepoint)

# Colour codes ------------------------------------------------------------
# Colours used for different factors in this study
muscleColour <- "#BEAED4"
lnColour <- "#7FC97F"
blColour <- "#FDC086"

EngerixColour <- "#99E7FF"
IFAColour <- "#808080"
LPSColour <- "#00ABFF"
PertussisColour <- "#FFA500"
PolyICColour <- "#4CA64C"
TriFluColour <- "#9B9BF1"
TriFluMF59Colour <- "#CC0000"
SalineColour <- "#000000"


# PCA ---------------------------------------------------------------------
#Perform Principal Component Analysis (PCA)
pca.all <- prcomp( t(mice$E), scale.=T)
#Percentage of variance explained
PoV <- pca.all$sdev^2/sum(pca.all$sdev^2)

mice.pca <- as.data.frame(pca.all$x)
identical(row.names(mice.pca), row.names(mice$targets))

pdf("murinePCA_all.pdf")
ggplot(mice.pca, aes(PC1, PC2, color=mice$targets$tissue)) + geom_point() +
  scale_colour_manual(values = c(blColour, lnColour, muscleColour),
                      name = "Tissue",
                      breaks = c("blood", "MLN", "muscle"),
                      labels = c("Blood", "Lymph Nodes", "Muscle")) +
  xlab(paste0("Principal Component 1 ", "(%",sprintf("%1.2f",PoV[1]*100),")")) +
  ylab(paste0("Principal Component 2 ", "(%",sprintf("%1.2f",PoV[2]*100),")")) +
  geom_rangeframe() + theme_tufte() + theme(
    axis.text.x = element_blank(),axis.text.y = element_blank(),
    axis.ticks = element_blank(),axis.line = element_line(colour = "gray"))
dev.off()


# Separate objects per tissue type ----------------------------------------
mice.bl <- mice[ , mice$targets$tissue == "blood" ]
mice.ln <- mice[ , mice$targets$tissue == "MLN" ]
mice.mu <- mice[ , mice$targets$tissue == "muscle" ]

#Perform PCA
mypca <- function(df){
  pca.df <- prcomp( t(df$E), scale.=T)
  PoV <- pca.df$sdev^2/sum(pca.df$sdev^2)
  
  pca <- as.data.frame(pca.df$x)
  
  cols <- c("Pertussis" = PertussisColour,"TriFlu" = TriFluColour,
            "TriFluMF59" = TriFluMF59Colour, "saline" = SalineColour,
            "Engerix" = EngerixColour,"LPS" = LPSColour,
            "PolyIC" = PolyICColour,"IFA" = IFAColour)
  
  labs <- c("Pertussis" = "Pentavac SD","TriFlu" = "Agrippal",
            "TriFluMF59" = "Fluad", "saline" = "Saline",
            "Engerix" = "Engerix B","LPS" = "LPS",
            "PolyIC" = "Poly I:C","IFA" = "IFA")
  
  g <- ggplot(pca, aes(PC1, PC2, color=df$targets$treatment, shape = df$targets$tp)) + geom_point() +
    scale_colour_manual(values = cols,
                        name = "Immunisation",
                        labels = labs) +
    scale_shape_manual(values=c(0,1,2,3,4,5,6), name = "Timepoints") +
    xlab(paste0("Principal Component 1 ", "(%",sprintf("%1.2f",PoV[1]*100),")")) +
    ylab(paste0("Principal Component 2 ", "(%",sprintf("%1.2f",PoV[2]*100),")")) +
    geom_rangeframe() + theme_tufte() + theme(
      axis.text.x = element_blank(),axis.text.y = element_blank(),
      axis.ticks = element_blank(),axis.line = element_line(colour = "gray"))

  return(g)
}

pdf("murinePCA_blood.pdf")
mypca(mice.bl)
dev.off()

pdf("murinePCA_ln.pdf")
mypca(mice.ln)
dev.off()

pdf("murinePCA_muscle.pdf")
mypca(mice.mu)
dev.off()

