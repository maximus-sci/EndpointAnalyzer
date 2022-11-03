library(readr)
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stats)
library(scales)
library(ggpmisc)
library(ggrepel)
library(gridExtra)



#where is the data you want to work on?
setwd("C:/Users/Max/Drive/Random/Megan_Stuff/20220831_DNase/results/endpoint_bedgraph/Subtracted")
#This prefix will determine which files in the selected directory are used. Can accept regex
suffix   <- ".relCntrlSubtract"


filter   <- 15  #ignore this many points on either end of the array
span     <- 20 #number of points to use for smoothing. ~20 generates peaks in linker regions, ~3-5 finds rotational positions
peakSpan <- 4  #coefficient that determines the minimum distance between identified peaks
vertSpan <- 2.5 #how much of the graph to show over the max/min smoothed peaks

dir.create("Plots") #creates a directory to dump plots into
targets  <- list.files(paste0(getwd()),pattern=paste0("*",suffix,"*"),full.names = T, recursive = T)



####################
#paired plots of plus/minus strand
##searches the target directory for bedgraphs of the plus/minus stranded data from a plot and graphs them

bedgraphs <- lapply(targets, function(x){
  read_tsv(x,col_names = c("chromosome", "start","stop","score"), show_col_types = F) %>%
                     filter(between(row_number(),filter,nrow(.)-filter))})

names(bedgraphs) <- gsub(paste0(suffix,".norm.filt.bg"),"",gsub("/","_",list.files(paste0(getwd()),pattern=paste0("*",suffix,"*"), recursive = T)))
#pairs <- unique(gsub("Minus|Plus","",names(bedgraphs)))

plotList <- list()
plotListPeaks <- list()

for (i in 1:length(bedgraphs)){
  kfit <- with(bedgraphs[[i]],
                ksmooth(stop,score, kernel = "normal", bandwidth = span))
  
  bg <- bedgraphs[[i]] %>%
    select(start, stop, score) %>%
    mutate(smoothed = kfit$y) %>%
    pivot_longer(score:smoothed, names_to = "type")
  
  maxPeak <- max(subset(bg, type == "smoothed")[["value"]])
  minPeak <- min(subset(bg, type == "smoothed")[["value"]])
  
  p <- bg %>%
    ggplot(aes(x=stop,y = value)) +
    ggtitle(paste(names(bedgraphs)[i], collapse = " " ))+
    xlab("position")+
    ylab("proportion of cleavage") +
    ylim(minPeak*vertSpan,maxPeak*vertSpan)+
    geom_point(data = subset(bg, type =="score"),
               size = 0.5,alpha = 0.4, color = "grey30") +
    #scale_color_manual(values = c("#4c4e87","#0008ff","#ad5b55","#ff1100"))+
    geom_line(data = subset(bg, type == "smoothed"),
              aes(stop, value), alpha = 0.8, size = 0.8, color = "red")+
    geom_hline(yintercept = 0)+
    stat_peaks(data = subset(bg, type == "smoothed"), span=(span*peakSpan), geom = "point", size = 2, shape = 21, ignore_threshold = .4, color = "black", fill = "red")+
    
    theme_bw()
    
  plotList[[names(bedgraphs)[i]]] <- p
  plotListPeaks[[names(bedgraphs)[i]]] <- p +
    stat_peaks(data = subset(bg, type == "smoothed"), span=(span*peakSpan), geom = "label_repel", 
             ignore_threshold = 0.4, label.padding = .1, label.size = .2, size = 4, box.padding = .4, color = "black", ylim = c(maxPeak*1.5,Inf))
}

grid.arrange(grobs = plotList)
dev.copy(png,paste0(getwd(),"/Plots/NucleasePeaks_noLabs.png"), width = 4000, height = 4000,res = 400)
dev.off()
dev.copy(pdf,paste0(getwd(),"/Plots/NucleasePeaks_noLabs.pdf"), width = 8, height = 8)
dev.off()

grid.arrange(grobs = plotListPeaks)
dev.copy(png,paste0(getwd(),"/Plots/NucleasePeaks_withLabs.png"), width = 4000, height = 4000,res = 400)
dev.off()
dev.copy(pdf,paste0(getwd(),"/Plots/NucleasePeaks_withLabs.pdf"), width = 8, height = 8)
dev.off()





