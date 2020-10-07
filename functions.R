library(tidyverse)
library(purrr)
library(data.table)
library(knitr)

library(mixtools)
library(permute)
library(PRROC)

library(biomaRt)
library(org.Hs.eg.db)

library(ggthemes)
library(gplots)
library(circlize)
library(cowplot)
library(ggbeeswarm)
library(annotate)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggpubr)
library(grid)

theme_Publication <- function(base_size=40, base_family="helvetica") {
  (theme_foundation(base_size=base_size)#, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2,size = 24),
           axis.title.x = element_text(vjust = -0.2,size = 24),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.key.size= unit(0.2, "cm"),
           legend.spacing= unit(0, "cm"),
           legend.title = element_text(face="italic"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}

'%!in%' <- function(x,y)!('%in%'(x,y))

fisher_edges <- function(GENE1,GENE2,data_feed){
  
  temp <- data_feed %>% filter(GENE %in% c(GENE1,GENE2)) %>% unique() %>%
    dplyr::select(GENE,Cell_Line, hit) %>% spread(Cell_Line,hit)
  
  rownames(temp) = temp$GENE
  temp$GENE = NULL
  temp <- t(temp)
  tab <- table(temp[,1],temp[,2])
  
  fish_test <- fisher.test(tab)
  return(fish_test$p.value)
}

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

fix_count <- function(dat) { #made to fix PR curves examples for a couple broken dataframes - scoring in curve breaks for one or two points
  store = c()
  dat <- dat %>% arrange(desc(Precision))
  start = 0
  for (i in 1:nrow(dat)){
    test <- dat$count[i]
    if(test < start){
      store = c(store,i)
      start = test
    }
    start = test
  }
  if (!is_empty(store)){
    dat <- dat[-store,]
  }
  return(dat)
}

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

limit_genes_tab <- function(dat) {
  tab <- table(dat$GENE) %>% sort(decreasing = TRUE)
  tab <- tab[tab > 6]
  #tab <- tab[1:25]
  dat <- dat %>% filter(GENE %in% names(tab)) %>% droplevels()
  return(dat)
}

color.bar <- function(cancer_types) {
  p <- brewer.pal(12, "Paired") #tissue colors, picking the palette in order to standardize colors across
  p2 <- brewer.pal(8, "Dark2")
  p3 <- brewer.pal(9, "Set1")
  p4 <- brewer.pal(8, "Set2")
  p5 <- brewer.pal(12, "Set3")
  
  pal <-  c(p, p2, p3, p4, p5)
  rm(p, p2, p3, p4, p5)
  #All possible tissues (not all are used)
  tissues <-
    c(
      "BLADDER CANCER",
      "BONE CANCER",
      "BRAIN CANCER",
      "BREAST CANCER",
      "COLON/COLORECTAL CANCER",
      "ENDOMETRIAL/UTERINE CANCER",
      "ESOPHAGEAL CANCER",
      "FIBROBLAST",
      "GASTRIC CANCER",
      "HEAD AND NECK CANCER",
      "KIDNEY CANCER",
      "LEUKEMIA",
      "LIVER CANCER",
      "LYMPHOMA",
      "NEUROBLASTOMA",
      "OVARIAN CANCER",
      "PANCREATIC CANCER",
      "SARCOMA",
      "SKIN CANCER",
      "LUNG CANCER",
      "THYROID CANCER",
      "PROSTATE CANCER",
      "RHABDOID",
      "CERVICAL CANCER",
      "IMMORT",
      "IMMORTALIZED",
      "BILE DUCT CANCER",
      "MYELOMA",
      "EYE CANCER",
      "LIPOSARCOMA",
      "NON-CANCEROUS",
      "ADRENAL CANCER",
      "EMBRYONAL CANCER",
      "GALLBLADDER CANCER",
      "PRIMARY CELLS"
    )
  
  colors <- pal[1:length(tissues)] #tissue colors, code is necessary to keep colors consistent
  names(colors) <- tissues
  cancer_color_key <- colors[cells] #assigning the colors to the tissue type
  return(cancer_color_key)
}

