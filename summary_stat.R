###  SUMMARY STATISTICS -- ALL MARKERS ###

## ----- install packages ---------
install.packages("tidyverse", dependencies = TRUE)
install.packages("dplyr", dependencies = TRUE)
install.packages("ggpubr", dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)
install.packages("rstatix", dependencies = TRUE)  # statistical analysis
install.packages("car", dependencies = TRUE)      # MANOVA analysis
install.packages("broom", dependencies = TRUE)    # summary of statistical tests as dataframes
install.packages("gridExtra")                     # multiple plots on the same page
update.packages()

library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("rstatix")       ## for mode    
library("RColorBrewer")
display.brewer.all(colorblindFriendly = TRUE)
library("sf")
library("ggspatial")  ## for scale bar

## ----- summary statistics ------
## loading dataset
setwd("C:/Users/Camilla/Documents/UZH.internship/sum_data/updated")
dat <- read.csv("marker_sum.csv", sep = ',', header = T)
## by marker
marker <- count(group_by(dat, marker), He = mean(He, na.rm = T), HeSD = sd(He, na.rm = T),
                Ar = mean(Ar, na.rm = T), Arsd = sd(Ar, na.rm = T),
                Fst = mean(Fst, na.rm = T),Fstsd = sd(Fst, na.rm=T))
npap <- count(group_by(dat, reference), year = year, marker = marker)  ## study summary

myCols <- brewer.pal(9, "Paired") # Set the color palette (RColorBrewer library)
names(myCols) <- levels(dat$marker) # assign specific colors to markers

## Figure 2A: No studies per year
F2a <- ggplot(npap, mapping = aes(x = year, y = frequency(year), fill = marker)) +
      geom_bar(stat = "identity", show.legend = FALSE) +
      scale_fill_manual("Molecular Marker", values = myCols)+
      guides(fill = FALSE)+
      scale_x_continuous(breaks = seq(1980, 2020, by = 5))+
      scale_y_continuous(breaks = seq(0, 12, by = 2),limits = c(0,12),expand=c(0,0))+
      theme_classic()+
      theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 8, margin = margin(0, 0, 2, 0, "mm"), colour = "black"), 
            axis.text.y = element_text(vjust = 0.5, hjust=1, size = 8, margin = margin(0, 0, 2, 0, "mm"), colour = "black"), 
            axis.line = element_line(size = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title = element_text(hjust = 0.5,face = "bold", size = 8))+
           # panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+  
      labs(x = "Year", y = "Number of studies")
F2a
   # save plot
    setwd("C:/Users/Camilla/Dropbox/FagusGenDiv/manuscript/Figures")
    jpeg(filename = "Fig2a.jpg", width =80, height = 50, res = 600, units = "mm", quality =  100) 
    plot(F2a)
    dev.off()
    
pop <- unique(dat$population)  # No. populations
length(pop)
countries <- unique(dat$country)  # No. countries
countries

npap <- count(group_by(dat, reference), year = year, marker = marker, country = country)
table(npap$country)             # No. studies per country
    write.csv("Nstudies_xcountry.csv", npap, sep = ",")
npop <- count(group_by(dat, population), country= country)
table(npop$country)             # No. populations per country

npop <- count(group_by(dat, reference), marker = marker) %>%  distinct(reference, .keep_all = T)
summary(npop$n)                 # No. populations per study 
hist(npop$n)
get_mode(npop$n)
npop[which.min(npop$n),]        # study with smaller No. pops
npop[which.max(npop$n),]        # study with biggest No. pops

# sample size
ntree <- count(group_by(dat, reference), marker = marker, population = population, N = N) %>%
  distinct(population, .keep_all = T)
hist(ntree$N)               # No. trees per population per study (sample size)
summary(ntree$N)
get_mode(ntree$N)
ntree[which.min(ntree$N),]
ntree[which.max(ntree$N),]
unique(ntree[which(ntree$N<10), "reference"]) # studies with N < 10 (Lander, ind grouped in 3 big groups, Emiliani, 5 just for RFLP analysis)

# No studies per marker type
npap <- count(group_by(dat, reference), year = year, marker = marker)
  df <- as.data.frame(table(npap$marker))
  df <- df %>% mutate(Perc = (Freq * 100)/sum(Freq)) %>% arrange(desc(Perc))
  df$Var1 <- factor(df$Var1, levels = rev(as.character(df$Var1)))
  col <- union(dat$marker, df$Var1)
  names(myCols) <- col

  ## Figure 2B: stacked bar of percentages and number of studies per molecular marker  
F2b <- ggplot(df, aes(x = 1, y = Perc, fill = Var1))+
        geom_bar(stat='identity', position  = "stack", color = "white", size = 0.5)+
        scale_fill_manual(values = myCols)+
        theme_classic() +  coord_flip()+
        guides(fill = FALSE)+                  # legend in Fig 2d
        labs(x = NULL, y =NULL) +
        theme(axis.line = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),axis.title = element_text(hjust = 0.5, color = "#666666"))+
  annotate("text", x = 1, y = 20, label = "40.98%
(25)", size = 2, fontface = "bold")+
  annotate("text", x = 1, y = 53, label = "24.59%
(15)", size = 2, fontface = "bold")+
  annotate("text", x = 1, y = 71, label = "9.83%
(6)", size = 2, fontface = "bold")+
  annotate("text", x = 1, y = 80, label = "8.19%
(5)", size = 2, fontface = "bold")+
  annotate("text", x = 1.34, y = 89, label = "4.91%
(3)", size = 2, fontface = "bold")+
  annotate("text", x = 0.68, y = 95, label = "3.27%
(2)", size = 2, fontface = "bold")+
  annotate("text", x = 1.34, y = 100, label = "1.63%
(1)", size =2, fontface = "bold")+
  geom_segment(aes(y = 84, yend = 93, x = 1.47, xend = 1.47), lineend = "square", size = 0.6)+
  geom_segment(aes(y = 97, yend = 99.5, x = 1.47, xend = 1.47), lineend = "square", size = 0.6)+
  geom_segment(aes(y = 93.5, yend = 96.5, x = 0.53, xend =0.53), lineend = "square", size = 0.6)

  F2b
    jpeg(filename = "Fig2b.jpg", width = 80, height =30, res = 400, units = "mm", quality =  100) 
    plot(F2b)
    dev.off()

## Figure 2C: No studies per country 
map.world<-map_data(map="world")
stud_count <- read.csv("C:/Users/Camilla/Documents/UZH.internship/sum_data/updated/Nstudies_xcountry.csv", sep = ",", header = T)
fagusrange <- st_read("C:/Users/Camilla/Documents/UZH.internship/Fagus_map/Fagus_sylvatica/Fagus_sylvatica_EUFORGEN.shp")

F2c <- ggplot() + 
  geom_map(data=map.world, map=map.world, aes(map_id=region), fill="grey40", alpha = 0.3, colour="white", size=0.15)+
  theme_void()+
  coord_quickmap(ylim=c(30,60), xlim=c(-15,35))+
  geom_sf(fagusrange, mapping =aes(fill= "Fagus sylvatica range"), color = NA, fill = "grey40", alpha = 0.7)+
  geom_point(data=stud_count, aes(x=X, y=Y, size = Nstudies), shape = 21, colour = "#FF9933",fill = "#FF9933")+
  geom_text(data = stud_count, aes(x=X, y=Y,label = Nstudies), size = 5, fontface = "bold")+
  scale_size(range=c(3,15),breaks=c(5,10,15,20), guide = FALSE)+
  guides(fill = FALSE)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  labs(x = "", y = "") +
  ggsn::scalebar(fagusrange, dist = 200, st.size=3, height=0.02, transform = TRUE, model = 'WGS84', dist_unit = "km", location = "bottomleft")
F2c
    #save plot
   jpeg(filename = "Fig2c.jpg", width = 170, height = 160, res = 600, units = "mm", quality =  100) 
   plot(F2c)
  dev.off()

## Figure 2D: summary map of populations, molecular markers and sample size
install.packages("ggsn")
F2d <- ggplot() +
  theme_void()+
        geom_map(data=map.world, map=map.world, aes(map_id=region), fill="grey40", alpha = 0.3, colour="white", size=0.15)+
        coord_quickmap(ylim=c(30,60), xlim=c(-15,40))+
        geom_sf(fagusrange, mapping =aes(fill= "Fagus sylvatica range"), color = NA, fill = "grey40", alpha = 0.7)+
        geom_point(data=dat, aes(x=Longitude, fill = marker,y=Latitude,  size=N), pch = 21, colour = "white",stroke = 1, alpha=I(0.7))+
        scale_fill_manual("Molecular Marker", values = myCols)+
        scale_size(range=c(1,15),name ="Sample size",
                   breaks=c(10,50,100,150,200,300), labels=c(">10",">50",">100",">150",">200",">300"))+
  guides(size = guide_legend(override.aes = list(alpha = 1,shape = 21, fill = "black")), 
         fill = guide_legend(override.aes = list(size = 6, alpha = 1,shape = 22), reverse =T))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.title = element_text(size=12, face = "bold"), 
        legend.text = element_text(size= 8, face = "bold"), 
        legend.spacing.x = unit(0.02, 'cm'),
        legend.position = c(0.885, 0.53))+
  labs(x = " ", y = " ", fill = "Level")+
  ggsn::scalebar(fagusrange, dist = 200, st.size=3, height=0.02, transform = TRUE, model = 'WGS84', dist_unit = "km", location = "bottomleft")

   
F2d
 
# save plot  
  setwd("C:/Users/Camilla/Dropbox/FagusGenDiv/manuscript/Figures")
  jpeg(filename = "Fig2d.jpg", width = 170, height = 160, res = 600, units = "mm", quality =  100) 
  plot(F2d)
  dev.off()

# Fig. S1 - distribution of He across the markers
## density distribution
install.packages("ggridges")
library(ggridges)

FS1 <- ggplot(drop_na(dat, He), aes(He, y = marker, fill = factor(marker), height = ..density..))+
  geom_density_ridges2(scale = 3, stat = "density",alpha = 0.7, size=0.6)+
  scale_fill_manual("Molecular Marker", values = myCols)+
  scale_x_continuous(limits=c(0,1))+
  theme_ridges()+
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 12, margin = margin(0, 0, 0.4, 0, "cm")), 
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 12, margin = margin(0, 0, 0.4, 0, "cm")), 
        axis.title.x = element_text(vjust = 0.5, hjust=0.5, size = 12,face="bold"),
        axis.title.y = element_text(vjust = 0.5, hjust=0.5, size = 12,face="bold"),
        axis.ticks.y=element_blank(),
        legend.position='none')+
  labs(x = "Expected Heterozygosity", y = "Molecular Marker")
FS1
# save plot  
  setwd("C:/Users/Camilla/Dropbox/FagusGenDiv/manuscript/Figures")
  jpeg(filename = "FigS1.jpg", width = 174, height = 95, res = 600, units = "mm", quality =  100)
  plot(FS1)
  dev.off()

## ---- Spearman's correlation of He based on different markers ------ 
## BILELA 2012
bil12 <- filter(dat, reference=="Bilela et al., 2012")
bil12_Hi <- bil12[bil12$marker=="isozymes", "He"] 
bil12_Hn <- bil12[bil12$marker=="nuclear microsatellites", "He"]

cor.test(bil12_Hi, bil12_Hn,          ## Pearson's correlation test He-isozymes vs He-SSR
         alternative = "greater", 
         method = c("spearman"),
         exact = NULL, conf.level = 0.95, continuity = F)
qplot(x = bil12_Hi, y = bil12_Hn, xlab = "He isozymes", ylab = "He nuclear microsatellites")

##CUERVO-ALARCON 2018
cue18 <- filter(dat, reference=="Cuervo-Alarcon et al., 2018")
cue18_Hn <- cue18[cue18$marker=="nuclear microsatellites", 13]
cue18_Hs <- cue18[cue18$marker=="SNP", 13]
cor.test(cue18_Hn, cue18_Hs,          ## Pearson's correlation test He-isozymes vs He-SSR
         alternative = "two.sided", 
         method = c("spearman"),
         exact = NULL, conf.level = 0.95, continuity = F)
qplot(x = cue18_Hn, y = cue18_Hs, xlab = "He nuclear microsatellites", ylab = "He SNP")
cue18$Henm <- cue18[cue18$marker=="nuclear microsatellites", 13]
cue18$Hesn <- cue18[cue18$marker=="SNP", 13]
cue18 = cue18[1:12,] 
ggplot(cue18, aes(x=Henm, y=Hesn))+
  geom_point(size = 2, alpha=1)+
  geom_smooth(method = lm, color="darkblue")+
  theme_bw()+
  labs(x ="He nuclear microsatellites", y = "He SNP", title = "He Cuervo-Alarcon et al., 2018")

## MULLER 2018
mul18 <- filter(dat, reference=="Muller et al., 2018")
mul18$Henm <- mul18[mul18$marker=="nuclear microsatellites", 13]
mul18$Hesn <- mul18[mul18$marker=="SNP", 13]
mul18 = mul18[1:24,]
cor.test(mul18$Henm,mul18$Hesn,          ## Pearson's correlation test He-SNP vs He-SSR
         alternative = "greater", 
         method = c("spearman"),
         exact = NULL, conf.level = 0.95, continuity = F)
ggplot(mul18, aes(x=Henm, y=Hesn))+
  geom_point(size = 2, alpha=1)+
  geom_smooth(method = lm, color="darkblue")+
  theme_bw()+
  labs(x ="He nuclear microsatellites", y = "He SNP", title = "He Muller et al., 2018")

## PAFFETTI 2012
paf12 <- filter(dat, reference=="Paffetti et al., 2012")
paf12_Hn <- paf12[paf12$marker=="nuclear microsatellites", 13]
paf12_Hr <- paf12[paf12$marker=="RAPD", 13]
qplot(x = paf12_Hn, y = paf12_Hr, xlab = "He nuclear microsatellites", ylab = "He RAPD")

## ---- summary of genetic metrics -----
### allelic richness 
  summary(dat$Ar)
  shapiro.test(dat$Ar) 
  kruskal.test(Ar~marker, data = dat) 
  
  ggplot(drop_na(dat, Ar), aes(x = marker, y = Ar, fill = marker)) + ## removing markers with no data
    geom_boxplot(width = 0.4) +
    geom_jitter(color="black", size = 1, alpha = 1)+theme(legend.position = c(0.25, 2))+
    theme(legend.position = c(0.25, 4),panel.border = element_rect(colour = "black", fill = NA, size = 1))+
    annotate("text", x = 1.5, y = 27, label = "Kruskal-Wallis chi squared = 592.54, p = < 2.2e-16", size = 3)

  summary(with(dat, Ar[marker == "isozymes"]))
  sd(dat[which(dat$marker == "isozymes"), ]$Ar, na.rm = T)
  summary(with(dat, Ar[marker == "nuclear microsatellites"]))
  sd(dat[which(dat$marker == "nuclear microsatellites"), ]$Ar, na.rm = T)
 
  ##He
  summary(dat$He)
    ## average and sd, He from nuclear markers
    summary(with(dat, He[marker == "isozymes" | marker == "nuclear microsatellites"| marker == "SNP" | marker == "RAPD" | marker == "AFLP"]))
    sd(dat[which(marker == "isozymes" | marker == "nuclear microsatellites"| marker == "SNP" | marker == "RAPD" | marker == "AFLP"),]$He, na.rm = T)
    ## average and sd, He from chloroplast markers 
    summary(with(dat, He[marker == "cpDNA microsatellites" | marker == "cpDNA PCR-RFLP"| marker == "cpDNA RAPD" | marker == "cpDNA SNP"]))
    sd(dat[which(marker == "cpDNA microsatellites" | marker == "cpDNA PCR-RFLP"| marker == "cpDNA RAPD" | marker == "cpDNA SNP"),]$He, na.rm = T)
    
  kruskal.test(He~marker, data = dat)
  
  ggplot(dat[which(dat$marker != "cpDNA PCR-RFLP"&dat$marker != "cpDNA RAPD"),], ## removing markers with no data
         aes(x =marker, y = He, fill = marker)) +
    geom_boxplot()+
    geom_jitter(color="black", size = 1, alpha = 1)+theme(legend.position = c(0.25, 2))+
    labs(y = "Expected Heterozygosity", x = "Molecular Marker")
 
 ## Fst
  summary(dat[grep("cpDNA", dat$marker),]$Fst) ## Fst for cpDNA molecular markers
  sd(dat[grep("cpDNA", dat$marker),]$Fst, na.rm = T)
  summary(dat[-grep("cpDNA", dat$marker),]$Fst) ## Fst for nuclear molecular markers
  sd(dat[-grep("cpDNA", dat$marker),]$Fst, na.rm = T)
  kruskal.test(Fst ~ marker, data = dat)
  
  ggplot(drop_na(dat, Fst), aes(x = marker, y = Fst)) + 
    geom_boxplot(width = 0.4) +
    theme_classic2()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))+
    scale_y_continuous(limits = c(0,1),expand=c(0,0))+
    labs(x = "Molecular Marker", y = "Fst")+
    annotate("text", x = 1.5, y = 0.9, label = "Kruskal-Wallis chi squared = 592.54, p = < 2.2e-16", size = 3)

### Fis 
  Fis <- c(-0.044,-0.013,0.003,0.015,0.01684,0.02,0.084,0.224,
          0.09,0.005,0.003,-0.003,0.016,0.024,0.035,0.064,
          0.038,0.022,0.042,0.05,0.065,0.239)
  mean(Fis)
  sd(Fis)
  
  ggplot(dat[which(dat$marker != "cpDNA PCR-RFLP"&dat$marker != "cpDNA RAPD"),], 
         aes(x =marker, y = Fst, fill = marker)) +
    geom_boxplot()+
    geom_jitter(color="black", size = 1, alpha = 1)+theme(legend.position = c(0.25, 2))+
    labs(y = "Fst", x = "Molecular Marker")

## ---- SSR ----
SSR <- filter(dat, marker == "nuclear microsatellites")
SSRsum <- SSR %>% 
  group_by(reference) %>% 
  summarise(HeSD = sd(He, na.rm = T), He = mean(He, na.rm = T)) 
summary(SSRsum) ## average SD for the He 

   ## detailed analysis in SSR script

## ---- SNP ----
    ## analysis Capblancq dataset
    setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SNP/Capblancq_2020_dataset/vcftools_analysis/all_pop")
    #genetic diversity per individual
    het <- read.table("het_ind.het", header = T)    
    str(het)
    het$E.HET <- het$N_SITES-het$E.HOM. # Expected No. heterozygous sites per ind
    het$He <-  het$E.HET/het$N_SITES  # expected heterozygosity per ind
    
    # + populations, latitude and longitude
    info <- read.csv("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SNP/Capblancq_2020_dataset/Infos_Individus.csv", header = T, sep = ",")
    cap.pop <- info[, c("Ind", "pop", "X", "Y")]
    hist(cap.pop$He)
    cap.pop$INDV <- cap.pop$Ind
    cap.pop <- cap.pop[, -1]
    cap.pop <- full_join(cap.pop, het, by = NULL, copy = T, keep = F)
    summary(cap.pop$He)
    sd(cap.pop$He, na.rm = T)
    cap.pop$He_mean <- cap.pop %>% group_by(pop) %>% mutate(He_mean = mean(He))
    
    ## He per population
    ggplot(cap.pop, aes(pop,He))+
      theme_classic()+
      geom_boxplot()+
      theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
      labs(title = "heterozygosity per population")
    setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SNP/Capblancq_2020_dataset")
    write.csv(cap.pop, "Capblancq_He.csv", sep = ",")
        ## dataset added to main data

### SNP data
SNP <- subset(dat, marker == "SNP")
head(SNP)
hist(SNP$He)
summary(SNP$He)
sd(SNP$Nuc_div, na.rm = T)
summary(SNP$Nuc_div)
sd(SNP$He, na.rm = T)

SNPlong <- gather(SNP, key = "param", value = "value", -marker, -reference, -year, 
                  -country, -population, -Latitude, -Longitude, -altitude, -N, -Nloci)
SNPlong <- SNPlong[which(SNPlong$param == "He"|SNPlong$param == "Nuc_div"),]
SNPlong <- drop_na(SNPlong, value)

snp.he <- drop_na(SNP, He)        ## Heterozygosity
ggplot()+ geom_boxplot(snp.he, mapping = aes(x = reference, y = He))
hist(snp.he$He)
snp.nd <- drop_na(SNP, Nuc_div)   ## Nucleotide diversity
ggplot()+ geom_boxplot(snp.nd, mapping = aes(x = reference, y = Nuc_div))

