## analysis of SSR datasets ###
#### install packages and libraries ####

.libPaths("C:/Program Files/R/R-3.6.3patched/library")
install.packages("tidyverse")
install.packages("gplots")
install.packages("dplyr")
install.packages("pegas")   ## to read Genepop file 
install.packages("adegenet")
install.packages("genepop")
install.packages("hierfstat")

library("dplyr")       ## for tidying data
library("tidyr")
library("ggplot2")
library("pegas")       ## for genetic analyses
library("adegenet")
library("hierfstat")
library("ggpubr")       ## for ggarrange
library("PopGenReport") ## allel.rich function
library("maps")         ## maps
library("sf")           ## reading shapefile
library("ggrepel")
library("RColorBrewer")
library("viridis")
library("ggspatial")  ## for scale bar
library("rgdal")      ## for interpolation
library("gstat")
library("sp")

## --------- Lander 2011 (no populations specified in the original dataset--- dividing by area)----------
setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SSR/Lander_2011_dataset")
lan11 <- read.loci("Lander_SSRdataset.csv", header = T, loci.sep = ",", allele.sep = "/", col.pop = 3, col.loci = 4:16)
lan11_gin <- loci2genind(lan11, ploidy = 2, na.alleles = "0", unphase = T)
summary(lan11_gin)

## converting to genepop 
lan11_gpop = genind2genpop(lan11_gin, quiet = T, process.other = T)
lan11_gin$loc.n.all
lan11_gpop$tab                                             ## N of individuals per locus.allele per population
lan11_allfr = as.data.frame(makefreq(lan11_gpop))                        ## allelic frequencies 
  write.csv(lan11_allfr, "lan11_allfr.csv")
lan11_he = data.frame(t(basic.stats(lan11_gin)$Hs))        ## He per locus per population
lan11_ar <- data.frame(allel.rich(lan11_gin, min.alleles = 25)$mean.richness)  ## allelic richness (N = 50)

lan_split = seppop(lan11_gin, lan11_gin$pop)               # separating populations to find Na for each one
East= c(nAll(lan_split$East, onlyObserved = T))
South= c(nAll(lan_split$South, onlyObserved = T))
West = c(nAll(lan_split$West, onlyObserved = T))
lan_na = data.frame(East, South, West)                     
lan_na = t(lan_na)                                         ## Na per population per locus

## --------- Lefevre 2012----------------------------
setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SSR/Lefevre_2012_dataset")
lef12 <- read.loci("Lefevre_SSRdata.csv", header = T, loci.sep = ",", allele.sep = "/", col.pop = 2, col.loci = 3:18)
lef12_gin <- loci2genind(lef12, ploidy = 2, na.alleles = "0", unphase = T)
summary(lef12_gin)

## converting to genepop 
lef12_gpop = genind2genpop(lef12_gin, quiet = T, process.other = T)

lef12_gpop$tab                                              ## N of individuals per locus.allele per population
lef12_allfr = as.data.frame(makefreq(lef12_gpop))                          ## allelic frequencies
  write.csv(lef12_allfr, "lef12_allfr.csv")
lef12_he <- data.frame(t(basic.stats(lef12_gin)$Hs))         ## expected heterozygosity (replaced to the reported one -- only average)
lef12_ar <- data.frame(allel.rich(lef12_gin, min.alleles = 25)$mean.richness)   ## allelic richness (N = 50)

lef_split = seppop(lef12_gin, lef12_gin$pop)               ## separating populations to find Na for each one
pop1= c(nAll(lef_split$pop1, onlyObserved = T))
pop2= c(nAll(lef_split$pop2, onlyObserved = T))
pop3= c(nAll(lef_split$pop3, onlyObserved = T))
pop4= c(nAll(lef_split$pop4, onlyObserved = T))
lef_na = data.frame(pop1,pop2,pop3, pop4)                
lef_na = t(lef_na)                                         ## Na per locus per population

## --------- Piotti 2012-----------------------------
setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SSR/Piotti_2012_dataset")
pio12 <- read.loci("Piotti_SSRdataset.csv", header = T, loci.sep = ",", allele.sep = "/", col.pop = 2, col.loci = 3:6, na.strings = "NA")
pio12_gin <- loci2genind(pio12, ploidy = 2, na.alleles = "NA", unphase = T)
summary(pio12_gin)

## converting to genepop 
pio12_gpop = genind2genpop(pio12_gin, quiet = T, process.other = T)

pio12_gpop$tab                                              ## N of individuals per locus.allele per population
pio12_allfr = as.data.frame(makefreq(pio12_gpop))                       ## allelic frequencies
  write.csv(pio12_allfr, "pio12_allfr.csv")
pio12_he = data.frame(t(basic.stats(pio12_gin)$Hs))        ## expected heterozygosity (replaced to the reported one -- only average)
pio12_ar <- data.frame(allel.rich(pio12_gin, min.alleles = 25)$mean.richness)  ## Allelic richness rarefacted to standard N = 50

pio_split = seppop(pio12_gin)
D1 = c(nAll(pio_split$D1, onlyObserved = T))
D2= c(nAll(pio_split$D2, onlyObserved = T))
SB= c(nAll(pio_split$SB, onlyObserved = T))
VE= c(nAll(pio_split$VE, onlyObserved = T))
pio_na = data.frame(D1,D2,SB,VE)
pio_na = t(pio_na)                                          ## Na per locus per population

## --------- Gauzere 2013----------------------------
setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SSR/Gauzere_2013_dataset")
gau13_csv <- read.loci("GeneticData.csv",header = T, loci.sep = ",", allele.sep = "/", col.pop = 2, col.loci = 3:15, na.strings = "1")
gau13_gin = loci2genind(gau13_csv, ploidy = 2, na.alleles = "1", unphase = T)
summary(gau13_gin)

# create genepop file
gau13_gpop = genind2genpop(gau13_gin, quiet = T, process.other = T)  

gau13_gpop$tab                                            ## N of individuals per locus.allele per population
gau13_allfr = as.data.frame(makefreq(gau13_gpop))                        ## allelic frequencies
  write.csv(gau13_allfr, "gau13_allfr.csv")
gau13_he = data.frame(t(basic.stats(gau13_gin)$Hs))       ## Expected Heterozygosity for the populations
gau13_ar <- data.frame(allel.rich(gau13_gin, min.alleles = 25)$mean.richness)   ## Allelic richness rarefacted to standard N = 50

nAll(gau13_gin)                                           ## Na per population
gau_split = seppop(gau13_gin, gau13_gin$pop)
N1 = c(nAll(gau_split$N1, onlyObserved = T))
N2 =  c(nAll(gau_split$N2, onlyObserved = T))
N4 =  c(nAll(gau_split$N4, onlyObserved = T))
gau_na = data.frame(N1, N2,N4)
gau_na = t(gau_na)                                        ## Na per locus per population

## --------- DeLafontaine 2013  -------------------------------------
setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SSR/DeLafontaine_2013_dataset")
del13 <- read.loci("SSR_Dataset.csv", header = T, loci.sep = ",", allele.sep = "/", col.pop = 1, col.loci = 3:18, na.strings = "0")
del13_gin <- loci2genind(del13, ploidy = 2, na.alleles = "0", unphase = T)
summary(del13_gin)

##converting to genepop file
del13_gpop = genind2genpop(del13_gin, quiet = T, process.other = T)
del13_gpop$other
                                              ## N of individuals per locus.allele per population
del13_allfr = makefreq(del13_gpop)                          ## allelic frequencies
  write.csv(del13_allfr, "del13_allfr.csv")
del13_he = data.frame(t(basic.stats(del13_gin)$Hs))         ## expected heterozygosity (slightly different from the reported by the study)
del13_ar <- data.frame(allel.rich(del13_gin, min.alleles = 25)$mean.richness)  ## Allelic richness rarefacted to standard N = 50

del13_ar <- drop_na(del13_ar)                               # 3 populations Ar is NaN (COM, GHE, MOI) ---> removed

del_ho = as.data.frame(basic.stats(del13_gin,diploid=TRUE,digits=4)$Ho)   ## Ho per locus per population
del_he = as.data.frame(basic.stats(del13_gin,diploid=TRUE,digits=4)$Hs)   ## He per locus per population

del_split = seppop(del13_gin)
length(del_split)
 for(loc.n.all in 1:length(del_split)){
   
 }  ##?

AGE= c(nAll(del_split$AGE, onlyObserved = T))                   
AIG= c(nAll(del_split$AIG, onlyObserved = T))
ARG= c(nAll(del_split$ARG, onlyObserved = T))
AUB= c(nAll(del_split$AUB, onlyObserved = T))
BAG= c(nAll(del_split$BAG, onlyObserved = T))
BAI= c(nAll(del_split$BAI, onlyObserved = T))
BAS= c(nAll(del_split$BAS, onlyObserved = T))
BAU= c(nAll(del_split$BAU, onlyObserved = T))
BRI= c(nAll(del_split$BRI, onlyObserved = T))
CAZ= c(nAll(del_split$CAZ, onlyObserved = T))
CHI= c(nAll(del_split$CHI, onlyObserved = T))
CIRA= c(nAll(del_split$CIRA, onlyObserved = T))
CIRB= c(nAll(del_split$CIRB, onlyObserved = T))
CIRC= c(nAll(del_split$CIRC, onlyObserved = T))
COL= c(nAll(del_split$COL, onlyObserved = T))
COM= c(nAll(del_split$COM, onlyObserved = T))
CRE= c(nAll(del_split$CRE, onlyObserved = T))
CUR= c(nAll(del_split$CUR, onlyObserved = T))
ELS= c(nAll(del_split$ELS, onlyObserved = T))
EPI= c(nAll(del_split$EPI, onlyObserved = T))
ESC= c(nAll(del_split$ESC, onlyObserved = T))
FOU= c(nAll(del_split$FOU, onlyObserved = T))
GAR= c(nAll(del_split$GAR, onlyObserved = T))
GAVR= c(nAll(del_split$GAVR, onlyObserved = T))
GHE= c(nAll(del_split$GHE, onlyObserved = T))
GOB= c(nAll(del_split$GOB, onlyObserved = T))
HAU= c(nAll(del_split$HAU, onlyObserved = T))
HAY= c(nAll(del_split$HAY, onlyObserved = T))
HES= c(nAll(del_split$HES, onlyObserved = T))
HET= c(nAll(del_split$HET, onlyObserved = T))
ISS= c(nAll(del_split$ISS, onlyObserved = T))
JAU= c(nAll(del_split$JAU, onlyObserved = T))
JOS= c(nAll(del_split$JOS, onlyObserved = T))
JOY= c(nAll(del_split$JOY, onlyObserved = T))
LAG= c(nAll(del_split$LAG, onlyObserved = T))
LAS= c(nAll(del_split$LAS, onlyObserved = T))
LAV= c(nAll(del_split$LAV, onlyObserved = T))
LEO= c(nAll(del_split$LEO, onlyObserved = T))
LOU= c(nAll(del_split$LOU, onlyObserved = T))
LUR= c(nAll(del_split$LUR, onlyObserved = T))
MAD= c(nAll(del_split$MAD, onlyObserved = T))
MAS= c(nAll(del_split$MAS, onlyObserved = T))
MAU= c(nAll(del_split$MAU, onlyObserved = T))
MIS= c(nAll(del_split$MIS, onlyObserved = T))
MOI= c(nAll(del_split$MOI, onlyObserved = T))
MOU= c(nAll(del_split$MOU, onlyObserved = T))
NOH= c(nAll(del_split$NOH, onlyObserved = T))
NOI= c(nAll(del_split$NOI, onlyObserved = T))
OLO= c(nAll(del_split$OLO, onlyObserved = T))
OME= c(nAll(del_split$OME, onlyObserved = T))
ORI= c(nAll(del_split$ORI, onlyObserved = T))
PAI= c(nAll(del_split$PAI, onlyObserved = T))
PRA= c(nAll(del_split$PRA, onlyObserved = T))
PYM= c(nAll(del_split$PYM, onlyObserved = T))
REY= c(nAll(del_split$REY, onlyObserved = T))
ROL= c(nAll(del_split$ROL, onlyObserved = T))
ROQA= c(nAll(del_split$ROQA, onlyObserved = T))
ROQB= c(nAll(del_split$ROQB, onlyObserved = T))
SYM= c(nAll(del_split$SYM, onlyObserved = T))
VAC= c(nAll(del_split$VAC, onlyObserved = T))
VAL= c(nAll(del_split$VAL, onlyObserved = T))
VAU= c(nAll(del_split$VAU, onlyObserved = T))
VIL= c(nAll(del_split$VIL, onlyObserved = T))
VIR= c(nAll(del_split$VIR, onlyObserved = T))
XOK = c(nAll(del_split$XOK, onlyObserved = T))

del_na = data.frame(AGE,AIG,ARG,AUB,BAG,BAI,BAS,BAU,BRI,CAZ,CHI,CIRA,CIRB,CIRC,COL,COM,CRE,CUR,ELS,EPI,ESC,FOU,
                    GAR,GAVR,GHE,GOB,HAU,HAY,HES,HET,ISS,JAU,JOS,JOY,LAG,LAS,LAV,LEO,LOU,LUR,MAD,MAS,MAU,MIS,
                    MOI,MOU,NOH,NOI,OLO,OME,ORI,PAI,PRA,PYM,REY,ROL,ROQA,ROQB,SYM,VAC,VAL,VAU,VIL,VIR,XOK)             
del_na = t(del_na)                                   ## Na per locus per population
  write.csv(del_na, file = "del13_Na.csv")
  
## --------- Rajendra 2014----------------------------
setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SSR/Rajendra_2014_dataset")
raj14_csv <- read.loci("Rajendra_2014_dataset.csv",header = T, loci.sep = ",", allele.sep = "/", col.pop = 2, col.loci = 3:11,na.strings ="") 
raj14_gin = loci2genind(raj14_csv, ploidy = 2, unphase = T, na.alleles = " ")    #creating genind file
summary(raj14_gin)
raj14_gpop = genind2genpop(raj14_gin, quiet = T, process.other = T)              #creating genepop file

raj14_gpop$tab                                                                   ## N of individuals per locus.allele per population
raj14_af = as.data.frame(makefreq(raj14_gpop))                                                 # allelic frequencies                                              
  write.csv(raj14_af, "raj_allfreq.csv")
raj_he = as.data.frame(basic.stats(raj14_gin,diploid=TRUE,digits=4)$Hs)          # expected heterozygosity per locus per population
raj14_ar <- data.frame(allel.rich(raj14_gin, min.alleles = 25)$mean.richness)    # Allelic richness (N = 50)

raj_split = seppop(raj14_gin)
AEW19 = c(nAll(raj_split$AEW19, onlyObserved = T))
AEW23= c(nAll(raj_split$AEW23, onlyObserved = T))
AEW39= c(nAll(raj_split$AEW39, onlyObserved = T))
AEW4= c(nAll(raj_split$AEW4, onlyObserved = T))
AEW40= c(nAll(raj_split$AEW40, onlyObserved = T))
AEW5= c(nAll(raj_split$AEW5, onlyObserved = T))
AEW6= c(nAll(raj_split$AEW6, onlyObserved = T))
AEW7= c(nAll(raj_split$AEW7, onlyObserved = T))
AEW8= c(nAll(raj_split$AEW8, onlyObserved = T))
AEW9= c(nAll(raj_split$AEW9, onlyObserved = T))
HEW10= c(nAll(raj_split$HEW10, onlyObserved = T))
HEW11= c(nAll(raj_split$HEW11, onlyObserved = T))
HEW12= c(nAll(raj_split$HEW12, onlyObserved = T))
HEW18= c(nAll(raj_split$HEW18, onlyObserved = T))
HEW21= c(nAll(raj_split$HEW21, onlyObserved = T))
HEW5= c(nAll(raj_split$HEW5, onlyObserved = T))
HEW6= c(nAll(raj_split$HEW6, onlyObserved = T))
HEW7= c(nAll(raj_split$HEW7, onlyObserved = T))
HEW8= c(nAll(raj_split$HEW8, onlyObserved = T))
HEW9= c(nAll(raj_split$HEW9, onlyObserved = T))
SEW30= c(nAll(raj_split$SEW30, onlyObserved = T))
SEW31= c(nAll(raj_split$SEW31, onlyObserved = T))
SEW38= c(nAll(raj_split$SEW38, onlyObserved = T))
SEW4= c(nAll(raj_split$SEW4, onlyObserved = T))
SEW46= c(nAll(raj_split$SEW46, onlyObserved = T))
SEW5= c(nAll(raj_split$SEW5, onlyObserved = T))
SEW6= c(nAll(raj_split$SEW6, onlyObserved = T))
SEW7= c(nAll(raj_split$SEW7, onlyObserved = T))
SEW8= c(nAll(raj_split$SEW8, onlyObserved = T))
SEW9= c(nAll(raj_split$SEW9, onlyObserved = T))

raj_na = data.frame(AEW19,AEW23,AEW39,AEW4,AEW40,AEW5,AEW6,AEW7,AEW8,AEW9,HEW10,HEW11,
                    HEW12,HEW18,HEW21,HEW5,HEW6,HEW7,HEW8,HEW9,SEW30,SEW31,SEW38,SEW4,
                    SEW46,SEW5,SEW6,SEW7,SEW8,SEW9)                              ## Na for each population
raj_na = t(raj_na)

## --------- Sjolund 2014 (2 df: GB populations and Continental populations (Fr, Ge, It) "others")-----
setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SSR/Sjolund_dataset")
sjo14_gb <- read.loci("Sjolund_SSR_GB.csv",header = T, loci.sep = ",", allele.sep = "/", col.pop = 3, col.loci = 8:20)
sjo14_ot <- read.loci("Sjolund_SSR_others.csv",header = T, loci.sep = ",", allele.sep = "/", col.pop = 3, col.loci = 10:19)
sjo14_gbgin = loci2genind(sjo14_gb, ploidy = 2, unphase = T, na.alleles = "NA")            #creating genind files
sjo14_otgin = loci2genind(sjo14_ot, ploidy = 2, unphase = T, na.alleles = "NA")

sjo14_gbpop = genind2genpop(sjo14_gbgin, quiet = T, process.other = T)                     # create genepop files
sjo14_otpop = genind2genpop(sjo14_otgin, quiet = T, process.other = T)

sjo14_gbgin$tab
sjo14_otgin$tab
sjo14_gbgin$pop
sjo14_gbaf = as.data.frame(makefreq(sjo14_gbpop))                                                         # allelic frequencies
sjo14_otaf = as.data.frame(makefreq(sjo14_otpop))
  write.csv(sjo14_gbaf, "sjo14_allfr_gb.csv")
  write.csv(sjo14_otaf, "sjo14_allfr_eu.csv")

sjo14_gbhe = data.frame(t(basic.stats(sjo14_gbgin)$Hs))                                    # expected heterozygosity
  write.csv(sjo14_gbhe, "sjo_gb_he.csv")
sjo14_othe = data.frame(t(basic.stats(sjo14_otgin)$Hs))
  write.csv(sjo14_othe, "sjo_ot_he.csv")

sjo14_gbar <- data.frame(allel.rich(sjo14_gbgin, min.alleles = 25)$mean.richness)          # allelic richness (N = 50)
sjo14_otar <- data.frame(allel.rich(sjo14_otgin,min.alleles = 25)$mean.richness)

sjo14_gb_split = seppop(sjo14_gbgin)
`Applecross Wood` = c(nAll(sjo14_gb_split$`Applecross Wood`, onlyObserved = T))
`Barons Haugh`= c(nAll(sjo14_gb_split$`Barons Haugh`, onlyObserved = T))
`Bedford Purleius`= c(nAll(sjo14_gb_split$`Bedford Purleius`, onlyObserved = T))
`Beech Hill Wood`= c(nAll(sjo14_gb_split$`Beech Hill Wood`, onlyObserved = T))
`Blean Woods`= c(nAll(sjo14_gb_split$`Blean Woods`, onlyObserved = T))
`Bridford Wood`= c(nAll(sjo14_gb_split$`Bridford Wood`, onlyObserved = T))
`Buckholt Wood`= c(nAll(sjo14_gb_split$`Buckholt Wood`, onlyObserved = T))
`Burnham Beeches`= c(nAll(sjo14_gb_split$`Burnham Beeches`, onlyObserved = T))
`Carstramon Wood`= c(nAll(sjo14_gb_split$`Carstramon Wood`, onlyObserved = T))
`Clerkhill Wood`= c(nAll(sjo14_gb_split$`Clerkhill Wood`, onlyObserved = T))
`Craig Wood`= c(nAll(sjo14_gb_split$`Craig Wood`, onlyObserved = T))
`Cwm Clydach (east)`= c(nAll(sjo14_gb_split$`Cwm Clydach (east)`, onlyObserved = T))
`Cwm Clydach (west)`= c(nAll(sjo14_gb_split$`Cwm Clydach (west)`, onlyObserved = T))
`Denny Wood`= c(nAll(sjo14_gb_split$`Denny Wood`, onlyObserved = T))
`Devachoys Wood`= c(nAll(sjo14_gb_split$`Devachoys Wood`, onlyObserved = T))
`Drumneil House`= c(nAll(sjo14_gb_split$`Drumneil House`, onlyObserved = T))
`Dunnottar Wood`= c(nAll(sjo14_gb_split$`Dunnottar Wood`, onlyObserved = T))
`Ecclesall Woods`= c(nAll(sjo14_gb_split$`Ecclesall Woods`, onlyObserved = T))
`Felbrigg Great Wood`= c(nAll(sjo14_gb_split$`Felbrigg Great Wood`, onlyObserved = T))
`Friary Wood`= c(nAll(sjo14_gb_split$`Friary Wood`, onlyObserved = T))
`Gelt Wood`= c(nAll(sjo14_gb_split$`Gelt Wood`, onlyObserved = T))
`Golitha Wood`= c(nAll(sjo14_gb_split$`Golitha Wood`, onlyObserved = T))
`Greenfield Copse`= c(nAll(sjo14_gb_split$`Greenfield Copse`, onlyObserved = T))
`Hembury Wood`= c(nAll(sjo14_gb_split$`Hembury Wood`, onlyObserved = T))
`Kinnoul Hill Woodland Park`= c(nAll(sjo14_gb_split$`Kinnoul Hill Woodland Park`, onlyObserved = T))
`Lady Park Wood`= c(nAll(sjo14_gb_split$`Lady Park Wood`, onlyObserved = T))
`Lullingstone Country Park`= c(nAll(sjo14_gb_split$`Lullingstone Country Park`, onlyObserved = T))
`Mabie Forest`= c(nAll(sjo14_gb_split$`Mabie Forest`, onlyObserved = T))
`Monk Wood`= c(nAll(sjo14_gb_split$`Monk Wood`, onlyObserved = T))
`Plora Wood`= c(nAll(sjo14_gb_split$`Plora Wood`, onlyObserved = T))
`Savernake Forest`= c(nAll(sjo14_gb_split$`Savernake Forest`, onlyObserved = T))
`Seckley Wood`= c(nAll(sjo14_gb_split$`Seckley Wood`, onlyObserved = T))
`Strid Wood`= c(nAll(sjo14_gb_split$`Strid Wood`, onlyObserved = T))
Talhenbont= c(nAll(sjo14_gb_split$Talhenbont, onlyObserved = T))
`Tan-y-Coed`= c(nAll(sjo14_gb_split$`Tan-y-Coed`, onlyObserved = T))
`Tongue Wood`= c(nAll(sjo14_gb_split$`Tongue Wood`, onlyObserved = T))
`Two Mile Bottom`= c(nAll(sjo14_gb_split$`Two Mile Bottom`, onlyObserved = T))
`Wallington East Woods`= c(nAll(sjo14_gb_split$`Wallington East Woods`, onlyObserved = T))
`Wealden Edge Hangars`= c(nAll(sjo14_gb_split$`Wealden Edge Hangars`, onlyObserved = T))
`Wychwood (Cornbury Park)`= c(nAll(sjo14_gb_split$`Wychwood (Cornbury Park)`, onlyObserved = T))
`Wytham Wood`= c(nAll(sjo14_gb_split$`Wytham Wood`, onlyObserved = T))
`Yellowcraig Wood`= c(nAll(sjo14_gb_split$`Yellowcraig Wood`, onlyObserved = T))

sjo_gbna = data.frame(`Applecross Wood`, `Barons Haugh`, `Bedford Purleius`, `Beech Hill Wood`, `Blean Woods`, `Bridford Wood`, `Buckholt Wood`, 
                      `Burnham Beeches`, `Carstramon Wood`, `Clerkhill Wood`, `Craig Wood`, `Cwm Clydach (east)`, `Cwm Clydach (west)`, `Denny Wood`, 
                      `Devachoys Wood`, `Drumneil House`, `Dunnottar Wood`, `Ecclesall Woods`, `Felbrigg Great Wood`, `Friary Wood`, `Gelt Wood`, 
                      `Golitha Wood`, `Greenfield Copse`, `Hembury Wood`, `Kinnoul Hill Woodland Park`, `Lady Park Wood`, `Lullingstone Country Park`,
                      `Mabie Forest`, `Monk Wood`, `Plora Wood`, `Savernake Forest`, `Seckley Wood`, `Strid Wood`, Talhenbont, `Tan-y-Coed`, `Tongue Wood`, 
                      `Two Mile Bottom`, `Wallington East Woods`, `Wealden Edge Hangars`, `Wychwood (Cornbury Park)`, `Wytham Wood`, `Yellowcraig Wood`)
sjo_gbna = t(sjo_gbna)                                                   # Na per locus per population

sjo_ot_split = seppop(sjo14_otgin)
France  = c(nAll(sjo_ot_split$France, onlyObserved = T))
Germany = c(nAll(sjo_ot_split$Germany, onlyObserved = T))
Italy = c(nAll(sjo_ot_split$Italy, onlyObserved = T))
sjo_otna = data.frame(France, Germany, Italy)
sjo_otna = t(sjo_otna)                                                   # Na per locus per population

## --------- Cvrckova 2017 ------
setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SSR/Cvrckova_2017dataset")
cvr17 <- read.loci("Cvr17_gen.csv", header = T, loci.sep = ",", allele.sep = "/", col.pop = 2, col.loci = 3:14, na.strings = "0/0")
cvr17_gin <- loci2genind(cvr17, ploidy = 2, na.alleles = "0", unphase = T)
summary(cvr17_gin)

## converting to genepop 
cvr17_gpop = genind2genpop(cvr17_gin, quiet = T, process.other = T)

cvr17_gpop$tab                                             ## N of individuals per locus.allele per population
cvr17_allfr = as.data.frame(makefreq(cvr17_gpop))                        ## allelic frequencies
  write.csv(cvr17_allfr, "cvr17_allfr.csv")
cvr17_he = data.frame(t(basic.stats(cvr17_gin)$Hs))        ## He per locus per population
cvr17_ar <- data.frame(allel.rich(cvr17_gin, min.alleles = 25)$mean.richness)  ## allelic richness (N = 50)

## Na for each locus for each population TO COMPUTE

## --------- Oddou-Muratorio 2018 -------------------
setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SSR/Oddou-Muratorio_2018_dataset")  
odd18_loci = read.loci("Oddou_all_loci.csv", loci.sep = ",", allele.sep = "/", col.loci= 3:15, col.pop = 2, row.names = 1, na.strings = c("0/0", "1"))
odd18_gin <- loci2genind(odd18_loci, ploidy = 2, na.alleles = "1", unphase = T)      # genind file 
summary(odd18_gin) ##  %of missing data > 3%

# genepop file 
odd18_gpop = genind2genpop(odd18_gin, quiet = T, process.other = T)    
  
odd18_gin$tab                                                                                   ## numbver of individuals per locus.allele per pop
odd18_allfreq = as.data.frame(makefreq(odd18_gpop))                                                           ## allelic frequencies
  write.csv(odd18_allfreq, file = "odd18_allfreq.csv", na = "NA", row.names = T, col.names = T) 
odd18_he = data.frame(t(basic.stats(odd18_gin)$Hs))                                             ## Expected Heterozygosity for each locus and each pop
odd18_ar <- data.frame(allel.rich(odd18_gin,min.alleles = 25)$mean.richness)                    # allelic richness (N = 50)

odd_split = seppop(odd18_gin, odd18_gin$pop)
N1 = c(nAll(odd_split$N1, onlyObserved = T))
N2 = c(nAll(odd_split$N2, onlyObserved = T))
N4  = c(nAll(odd_split$N4, onlyObserved = T))
odd_na = data.frame(N1, N2, N4)
odd_na = t(odd_na)                                                                              ## Na per locus per population

## --------- Muller 2018 (2 populations) -------------------
setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SSR/Muller_2018_dataset")
mul18 <- read.loci("Muller_dataset_SSR.csv", header = T, loci.sep = ",", allele.sep = "/", col.pop = 2, col.loci = 3:11, na.strings = "1")
mul18_gin <- loci2genind(mul18, ploidy = 2, na.alleles = "1", unphase = T)
summary(mul18_gin)

## converting to genepop 
mul18_gpop = genind2genpop(mul18_gin, quiet = T, process.other = T)

mul18_gin$tab                                              ## N ind per loucs.allele per pop
mul18_allfr = as.data.frame(makefreq(mul18_gpop))                         ## allelic frequencies
  write.csv(mul18_allfr, "mul18_allfr.csv")
mul18_ar <- data.frame(allel.rich(mul18_gin,min.alleles = 25)$mean.richness)
mul18_he = data.frame(t(basic.stats(mul18_gin)$Hs))        # expected heterozygosity

mul_split = seppop(mul18_gin, mul18_gin$pop)
`DE-AEW5`= c(nAll(mul_split$`DE-AEW5`, onlyObserved = T))
`DE-HEW10` = c(nAll(mul_split$`DE-HEW10`, onlyObserved = T))
mul_na = data.frame(`DE-AEW5`,`DE-HEW10`)
mul_na = t(mul_na)                                        ## Na for each population

## --------- Cuervo 2018 ---------------
setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SSR/Cuervo_2018_dataset")
cue18 <- read.loci("Cuervo_SSR.csv", col.pop = 1, col.loci = 6:18, header = T, loci.sep = ",", allele.sep = "/", na.strings = "NA/NA")
cue18 = filter(cue18, Stage == "Adult")                    ## removing saplings from the dataset
cue18_gin <- loci2genind(cue18, ploidy = 2, na.alleles = "NA", unphase = T)

cue18_gpop = genind2genpop(cue18_gin, quiet = T, process.other = T)

cue18_gin$tab                                                  ## N ind per loucs.allele per pop
cue_allfr = as.data.frame(makefreq(cue18_gpop))                              ## allelic frequencies
  write.csv(cue_allfr, "cue18_allfr.csv")
cue18_he = data.frame(t(basic.stats(cue18_gin)$Hs))             # expected heterozygosity
cue18_ar <- data.frame(allel.rich(cue18_gin,min.alleles = 25)$mean.richness)   ## allelic richness (N=50)

cue_split = seppop(cue18_gin)
Ardon = c(nAll(cue_split$Ardon, onlyObserved = TRUE))
Chamoson = c(nAll(cue_split$Chamoson, onlyObserved = TRUE))
Chur = c(nAll(cue_split$Chur, onlyObserved = TRUE))
Collombey= c(nAll(cue_split$Collombey, onlyObserved = TRUE))
Felsberg= c(nAll(cue_split$Felsberg, onlyObserved = TRUE))
Malans= c(nAll(cue_split$Malans, onlyObserved = TRUE))
Martigny= c(nAll(cue_split$Martigny, onlyObserved = TRUE))
Mastrils= c(nAll(cue_split$Mastrils, onlyObserved = TRUE))
Mels= c(nAll(cue_split$Mels, onlyObserved = TRUE))
Ollon= c(nAll(cue_split$Ollon, onlyObserved = TRUE))
Sargans= c(nAll(cue_split$Sargans, onlyObserved = TRUE))
Saxon= c(nAll(cue_split$Saxon, onlyObserved = TRUE))

cue_na = data.frame(Ardon,Chamoson,Chur,Collombey,Felsberg,Malans,Martigny,Mastrils,Mels,Ollon,Sargans,Saxon)
cue_na = t(cue_na)                                             ##  Number of alleles per locus per population

## --------- Ulaszewski 2020 ---------
setwd("C:/Users/Camilla/Documents/UZH.internship/papers/Fagus sylvatica genetics and environmental stress/final.dataset/SSR/Ulaszewski_2020")
ula <- read.loci("Ula_SSR.csv", header = T, loci.sep = ",", allele.sep = "/", col.pop = 2, col.loci = 5:24, na.strings = "0/")
ula_gin <- loci2genind(ula, ploidy = 2, na.alleles = "0", unphase = T)

summary(ula_gin)
# create genepop file
ula_gpop = genind2genpop(ula_gin, quiet = T, process.other = T)  

ula_gpop$tab                                            ## N of individuals per locus.allele per population
ula_allfr = as.data.frame(makefreq(ula_gpop))                        ## allelic frequencies
write.csv(ula_allfr, "ula_allfr.csv")
ula_he = data.frame(t(basic.stats(ula_gin)$Hs))       ## Expected Heterozygosity for the populations
ula_ar <- data.frame(allel.rich(ula_gin, min.alleles = 25)$mean.richness)   ## Allelic richness rarefacted to standard N = 50

SSR.fst <- pairwise.fstb(ula_gin)  ## PopGenReport Pairwise Fst to compare with SNP-Fst

find.clusters(ula_gin)  ## interesting option to find clusters
inbreeding(ula_gin)   ## coefficient of inbreeding F to compare with SNP-F


## --------- distribution of allelic frequencies ------
## merging all.fr from all the studies
all_freq = full_join(lan11_allfr, lef12_allfr,  by = NULL, copy = T, keep = F)
all_freq = full_join(all_freq , pio12_allfr, by = NULL, copy = T, keep = F)
all_freq = full_join(all_freq , gau13_allfr, by = NULL, copy = T, keep = F)
all_freq = full_join(all_freq , del13_allfr, by = NULL, copy = T, keep = F)
all_freq = full_join(all_freq , raj14_af, by = NULL, copy = T, keep = F)
all_freq = full_join(all_freq , sjo14_gbaf, by = NULL, copy = T, keep = F)
all_freq = full_join(all_freq , sjo14_otaf, by = NULL, copy = T, keep = F)
all_freq = full_join(all_freq , cvr17_allfr, by = NULL, copy = T, keep = F)
all_freq = full_join(all_freq , odd18_allfreq, by = NULL, copy = T, keep = F)
all_freq = full_join(all_freq , mul18_allfr, by = NULL, copy = T, keep = F)
all_freq = full_join(all_freq , cue_allfr, by = NULL, copy = T, keep = F)
all_freq = full_join(all_freq , ula_allfr, by = NULL, copy = T, keep = F)
sort(colnames(all_freq))                  # rearrange in alphabetic order
setwd("C:/Users/Camilla/Documents/UZH.internship/sum_data/updated")
write.csv(all_freq, file = "allelic.freq_SSR.csv", na = "NA")

## plot of allelic frequencies
  all_f <- read.csv("C:/Users/Camilla/Documents/UZH.internship/sum_data/updated/allelic_freq_SSR.csv", sep = ",", header = T)
 # colnames(all_f) <- gsub(".1", "", colnames(all_f))
  all_f.long <- all_f %>% 
    tidyr::gather(locus, frequency, -X, -Reference, -Population, -Latitude, -Longitude) %>%
    separate(locus, into = c('locus', 'allele'), sep = "__") 
## fix allele problem!
  unique(all_f.long[which(all_f.long$locus == "DUKCT_A_0"),"allele"])
  unique(all_f.long[which(all_f.long$locus == "EEU75_A_0"),"allele"])
  unique(all_f.long[which(all_f.long$locus == "FS1_15"),"allele"])
  
  ## removing ".1" from duplicate allele columns
  all_f.long$allele <- gsub("87.1", "87", all_f.long$allele)
  all_f.long$allele <- gsub("93.1", "93", all_f.long$allele)
  all_f.long$allele <- gsub("85.1", "85", all_f.long$allele)
  all_f.long$allele <- gsub("83.1", "83", all_f.long$allele)
  all_f.long$allele <- gsub("95.1", "95", all_f.long$allele)
  all_f.long$allele <- gsub("89.1", "89", all_f.long$allele)
  all_f.long$allele <- gsub("75.1", "75", all_f.long$allele)
  all_f.long$allele <- gsub("97.1", "97", all_f.long$allele)
  all_f.long$allele <- gsub("91.1", "91", all_f.long$allele)
  all_f.long$allele <- gsub("99.1", "99", all_f.long$allele)
  all_f.long$allele <- gsub("94.1", "94", all_f.long$allele)
  all_f.long$allele <- gsub("96.1", "96", all_f.long$allele)
  all_f.long$allele <- gsub("99.1", "99", all_f.long$allele)
  all_f.long$allele <- gsub("98.1", "98", all_f.long$allele)
  
  
  no.cl <- length(unique(all_f$Reference))                          ## set number of colors = number of references
  myCol <- colorRampPalette(brewer.pal(8, "Set1"))(no.cl)         ## get the colors
  library("ggforce")
  FS2a <- ggplot(all_f.long, aes(x = as.factor(allele), y = frequency, na.rm = T, group = locus))+
          facet_wrap_paginate(~ locus, nrow = 4, ncol = 3, scales = "free_x", 
                              page = 1) +
          geom_jitter(alpha = 0.5,size = 1, na.rm = T, aes(color = factor(Reference)))+ 
          scale_color_manual("Reference", values = myCol)+
          theme_classic()+
          theme(strip.background = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size = 4.8), 
                axis.text.y = element_text(size = 5.5),
                strip.text.y = element_blank(),
                axis.title = element_text(hjust = 0.5,face = "bold", size = 8), 
                panel.border = element_rect(colour = "black", fill = NA, size = 1), 
                legend.position = "top",
                legend.title = element_text(face = "bold"),
                legend.text = element_text(size = 8),
                text = element_text(size = 9))+  
    guides(color = guide_legend(override.aes = list(size = 7,shape = 15, alpha = 1), nrow = 4, ncol = 3))+
          labs(x = "Allele", y = "Allelic frequency")
  FS2a
   
  setwd("C:/Users/Camilla/Dropbox/FagusGenDiv/manuscript/Figures")
  jpeg(filename = "FigS2a.jpg", width = 174, height = 230, res = 700, units = "mm", quality =  100) 
  plot(FS2a)
  dev.off()
  
 FS2b <- ggplot(all_f.long, 
         aes(x = as.factor(allele), y = frequency, na.rm = T, color = Reference, group = locus))+
   facet_wrap_paginate(~ locus, nrow = 5, ncol = 3, scales = "free_x", page = 2) +
     geom_jitter(alpha = 0.5,size = 1, na.rm = T, aes(color = factor(Reference)))+ 
   scale_color_manual(values = myCol)+
   theme_classic()+
   theme(strip.background = element_blank(),
         #plot.title = element_text(vjust =  -10, size = 10, face = "bold"), 
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size = 4.6), 
         axis.text.y = element_text(size = 7),
         strip.text.y = element_blank(),
         axis.title = element_text(hjust = 0.5,face = "bold", size = 9), 
         panel.border = element_rect(colour = "black", fill = NA, size = 1), 
         legend.position = "top",
         legend.title = element_text(face = "bold"),
         legend.text = element_text(size = 8),
         text = element_text(size = 9))+  
   guides(color = FALSE)+
   labs(x = "Allele", y = "Allelic frequency")
 FS2b
 
 ## add Pag 2
 setwd("C:/Users/Camilla/Dropbox/FagusGenDiv/manuscript/Figures")
 jpeg(filename = "FigS2b.jpg", width = 174, height = 230, res = 500, units = "mm", quality =  100) 
 plot(FS2b)
 dev.off()
 
 FS2c <- ggplot(all_f.long, 
                aes(x = as.factor(allele), y = frequency, na.rm = T, color = Reference, group = locus))+
   facet_wrap_paginate(~ locus, nrow = 2, ncol = 3, scales = "free_x", page = 3) +
   geom_jitter(alpha = 0.5,size = 1, na.rm = T, aes(color = factor(Reference)))+ 
   scale_color_manual(values = myCol)+
   theme_classic()+
   theme(strip.background = element_blank(),
         #plot.title = element_text(vjust =  -10, size = 10, face = "bold"), 
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size = 4.8), 
         axis.text.y = element_text(size = 7),
         strip.text.y = element_blank(),
         axis.title = element_text(hjust = 0.5,face = "bold", size = 9), 
         panel.border = element_rect(colour = "black", fill = NA, size = 1), 
         legend.position = "top",
         legend.title = element_text(face = "bold"),
         legend.text = element_text(size = 8),
         text = element_text(size = 9))+  
   guides(color = FALSE)+
   labs(x = "Allele", y = "Allelic frequency")
 FS2c
 
 ## add Pag 3
 setwd("C:/Users/Camilla/Dropbox/FagusGenDiv/manuscript/Figures")
 jpeg(filename = "FigS2c.jpg", width = 174, height = 90, res = 500, units = "mm", quality =  100) 
 plot(FS2c)
 dev.off()
 
## --------- analysis by groups of references (Set of Loci) -------
 ## common loci between the studies
 bin_loci = read.csv("C:/Users/Camilla/Dropbox/FagusGenDiv/Data/sum_data/sum_data_02.2021/bin_loci.csv", sep= ",", header = T)  ## binary dataframe: 0 = the locus was not used, 1 = the locus was used
 bin_loci[is.na(bin_loci)] = 0
 ncol(bin_loci)
 row.names(bin_loci) = bin_loci[,1] ## set the references as row names
 bin_loci = bin_loci[, -1]
 mat_bin = data.matrix(bin_loci, rownames.force = NA)   #converting df to matrix
 
 ## Fig S2 (heatmap)
 setwd("C:/Users/Camilla/Dropbox/FagusGenDiv/manuscript/Figures")

 hclustfunc <- function(x) hclust(x, method="complete")
 distfunc <- function(x) dist(x,method="euclidean")
 d <- distfunc(mat_bin)
 fit <- hclustfunc(d)
 
 #plot dendogram only with groups of studies
 plot(fit)
 groups <- cutree(fit, k=5) 
 rect.hclust(fit, k=5, border="red")    # Add rectangle in cluster
 library("gplots")
 jpeg(filename = "FigS2.jpg", width = 174, height = 100, res = 600, units = "mm", quality =  100) 
 heatmap.2(as.matrix(mat_bin),dendrogram="none",trace="none", margin=c(8,9), 
           hclust=hclustfunc, distfun=distfunc, RowSideColors=as.character(groups),
           key = FALSE,col=c("#332288", "#DDCC77"), sepcolor="white",
           colsep = seq(1, 60),rowsep = seq(1, 60),
           adjRow = c(0.025,0.5),adjCol = c(0.8,0.5),
           sepwidth=c(0.05,0.05), margins = c(6,10),
           lmat=rbind(c(5,0,4,0,0), c(4,1,2,3,0)),lhei=c(0.3,18), 
           lwid = c(0.2,0.5,18,0.25,0.25))
 dev.off()

 ## gathering studies according to the set of loci used (creating 5 groups)
 SSR$LociSet <- ifelse(SSR$reference == "Gauzere et al., 2013" | SSR$reference == "Oddou-Muratorio et al., 2018"| 
                         SSR$reference == "Kempf et al., 2016"| SSR$reference == "Bontemps et al., 2016"| SSR$reference == "Lander et al., 2011", "1",
                       ## group 2
                       ifelse(SSR$reference == "Muller et al., 2018"| SSR$reference == "Rajendra et al., 2014"| 
                                SSR$reference == "Cuervo-Alarcon et al., 2018" , "2",
                        ## group 3      
                              ifelse(SSR$reference == "Sjolund et al., 2015"| SSR$reference == "Chybicky et al., 2009"| 
                                       SSR$reference == "Pastorelli et al., 2002"| SSR$reference == "Kraj et al., 2009"| 
                                       SSR$reference == "Bilela et al., 2012"| SSR$reference == "Paffetti et al., 2012"| 
                                       SSR$reference == "Buiteveld et al., 2007"| SSR$reference == "Piotti et al., 2012"| 
                                       SSR$reference == "Nyári 2010"| SSR$reference == "Nowakowska 2011"| 
                                       SSR$reference == "Cvrckova et al., 2017"| SSR$reference == "Vornam et al., 2004", "3",
                            # group 4         
                                     ifelse(SSR$reference == "Lefevre et al., 2012"|SSR$reference == "De Lafontaine et al., 2013"|
                                              SSR$reference == "Sandurska et al., 2017"|SSR$reference == "Ulaszewski et al., 2021", "4",
                                 ## group 5           
                                            ifelse(SSR$reference == "Pluess et al., 2016"|SSR$reference == "Pluess et al., 2013", "5", "something else")))))
 
 ## ANOVA test between the groups
 kruskal.test(He ~ LociSet, SSR)
 kruskal.test(Ar ~ LociSet, SSR)
 ##post-hoc test
 install.packages("dunn.test")
 library(dunn.test)
 #library(FSA)
 dunnTest(SSR$He, SSR$LociSet, method = "bonferroni")
 
 
## --------- scaling He per locus per study -----------
all_freq <- read.csv("C:/Users/Camilla/Documents/UZH.internship/sum_data/updated/all_freq_SSR11.csv")
  summary(all_freq) 

  gen <- all_freq[,-c(1:5)]
  all.loci <- unique(unlist(strsplit(names(gen), ".", fixed=T))[seq(1,2*ncol(gen),2)])
  Hall <- matrix(ncol=length(all.loci), nrow=nrow(gen))  ## empty matrix for He scaled
  Ho.all <- matrix(ncol=length(all.loci), nrow=nrow(gen)) ## empty matrix for He observed
  
  ## scaling across all studies and all populations
  for(loc in 1:length(all.loci)) {
    myloc <- all.loci[loc]                              # take the locus 
    mydat <- all_freq[, grep(myloc, names(all_freq))]   # take the allelic frequencies for that locus 
    num.alleles <- ncol(mydat)                          # take the number of alleles for that locus 
    M <- apply(mydat, 1, max, na.rm=T)                  # frequency of the most frequent allele for each locus 
    K <- unlist(apply(mydat, 1, function(a) length(na.omit(a)))) # No alleles for each locus 
    Hobs <- apply(mydat, 1, function(a)
      ifelse(sum(is.na(a)) == num.alleles, NA, sum(a^2, na.rm=T)))   # Homozygosity for each locus and each population 
    Hmin <- (K*M^2 - 2*M + 1)/(K-1)                     # lower theoretical boundary
      # ceiling function: https://en.wikipedia.org/wiki/Floor_and_ceiling_functions 
    Hmax <- (1-M*(ceiling(1/M)-1)*(2-ceiling(M^-1)*M))  # upper theoretical boundary
    Hrange <- Hmax-Hmin
    Ho.all[, loc] <- (1-Hobs)                           # matrix for observed He
    Hall[,loc] <- (Hobs-Hmin)/(Hmax-Hmin)               # scaling homozygosity (per locus per population)
  }
  
  hist(Hall) 
  hist(Hmax)
  hist(Hmin)
  hist(1-Hall) ## scaled He
  hist(Ho.all) ## raw He

## Plot of the scaled He
  He_sc <- as.data.frame(1-Hall)
  ## add information about Population, reference and coordinates
    colnames(He_sc) <- all.loci
    for(i in 1:ncol(He_sc)){
      He_sc[,i] <- as.numeric(He_sc[,i])
    }
    He_sc <- He_sc %>%
      mutate("scaled_he_mean" = rowMeans(He_sc, na.rm = T))
    summary(He_sc$scaled_He_mean)
    He_sc <- He_sc %>% mutate("Population" = all_freq$Population)
    He_sc <- He_sc %>% mutate("Reference" = all_freq$Reference)
    He_sc <- He_sc %>% mutate("Latitude" = all_freq$Latitude)
    He_sc <- He_sc %>% mutate("Longitude" = all_freq$Longitude)
    
He_sc.long <- He_sc %>% gather(variable, value, -Population, -N, -Reference, -Latitude, -Longitude, -scaled_he_mean)
 str(He_sc.long)
ggplot(He_sc.long, aes(x = variable, y = value, na.rm = T))+
    geom_boxplot(size=0.7, na.rm=T)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size = 7), 
          axis.title = element_text(hjust = 0.5,face = "bold", size = 9), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          legend.position = "none")+
    scale_y_continuous(limits = c(0,1),expand=c(0,0))+
    labs(x = "Locus", y = "scaled He", tag = "b")
  # write.csv(He_sc, "SSR_He_scaled.csv", sep = ",")

  He_obs <- as.data.frame(1-Hobs)
   He_obs <- as.data.frame(Ho.all)
   ## add information about Population, reference and coordinates
   colnames(He_obs) <- all.loci
   for(i in 1:ncol(He_obs)){
     He_obs[,i] <- as.numeric(He_obs[,i])
   }
   He_obs <- He_obs %>% mutate("he_mean" = rowMeans(He_obs, na.rm = T))
   He_obs <- He_obs %>% mutate("Population" = all_freq$Population)
   He_obs <- He_obs %>% mutate("Reference" = all_freq$Reference)
   He_obs <- He_obs %>% mutate("Latitude" = all_freq$Latitude)
   He_obs <- He_obs %>% mutate("Longitude" = all_freq$Longitude)
  # write.csv(He_obs, "SSR_He_obs.csv", sep = ",")
   
   He_obs.long <- He_obs %>% gather(variable, value, -Population, -Reference, -Latitude, -Longitude, -he_mean)
   ggplot(He_obs.long, aes(x = variable, y = value, na.rm = T))+
     geom_boxplot(size=0.7, na.rm=T)+
     theme_bw()+
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size = 7), 
           axis.title = element_text(hjust = 0.5,face = "bold", size = 9), 
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
           legend.position = "none")+
     scale_y_continuous(limits = c(0,1),expand=c(0,0))+
     labs(x = "Locus", y = "raw He", tag = "a")

   ## merging dataframes and plot
  merge <- merge( He_obs.long, He_sc.long, by = c("Reference", "Population", "variable", "Latitude", "Longitude"))
  colnames(merge)[colnames(merge) == "variable"] <- "locus"
  merge.long <-  merge %>% gather(variable, value, -Population, -Reference, -locus, -Latitude, -Longitude, -he_mean, -scaled_he_mean, -N)
  colnames(merge.long)[colnames(merge.long) == "value"] <- "He"
  FigS3 <- ggplot(merge.long, aes(x = locus, y = He, fill = variable, na.rm = T))+
    geom_boxplot(size=0.5, na.rm=T, outlier.size = 1, notch = T)+
    scale_fill_manual(name="", values = c("#332287", "#DDCC76"),labels = c("raw He", "scaled He"))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size = 7), 
          axis.title = element_text(hjust = 0.5,face = "bold", size = 9), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          legend.position = c(0.9, 0.15))+
    scale_y_continuous(limits = c(0,1),expand=c(0,0))+
    labs(x = "Locus")
  FigS3
    
  setwd("C:/Users/Camilla/Dropbox/FagusGenDiv/manuscript/Figures")
  jpeg(filename = "FigS3.jpg", width = 174, height = 100, res = 600, units = "mm", quality =  100) 
  plot(FigS3)
  dev.off()
  
## --------- repeating the analysis by group of study ------
  ## mean test between the loci set groups (raw He) 
  # assign one group number to each study
  He_obs$LociSet <- ifelse(He_obs$Reference == "Gauzere et al., 2013"| He_obs$Reference == "Oddou-Muratorio et al., 2018"| He_obs$Reference == "Lander et al., 2011", "1",
                           ifelse(He_obs$Reference == "Muller et al., 2018"|He_obs$Reference == "Rajendra et al., 2014"| He_obs$Reference == "Cuervo-Alarcon et al., 2018" , "2",
                                  ifelse(He_obs$Reference == "Sjolund et al., 2014 gb "| He_obs$Reference == "Sjolund et al., 2014"|He_obs$Reference == "Piotti et al., 2012"| He_obs$Reference == "Cvrckova et al., 2017", "3",
                                         ifelse(He_obs$Reference == "Lefevre et al., 2012"|He_obs$Reference == "De Lafontaine et al., 2013"|He_obs$Reference == "Ulaszewski et al., 2021", "4", "something else"))))
  
  ggplot(He_obs, aes(x = LociSet, y = he_mean, group = LociSet, fill = LociSet)) +  
    geom_boxplot() +
    geom_jitter(color="black", size = 1, alpha = 1)+
    labs(x= "Set of loci")
  
  shapiro.test(He_obs$he_mean)
  kruskal.test(he_mean ~ LociSet, He_obs)
  out <- boxplot.stats(He_obs.long$value)$out  #outliers vector
  pop_out <- which(He_obs.long$value %in% c(out)) # outliers rows
  out1 <- he_loc[pop_out, ]
  out1_loci <- count(group_by(out1, Locus)) # outlier data per locus
  
  
  ## mean test between the loci set groups (scaled He) 
  # assign one group number to each study
  He_sc$LociSet <- ifelse(He_sc$Reference == "Gauzere et al., 2013"| He_sc$Reference == "Oddou-Muratorio et al., 2018"| He_sc$Reference == "Lander et al., 2011", "1",
                        ifelse(He_sc$Reference == "Muller et al., 2018"|He_sc$Reference == "Rajendra et al., 2014"| He_sc$Reference == "Cuervo-Alarcon et al., 2018" , "2",
                               ifelse(He_sc$Reference == "Sjolund et al., 2014 gb "| He_sc$Reference == "Sjolund et al., 2014"|He_sc$Reference == "Piotti et al., 2012"| He_sc$Reference == "Cvrckova et al., 2017", "3",
                                      ifelse(He_sc$Reference == "Lefevre et al., 2012"|He_sc$Reference == "De Lafontaine et al., 2013"|He_sc$Reference == "Ulaszewski et al., 2021", "4", "something else"))))
 
  ggplot(He_sc, aes(x = LociSet, y = scaled_he_mean, group = LociSet, fill = LociSet)) +  
    geom_boxplot() +
    geom_jitter(color="black", size = 1, alpha = 1)+
    labs(x= "Set of loci")
  
  ## mean test between the loci set groups (scaled He)
  shapiro.test(He_sc$scaled_he_mean)
  kruskal.test(scaled_he_mean ~ LociSet, He_sc)
  ##post-hoc test
  install.packages("dunn.test")
  library(dunn.test)
  library(FSA)
  dunnTest(He_sc$scaled_he_mean, He_sc$LociSet, method = "bonferroni")
  
  colMeans(He_obs[, 1:42], na.rm = TRUE)
  
## --------- Correlation He - geographic variables ------
  setwd("C:/Users/Camilla/Documents/UZH.internship/sum_data/updated")
  
  ## original data
  He_obs <- read.csv("SSR_He_obs.csv", sep = ",")
  hist(He_obs$he_mean)
  sd(He_obs$he_mean)
  shapiro.test(He_obs$he_mean)
 
  # Correlation He-latitude (spearman)
      cor.test(He_obs$he_mean, He_obs$Latitude, alternative = "two.sided", 
               method = c("spearman"), 
               exact = NULL, conf.level = 0.95, continuity = F)
      # only Ulaszewski dataset
      Ula_obs <- He_obs[which(He_obs$Reference == "Ulaszewski and Burczyk, 2021"),]
      cor.test(Ula_obs$he_mean, Ula_obs$Latitude, alternative = "two.sided", 
               method = c("spearman"), exact = NULL, conf.level = 0.95, continuity = F)
  
      no.cl <- length(unique(He_obs$Reference))                       ## set number of colors = number of references
      myCol <- colorRampPalette(brewer.pal(8, "Set1"))(no.cl)         ## get the colors
      
  f4a <- ggplot(He_obs, aes(x=Latitude, y=he_mean))+
        geom_point(data = He_obs, size = 2.5, alpha=0.8,aes(color = factor(Reference)))+ 
        scale_color_manual(values = myCol)+
        geom_smooth(method = lm, color="black")+
        scale_y_continuous(breaks = seq(0.4, 0.8, by = 0.1),limits = c(0.4, 0.8),expand=c(0,0))+
        theme_classic()+
        theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 9, margin = margin(0, 0, 0.1, 0, "cm")), 
              axis.text.y = element_text(vjust = 0.5, hjust=0.5, size = 9, margin = margin(0, 0, 0.1, 0, "cm")), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.title = element_text(hjust = 0.5,face = "bold", size = 9), 
              panel.border = element_rect(colour = "black", fill = NA, size = 1), 
              legend.position = "right",
              legend.title = element_text(face = "bold"),
              legend.text = element_text(size = 8),
              text = element_text(size = 9))+  
        guides(color = guide_legend(override.aes = list(size = 7,shape = 15, alpha = 1), nrow = 4, ncol = 3))+
        annotate(geom = "text", label ="r = -0.31, p = 1.46e-06", x = 44, y = 0.42, size = 3)+
        labs(y ="Raw He", x = "Latitude", color = "Reference")
      f4a
      
    # Correlation He-longitude (spearman)
      cor.test(He_obs$he_mean, He_obs$Longitude, alternative = "two.sided", 
               method = c("spearman"), 
               exact = NULL,conf.level = 0.95,continuity = F)
      # only Ulaszewski
      cor.test(Ula_obs$he_mean, Ula_obs$Longitude, alternative = "two.sided", method = c("spearman"), exact = NULL, conf.level = 0.95, continuity = F)
    
      f4c <- ggplot(He_obs, aes(y=he_mean, x=Longitude))+
        geom_point(data = He_obs, size = 2.5, alpha=0.8,aes(color = factor(Reference)))+ 
        scale_color_manual(values = myCol)+
        geom_smooth(method = lm, color="black")+
        scale_y_continuous(breaks = seq(0.4, 0.8, by = 0.1),limits = c(0.4, 0.8),expand=c(0,0))+
        theme_classic()+
        theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 9, margin = margin(0, 0, 0.1, 0, "cm"), colour = "black"), 
              axis.text.y = element_text(vjust = 0.5, hjust=0.5, size = 9, margin = margin(0, 0, 0.1, 0, "cm"), colour = "black"), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.title = element_text(hjust = 0.5,face = "bold", size = 9), 
              panel.border = element_rect(colour = "black", fill = NA, size = 1))+
        guides(FALSE)+
        annotate(geom = "text", label ="r = -0.19, p = 0.004", x = 4, y = 0.42, size = 3)+
        labs(y ="Raw He", x = "Longitude", color = "Reference")
f4c
    
## correlation He - distance from refugia
      ## adding the coordinates of the core region
      He_obs$Lat_center1 <- 46.02
      He_obs$Lon_center1 <- 14.65
      
      #distance in km of each population from this point
      library(geosphere)
      mat1 <- distm(He_obs[, c("Longitude", "Latitude")], He_obs[,c("Lon_center1", "Lat_center1")], fun=distVincentyEllipsoid)
      He_obs$dist1 <- apply(mat1,1, min)
      He_obs$dist1 <- He_obs$dist1 / 1000

      ## correlation He - distance from the center (assumed to be in slovenia) 
      cor.test(He_obs$dist1,He_obs$he_mean,alternative = "two.sided", 
               method = c("spearman"),exact = NULL, 
               conf.level = 0.95, continuity = F)
     ## Ulaszewski only
      Ula_obs <- He_obs[which(He_obs$Reference == "Ulaszewski and Burczyk, 2021"),]
      cor.test(Ula_obs$dist1, Ula_obs$he_mean, alternative = "two.sided", 
               method = c("spearman"),exact = NULL, conf.level = 0.95, continuity = F)
      
      spain <- subset(He_obs, Latitude <= 43.50 & Longitude <= 3.2 & Longitude >= -7) # only Spanish populations
      french <- subset(He_obs, Latitude <= 45.8 & Latitude >= 43.2 & Longitude >= 4.86 & Longitude <= 6.61) # only french populations
      
    f4e <-  ggplot()+
        theme_bw()+
        geom_point(data = He_obs, aes(x=dist1, y=he_mean), size =2.5, alpha=0.8)+
        geom_smooth(data = He_obs, aes(x=dist1, y=he_mean), method = "lm", color="black")+
      geom_point(data = spain, aes(x=dist1, y=he_mean), colour = "#FF9999", size = 2.6, alpha=0.9)+            
      geom_point(data = french, aes(x=dist1, y=he_mean), colour = "#56B4E9", size = 2.6, alpha=0.9)+       
      scale_y_continuous(breaks = seq(0.4, 0.8, by = 0.1),limits = c(0.4, 0.8),expand=c(0,0))+
      theme_classic()+
      theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 9, margin = margin(0, 0, 0.1, 0, "cm"), colour = "black"), 
            axis.text.y = element_text(vjust = 0.5, hjust=0.5, size = 9, margin = margin(0, 0, 0.1, 0, "cm"), colour = "black"), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_text(hjust = 0.5,face = "bold", size = 9), 
            panel.border = element_rect(colour = "black", fill = NA, size = 1))+
      annotate(geom = "text", label ="r = 0.26, p = 5.82e-15", x = 500, y = 0.42, size = 3)+
      labs(x ="Distance from Origin (km)", y = "Raw He", size = 9)   
  f4e
      ## scaled data
      He_sc <- read.csv("SSR_He_scaled.csv", sep = ",")
      hist(He_sc$scaled_he_mean)
      sd(He_sc$scaled_he_mean)
      shapiro.test(He_sc$scaled_he_mean)
      str(He_sc)

  ## correlation He - latitude (Spearman) 
  cor.test(He_sc$Latitude, He_sc$scaled_he_mean,          
           alternative = "two.sided", 
           method = c("spearman"),
           exact = NULL, conf.level = 0.95, continuity = F)
   ## Ulaszewski only
  Ula_sc <- He_sc[which(He_sc$Reference == "Ulaszewski and Burczyk, 2021"),]
  cor.test(Ula_sc$Latitude, Ula_sc$scaled_he_mean,          
           alternative = "two.sided", 
           method = c("spearman"),
           exact = NULL, conf.level = 0.95, continuity = F)
  
  f4b <- ggplot(He_sc, aes(x=Latitude, y=scaled_he_mean))+
    geom_point(data = He_sc, size = 2.5, alpha=0.8,aes(color = factor(Reference)))+ 
    scale_color_manual(values = myCol)+
    geom_smooth(method = lm, color="black")+
    scale_y_continuous(breaks = seq(0.4, 0.8, by = 0.1),limits = c(0.4, 0.8),expand=c(0,0))+
    theme_classic()+
    theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 9, margin = margin(0, 0, 0.1, 0, "cm"), colour = "black"), 
          axis.text.y = element_text(vjust = 0.5, hjust=0.5, size = 9, margin = margin(0, 0, 0.1, 0, "cm"), colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(hjust = 0.5,face = "bold", size = 9), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1))+
    guides(FALSE)+
    annotate(geom = "text", label ="r = 0.22, p = 0.001", x = 43, y = 0.42, size = 3)+
    labs(y ="Scaled He", x = "Latitude", color = "Reference")
 f4b 

   ## correlation He - longitude 
  cor.test(He_sc$Longitude,He_sc$scaled_he_mean,          
           alternative = "two.sided", method = c("spearman"),
           exact = NULL, conf.level = 0.95, continuity = F)
  ##only Ulaszewski
  cor.test(Ula_sc$Longitude, Ula_sc$scaled_he_mean,          
           alternative = "two.sided", 
           method = c("spearman"),
           exact = NULL, conf.level = 0.95, continuity = F)
  
  f4d <- ggplot(He_sc, aes(y=scaled_he_mean, x=Longitude))+
    geom_point(data = He_sc, size = 2.5, alpha=0.8,aes(color = factor(Reference)))+ 
    scale_color_manual(values = myCol)+
    geom_smooth(method = lm, color="black")+
    scale_y_continuous(breaks = seq(0.4, 0.8, by = 0.1),limits = c(0.4, 0.8),expand=c(0,0))+
    theme_classic()+
    theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 9, margin = margin(0, 0, 0.1, 0, "cm"), colour = "black"), 
          axis.text.y = element_text(vjust = 0.5, hjust=0.5, size = 9, margin = margin(0, 0, 0.1, 0, "cm"), colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(hjust = 0.5,face = "bold", size = 9), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1))+
    guides(color = FALSE)+
    annotate(geom = "text", label ="r = 0.4, p = 2.4e-10", x = 4, y = 0.42, size = 3)+
    labs(y ="Scaled He", x = "Longitude", color = "Reference")
  f4d
  
## correlation He - distance from refugia
  ## adding the coordinates of the core region 1: Slovenia
  He_sc$Lat_center1 <- 46.02
  He_sc$Lon_center1 <- 14.65
  
  #distance in km of each population from this point
  library(geosphere)
  mat2 <- distm(He_sc[, c("Longitude", "Latitude")], He_sc[,c("Lon_center1", "Lat_center1")], fun=distVincentyEllipsoid)
  He_sc$dist1 <- apply(mat2,1, min)
  He_sc$dist1 <- He_sc$dist1 / 1000
  hist(He_sc$dist1)
 
  ## correlation between He and distance from the center (assumed to be in slovenia) 
  cor.test(He_sc$dist1,He_sc$scaled_he_mean,alternative = "two.sided", 
           method = c("spearman"),exact = NULL, conf.level = 0.95, continuity = F)
  # only Ulaszewski
  Ula_sc <- He_sc[which(He_sc$Reference == "Ulaszewski and Burczyk, 2021"),]
  cor.test(Ula_sc$dist1,Ula_sc$scaled_he_mean,alternative = "two.sided", 
           method = c("spearman"),exact = NULL, conf.level = 0.95, continuity = F)
  
  spain1 <- subset(He_sc, Latitude <= 43.50 & Longitude <= 3.2 & Longitude >= -7) # subsetting Spanish populations
  french1 <- subset(He_sc, Latitude <= 45.8 & Latitude >= 43.2 & Longitude >= 4.86 & Longitude <= 6.61)
  str(He_sc)
  
  f4f <- ggplot()+
    theme_bw()+
    geom_point(data = He_sc, aes(x=dist1, y=scaled_he_mean), size =2.5, alpha=0.8)+
    geom_smooth(data = He_sc, aes(x=dist1, y=scaled_he_mean), method = "lm", color="black")+
    geom_point(data = spain1, aes(x=dist1, y=scaled_he_mean),colour = "#FF9999", size = 2.6, alpha=0.9)+                #shape = 17?
    geom_point(data = french1, aes(x=dist1, y=scaled_he_mean),colour = "#56B4E9",size = 2.6, alpha=0.9)+              # shape = 15?
    scale_y_continuous(breaks = seq(0.4, 0.8, by = 0.1),limits = c(0.4, 0.8),expand=c(0,0))+
    theme_classic()+
    theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 9, margin = margin(0, 0, 0.1, 0, "cm"), colour = "black"), 
          axis.text.y = element_text(vjust = 0.5, hjust=0.5, size = 9, margin = margin(0, 0, 0.1, 0, "cm"), colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(hjust = 0.5,face = "bold", size = 9), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1))+
    guides(col = F)+
    annotate(geom = "text", label ="r = -0.43, p = 1.60e-11", x = 500, y = 0.42, size = 3)+
    labs(x ="Distance from Origin (km)", y = "Scaled He", size = 9) 
f4f 

F4 <- ggarrange(f4a, f4b, f4c, f4d, f4e, f4f,
           ncol = 2, nrow = 3,
           labels = c("a", "b", "c", "d", "e", "f"),
           font.label = list(size = 15, face = "bold"),
           common.legend = T, 
           align = "hv",
           legend = "top") 
F4
setwd("C:/Users/Camilla/Dropbox/FagusGenDiv/manuscript/Figures")
 jpeg(filename = "Fig4.jpg", width = 174, height = 200, res = 700, units = "mm", quality =  100) 
 plot(F4)
 dev.off()
 
## --------- He maps--------
 setwd("C:/Users/Camilla/Documents/UZH.internship/sum_data/updated")
 map.world<-map_data(map="world")
 fagusrange <- sf::st_read("C:/Users/Camilla/Documents/UZH.internship/Occurrences/occurrences/Fagus sylvatica.shp")
 
 fagusrange <- sf::st_read("C:/Users/Camilla/Documents/UZH.internship/Fagus_map/Fagus_sylvatica/Fagus_sylvatica_EUFORGEN.shp")
 myCols <- brewer.pal(4, "Set1") # Set the color palette (RColorBrewer library)

 ## OBSERVED HE
 He_obs <- read.csv("SSR_He_obs.csv", sep = ",")

 ## Figure 3A: observed He from SSR-based studies 
 F3a <- ggplot() + 
   geom_map(data=map.world, map=map.world, aes(map_id=region), fill="grey40", alpha = 0.3, colour="white", size=0.15)+
   theme_void()+
   geom_sf(fagusrange, mapping =aes(fill= "Fagus sylvatica range"), color = NA, fill = "grey40", alpha = 0.7)+
   geom_point(data=He_obs, aes(x=Longitude, y=Latitude, fill = he_mean), size=1.5, shape = 21, alpha = 0.8, color = "white", stroke = 0.2)+
   coord_sf(ylim=c(37,60), xlim=c(-10,28))+
   scale_fill_viridis_b(option = "viridis",direction = -1,
                        limits = c(0.50, 0.9), breaks = c(0.5, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90))+ 
   guides(fill = guide_colorbar(direction = "horizontal",label.position = "bottom",
                                title.position="top",title = "Expected heterozygosity"))+
   theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         legend.key.height = unit(0.3, 'cm'), 
         legend.key.width = unit(1.5, 'cm'), 
         legend.title = element_text(size=8, face = "bold"), 
         legend.text = element_text(size= 7, face = "bold"), 
         legend.position = "bottom")+
   labs(x = " ", y = " ", tag = "a")+
   ggsn::scalebar(fagusrange, dist = 400, st.size=2, border.size = 0.2, height=0.02, transform = TRUE, model = 'WGS84', dist_unit = "km", location = "bottomleft", size = 0.5)
 F3a
 
 ### SCALED HE
 He_sc <- read.csv("SSR_He_scaled.csv", sep = ",")
      ## removing Sjolund populations in GB out of the Fagus range and plot on the map
     #He_sc <- He_sc[!(He_sc$Reference == "Sjolund 2014 gb "& He_sc$Latitude>=53 | He_sc$Longitude < -3),]

 ## Figure 3B: scaled He from SSR-based studies 
F3b<- ggplot() + 
   geom_map(data=map.world, map=map.world, aes(map_id=region), fill="grey40", alpha = 0.3, colour="white", size=0.15)+
   theme_void()+
   geom_sf(fagusrange, mapping =aes(fill= "Fagus sylvatica range"), color = NA, fill = "grey40", alpha = 0.7)+
   geom_point(data=He_sc, aes(x=Longitude, y=Latitude, fill = scaled_he_mean), size=1.5, shape = 21, alpha = 0.8, color = "white", stroke = 0.2)+
  coord_sf(ylim=c(37,60), xlim=c(-10,28))+
  scale_fill_viridis_b(option = "viridis",direction = -1)+ 
  guides(fill = FALSE)+
  theme(#panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())+
  labs(x = " ", y = " ", tag = "b")
plot(F3b)

## --------- IDW -----------
setwd("C:/Users/Camilla/Documents/UZH.internship/sum_data/updated")
he <- read.csv("SSR_He_scaled.csv", sep = ",")
## trasform into sf object with WGS84 coordinate system
x <- SpatialPointsDataFrame(cbind(he$Longitude, he$Latitude), data.frame(he=he$scaled_he_mean), 
                            proj4string= CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Create an empty grid where n is the total number of cells
newdat              <- as.data.frame(spsample(x, "regular", n=50000))
names(newdat)       <- c("Longitude", "Latitude")
coordinates(newdat) <- c("Longitude", "Latitude")
gridded(newdat)     <- TRUE  # Create SpatialPixel object
fullgrid(newdat)    <- TRUE  # Create SpatialGrid object
proj4string(newdat) <- proj4string(x)

myidw <- gstat::idw(he ~ 1, x, newdata=newdat, idp=5)  # IDW with power value 5
library("raster")
r <- raster(myidw)

#interpolation map for ggplot
map.world<-map_data(map="world")
f <- rgdal::readOGR("C:/Users/Camilla/Documents/UZH.internship/Fagus_map/Fagus_sylvatica/Fagus_sylvatica_EUFORGEN.shp")
            
# mask the map to the Fagus sylvatica range
r2 <- crop(r, extent(f))
r3 <- mask(r2, f)
plot(r3)
plot(f, add=TRUE, lwd=2)
r3_df <- as.data.frame(r3, xy = TRUE)

F3c <- ggplot() +
  theme_void()+
  geom_map(data=map.world, map=map.world, aes(map_id=region), fill="grey40", alpha = 0.3, colour="white", size=0.3)+
  geom_raster(r3_df, mapping = aes(x = x, y = y, fill = var1.pred), interpolate = T)+
  scale_fill_viridis_c(option = "viridis", na.value = "transparent", direction = -1)+
  coord_sf(ylim=c(37,60), xlim=c(-10,28))+
  guides(fill = FALSE)+
  theme(
    axis.text.x = element_blank(),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.key.height = unit(0.5, 'cm'), 
    legend.key.width = unit(2, 'cm'), 
    legend.spacing = unit(0.01, "cm"),
    legend.spacing.y = unit(0.05, "cm"), 
    legend.title = element_text(size=15, face = "bold"), 
    legend.text = element_text(size= 13, face = "bold"))+ 
    labs(x = " ", y = " ", tag = "c")
F3c

colnames(map.world)[match(c("long","lat"),colnames(map.world))] = c("Longitude","Latitude")

#merging panels
fig3 <- ggarrange(F3a, F3b, F3c,
                  #labels = c("a", "b", "c"),
                  ncol = 3, nrow = 1,
                  align = "hv",
                  widths = c(2,2,2),
                  heights = c(2,2,2),
                  common.legend = T, 
                  legend = "bottom")
fig3
setwd("C:/Users/Camilla/Dropbox/FagusGenDiv/manuscript/Figures")
jpeg(filename = "Fig3.jpg", width = 174, height = 60, res = 600, units = "mm", quality =  100) 
plot(fig3)
dev.off()

## --------- IDW Cross Validation ---------
# Leave-one-out validation routine
IDW.out <- vector(length = length(x))
for (i in 1:length(x)) {
  IDW.out[i] <- idw(he ~ 1, x[-i,], x[i,], idp=5.0)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
setwd("C:/Users/Camilla/Dropbox/FagusGenDiv/manuscript/Figures")
jpeg(filename = "FigS4a.jpg", width = 174, height = 100, res = 400, units = "mm", quality =  100) 
plot(IDW.out ~ x$he, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ x$he), col="red", lw=2,lty=2)
abline(0,1)
par(OP)

dev.off()

# Compute RMSE
sqrt( sum((IDW.out - x$he)^2) / length(x))

## 95% CI map from IDW with power parameter = 5
# jackknife technique to estimate a confidence interval at each unsampled point

n   <- length(x)       # Creating the interpolated surface
Zi  <- matrix(nrow = length(myidw$var1.pred), ncol = n)

# Removing a point then interpolate (n times for each point)
st <- stack()
for (i in 1:n){
  Z1 <- gstat::idw(he~1, x[-i,], newdata=newdat, idp=5)
  st <- addLayer(st,raster(Z1,layer=1))
  # Calculated pseudo-value Z at j
  Zi[,i] <- n * myidw$var1.pred - (n-1) * Z1$var1.pred
}

# Jackknife estimator of parameter Z at location j
Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )

# Compute (Zi* - Zj)^2 (square of residuals)
c1 <- apply(Zi,2,'-',Zj)            # Compute the difference
c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference

# Compute the confidence interval
CI <- sqrt( 1/(n*(n-1)) * c1)

# Create (CI / interpolated value) raster
img.sig   <- myidw
img.sig$v <- CI /myidw$var1.pred 

# Clip the confidence raster to he table
r4 <- raster(img.sig, layer="v")
r.m <- mask(r4, f) 

str(r4)
r5 <- crop(r4, extent(f))
plot(r5)
plot(f, add=TRUE, lwd=2) ## add F. sylvatica range

library(tmap)
FigS4b <- tm_shape(r.m) + tm_raster(n=7,title="95% confidence interval") +
  tm_shape(x) + tm_dots(size=0.1) +
  tm_legend(legend.outside=TRUE)+tm_scale_bar(position=c("left", "bottom"))

setwd("C:/Users/Camilla/Dropbox/FagusGenDiv/manuscript/Figures")

tmap_save(tm = FigS4b, filename = "FigS4.jpg", width = 174, height = 100,
          units = "mm", dpi = 300)
  dev.off()
