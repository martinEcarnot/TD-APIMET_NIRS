library(rchemo)
library(nirsextra)
library(FactoMineR) # Data analysis (PCA)
library(factoextra)  # Plots for PCA outputs

d="/home/ecarnot/Documents/INRA/enseignement/TD-APIMET_NIRS/2025/"
# MicroNIR brut
# sp=read.table(file.path(d,"MicroNIR_HD.csv"), header = TRUE, sep=";",dec=",")
# sp$x=as.matrix(sp[,7:ncol(sp)])
# sp[,7:(ncol(sp)-1)]=NULL
# attributes(sp$x)[[2]][[2]]=round(as.numeric(gsub("X","",attributes(sp$x)[[2]][[2]])))
## MicroNIR pretraité
# sp=read.table(file.path(d,"03_MicroNIR_spectra.csv"), header = TRUE, sep=";",dec=",")
sp=read.table(file.path(d,"12_Raw_MicroNIR_spectra.csv"), header = TRUE, sep=";",dec=",")
## ASD brut
# sp=read.table(file.path(d,"15_Raw_ASD_spectra.csv"), header = TRUE, sep=";",dec=",")
## NeoSpectra brut
# sp=read.table(file.path(d,"13_Raw_NeoSpectra_spectra.csv"), header = TRUE, sep=";",dec=",")
sp[, 4:ncol(sp)] <- lapply(sp[, 4:ncol(sp)], function(col) {as.numeric(trimws(col)) })
sp$x=as.matrix(sp[,4:ncol(sp)])
sp[,4:(ncol(sp)-1)]=NULL
attributes(sp$x)[[2]][[2]]=round(as.numeric(gsub("NIRS.","",attributes(sp$x)[[2]][[2]])))
sp$Date=substr(sp$Leaf_ID,6,15)


# y=read.table(file.path(d,"01_Phenotyping_traits.csv"), header = TRUE, sep=";",dec=",")
y=read.table(file.path(d,"02_Design_Greenhouse_experiment.csv"), header = TRUE, sep=";",dec=",")

sp=merge(sp,y,by="Leaf_ID")
# sp=sp[sp$Water_treatment!="WD1",]
# sp$month="Juin"
# sp$month[substr(sp$Date,6,7)=="07"]="Juillet"
# sp=sp[sp$Date!="2022-06-27" & sp$Date!="2022-06-28" & sp$Date!="2022-06-30",]
# sp=sp[sp$Date=="2022-07-25" | sp$Date=="2022-07-26",]
# sp=sp[sp$month=="Juillet",]
# sp$jour=substr(sp$Date,9,10)
str(sp)

class=as.factor(sp$Water_treatment)
plotsp(sp$x, col=class, xlab="wavelength", ylab="Reflectance", main ="Spectre Brut")
legend(x="topright", legend=unique(class), col=1:nlevels(class), lty=1)

# pca <- PCA(sp$x, graph=FALSE)
# fviz_pca_ind(pca, axes=c(1,2), habillage =class, label= "none", addEllipses = TRUE )

p = rbind(list('snv',''))
# p = rbind(list('red',c(1,1,1)),list('snv',''),list('sder',c(2,3,15)))
# p = rbind(list('sder',c(2,3,15)))
sp$xp=pre(sp$x,p)
pca <- PCA(sp$xp, graph=FALSE, ncp = 20)
i=1
fviz_pca_ind(pca, axes=c(i,i+1), habillage =class, label= "none", addEllipses = TRUE, geom="text") + ggplot2::geom_text(aes(label = sp$jour, color = class),  check_overlap = FALSE,  vjust = -1 )

# fm=plsrda(sp$xp,sp$Scenario,nlv=10)
# pred <- predict(fm, sp$xp)$pred
# err(pred, sp$Scenario)

