library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)

setwd("D:/ESAPhiWeek/")
gc()

#dump.fld <- "./dumps"
dump.fld <- "./Out_Optim"

pt.size <- 1.5
text.size <- 15

##############################################
#Figure: optimization vs optimization 3 params
##############################################

par.df.3 <- read.csv(paste(dump.fld,"Optim_S2_3Param_outdata_BroadL_corrected.csv",sep="/"))


cab.3.SCE <- ggplot(par.df.3, aes(x=Cab)) + 
  geom_point(aes(y = Cab), color = "black",size=pt.size) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_Cab), color = "darkred",size=pt.size+.5)+
  geom_smooth(aes(x=Cab,y = SCE_Cab),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()+
  theme(axis.text=element_text(size=text.size),
        axis.title=element_text(size=text.size,face="bold"))

cab.3.LBF <- ggplot(par.df.3, aes(x=Cab))+
  geom_point(aes(y = Cab), color = "black",size=pt.size) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_Cab), color="steelblue",size=pt.size+.5) +
  geom_smooth(aes(x=Cab,y = LBFGSB_Cab),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()+
  theme(axis.text=element_text(size=text.size),
        axis.title=element_text(size=text.size,face="bold"))




tiff(paste(dump.fld,"Phi_Week_3Param_Cab.tif",sep="/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))

grid.arrange(cab.3.SCE, cab.3.LBF,
             nrow = 2,ncol=2,
             top=textGrob("3 Parameters simulation - Chlorophyll a+b", gp=gpar(fontsize=18)))
dev.off()


##############################################
#Figure: optimization vs optimization 4 params
##############################################

par.df.4 <- read.csv(paste(dump.fld,"Optim_S2_4Param_outdata_BroadL_corrected.csv",sep="/"))


cab.4.SCE <- ggplot(par.df.4, aes(x=Cab)) + 
  geom_point(aes(y = Cab), color = "black",size=pt.size) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_Cab), color = "darkred",size=pt.size+.5)+
  geom_smooth(aes(x=Cab,y = SCE_Cab),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()+
  theme(axis.text=element_text(size=text.size),
        axis.title=element_text(size=text.size,face="bold"))

cab.4.LBF <- ggplot(par.df.4, aes(x=Cab))+
  geom_point(aes(y = Cab), color = "black",size=pt.size) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_Cab), color="steelblue",size=pt.size+.5) +
  geom_smooth(aes(x=Cab,y = LBFGSB_Cab),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()+
  theme(axis.text=element_text(size=text.size),
        axis.title=element_text(size=text.size,face="bold"))


tiff(paste(dump.fld,"Phi_Week_4Param_Cab.tif",sep="/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))

grid.arrange(cab.4.SCE, cab.4.LBF,
             nrow = 2,ncol=2,
             top=textGrob("4 Parameters simulation - Chlorophyll a+b", gp=gpar(fontsize=18)))
dev.off()


##############################################
#Figure: optimization vs optimization 5 params
##############################################

par.df.5 <- read.csv(paste(dump.fld,"Optim_S2_5Param_outdata_BroadL_corrected.csv",sep="/"))


cab.5.SCE <- ggplot(par.df.5, aes(x=Cab)) + 
  geom_point(aes(y = Cab), color = "black",size=pt.size) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = SCE_Cab), color = "darkred",size=pt.size+.5)+
  geom_smooth(aes(x=Cab,y = SCE_Cab),method = "lm",color = "darkred",linetype = "dashed",size=.5)+
  theme_bw()+
  theme(axis.text=element_text(size=text.size),
        axis.title=element_text(size=text.size,face="bold"))

cab.5.LBF <- ggplot(par.df.5, aes(x=Cab))+
  geom_point(aes(y = Cab), color = "black",size=pt.size) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = LBFGSB_Cab), color="steelblue",size=pt.size+.5) +
  geom_smooth(aes(x=Cab,y = LBFGSB_Cab),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  theme_bw()+
  theme(axis.text=element_text(size=text.size),
        axis.title=element_text(size=text.size,face="bold"))


tiff(paste(dump.fld,"Phi_Week_5Param_Cab.tif",sep="/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))

grid.arrange(cab.5.SCE, cab.5.LBF,
             nrow = 2,ncol=2,
             top=textGrob("4 Parameters simulation - Chlorophyll a+b", gp=gpar(fontsize=18)))
dev.off()



##################################
## Figure: random forest vs ANN
##################################
dump.ml.fld <- "./Out_MachL"

#this is better performed on the script itself because its a big dataset and the structure is unfriendly


##################################
## Figure: SCE vs ANN
##################################
dump.ml.opti.fld <- "./Out_OptimVsML"

df.optVsML <- read.csv(paste(dump.ml.opti.fld,"S2_Opti_Vs_ML.csv",sep = "/"))

p.cab <- ggplot(df.optVsML, aes(x=Cab)) + 
  #line of values
  geom_point(aes(y = Cab), color = "black",size=pt.size+.5) + 
  geom_smooth(aes(x=Cab,y = Cab),method = "lm",color = "black",linetype = "dashed",size=pt.size-.7)+
  scale_y_continuous(name="Predicted")+
  scale_x_continuous(name="Reference")+
  #line of optim
  geom_point(aes(y = SCE_Cab), color = "darkred",size=pt.size+.5)+
  geom_smooth(aes(x=Cab,y = SCE_Cab),method = "lm",color = "darkred",linetype = "dashed",size=pt.size-.7)+
  #line of ann3
  geom_point(aes(y = ann3_Cab), color = "blue4",size=pt.size+.5)+
  geom_smooth(aes(x=Cab,y = ann3_Cab),method = "lm",color = "blue4",linetype = "dashed",size=pt.size-.7)+
  #line of ann4
  #geom_point(aes(y = out.df[,4]), color = "darkorchid3",size=2)+
  #geom_smooth(aes(x=Cab,y = out.df[,4]),method = "lm",color = "darkorchid3",linetype = "dashed",size=.5)+
  ggtitle("Chlorophyll a+b")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=text.size),
        axis.title=element_text(size=text.size,face="bold"))


p.cab

p.LAI <- ggplot(df.optVsML, aes(x=LAI)) + 
  #line of values
  geom_point(aes(y = LAI), color = "black",size=pt.size+.5) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=pt.size-.7)+
  scale_y_continuous(name="Predicted")+
  scale_x_continuous(name="Reference")+
  #line of optim
  geom_point(aes(y = SCE_LAI), color = "darkred",size=pt.size+.5)+
  geom_smooth(aes(x=LAI,y = SCE_LAI),method = "lm",color = "darkred",linetype = "dashed",size=pt.size-.7)+
  #line of ann3
  geom_point(aes(y = ann3_LAI), color = "blue4",size=pt.size+.5)+
  geom_smooth(aes(x=LAI,y = ann3_LAI),method = "lm",color = "blue4",linetype = "dashed",size=pt.size-.7)+
  #line of ann4
  #geom_point(aes(y = out.df[,4]), color = "darkorchid3",size=2)+
  #geom_smooth(aes(x=LAI,y = out.df[,4]),method = "lm",color = "darkorchid3",linetype = "dashed",size=.5)+
  ggtitle("Leaf Area Index")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=text.size),
        axis.title=element_text(size=text.size,face="bold"))
p.LAI

tiff(paste(dump.ml.opti.fld,"PhiWeek_Example_S2_Opti_Vs_ML.tif",sep = "/"),
     units="px", width = 2048, height = 1024, res=124,
     compression = c("lzw"))

grid.arrange(p.cab,p.LAI,
             nrow = 2,ncol=2,
             top=textGrob("SCE vs ANN", gp=gpar(fontsize=18)))

dev.off()



##################################
## Figure: SCE vs ANN
##################################
dump.esa.results <- "./Out_ESA"

sceVSann.df <- read.csv(paste(dump.esa.results,"ESA_Preds_new.csv",sep="/"))

names(sceVSann.df)
sceVSann.df <- sceVSann.df[,c(2,5,13)]
names(sceVSann.df)<- c("LAI","SCE_LAI","ANN_LAI")

p.sce <- ggplot(sceVSann.df, aes(x=LAI))+
  geom_point(aes(y = LAI), color = "black",size=pt.size+.5) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=pt.size-.7)+
  geom_point(aes(y = SCE_LAI), color = "darkred",size=pt.size+.5) + 
  geom_smooth(aes(x=LAI,y = SCE_LAI),method = "lm",color = "darkred",linetype = "dashed",size=pt.size-.7)+
  ylab("Predicted LAI")+xlab("")+
  ggtitle("Shuffled complex evolution")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=text.size),
        axis.title=element_text(size=text.size,face="bold"))



p.mlr <- ggplot(sceVSann.df, aes(x=LAI))+
  geom_point(aes(y = LAI), color = "black",size=pt.size+.5) + 
  geom_smooth(aes(x=LAI,y = LAI),method = "lm",color = "black",linetype = "dashed",size=.5)+
  geom_point(aes(y = ANN_LAI), color="steelblue",size=pt.size+.5)+
  geom_smooth(aes(x=LAI,y = ANN_LAI),method = "lm",color = "steelblue",linetype = "dashed",size=.5)+
  ylab("Predicted LAI")+ xlab("ESA LAI")+
  ggtitle("Artificial neural network")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=text.size),
        axis.title=element_text(size=text.size,face="bold"))




tiff(paste(dump.esa.results,"PHIWeek_OVP_LAI_Prediction_Horiz.tif",sep = "/"),
     units="px", width = 1024, height = 1024, res=124,
     compression = c("lzw"))


plot(grid.arrange(p.sce,p.mlr,
                  nrow = 2,ncol=1))

dev.off()