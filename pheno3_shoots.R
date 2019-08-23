library(plyr)
library(ggplot2)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(car)
library(reshape2)
library(grid)
library(stringr)
library(RColorBrewer)
library(scales)
library(multcomp)
library(lme4)
library(lsmeans)
library(gridExtra)
library(ggdendro)
library(dendextend)
library(dendextendRcpp)
library(survminer)
library(ggridges)
library(gstat)
library(sp)

setwd("/media/jberry/Extra Drive 1/Danforth/Sorghum/Pheno3")

img_to_barcode <- read.csv("SnapshotInfo.csv",header = T,stringsAsFactors = F)
img_to_barcode <- img_to_barcode[img_to_barcode$tiles != "",]
colnames(img_to_barcode)[3] <- "Barcodes"
img_to_barcode <- img_to_barcode[,c("id","Barcodes","timestamp")]
sv_shapes <- read.table("pheno3_shapes.txt",header = F,stringsAsFactors = F,sep = " ")
colnames(sv_shapes) <- c("meta","area","hull_area","solidity","perimeter","width","height","cmx","cmy","hull_verticies","ex","ey","emajor","eminor","angle","eccen","circ","round","ar","fd","oof", "det")
sv_shapes$id <- substring(as.character(sapply(sv_shapes$meta,function(i) strsplit(i,"/")[[1]][2])),9)
sv_shapes$imgname <- as.character(sapply(sv_shapes$meta,function(i) strsplit(strsplit(i,"/")[[1]][3],"[.]")[[1]][1]))
sv_shapes <- join(sv_shapes,img_to_barcode[,c("id","Barcodes","timestamp")],by="id")
assoc <- read.csv("assoc.csv",header=T,stringsAsFactors = F)
sv_shapes <- join(sv_shapes,assoc,by="Barcodes")
sv_shapes <- sv_shapes[!is.na(sv_shapes$Drought),]
sv_shapes$timestamp <- as.POSIXct(strptime(sv_shapes$timestamp,format = "%Y-%m-%d %H:%M:%S"))
beg <- min(sv_shapes$timestamp)
sv_shapes$DAP <- floor(as.numeric((sv_shapes$timestamp - beg)/60/60/24))+2
sv_shapes$Drought <- ordered(sv_shapes$Drought, levels=c("AAA","ABB","ABA"))
sv_shapes$Genotype <- "BTx623"
sv_shapes$hour <- lubridate::hour(sv_shapes$timestamp)
sv_shapes <- sv_shapes[sv_shapes$Microbes != "Blank",]
empties <- sv_shapes[sv_shapes$DAP == 24 & sv_shapes$area == 0,]
sv_shapes <- sv_shapes[!(sv_shapes$Barcodes %in% empties$Barcodes),]
sv_shapes <- sv_shapes[(read.csv("pheno3_outliers_bool.csv",header = T,stringsAsFactors = F)$x),]
sv_shapes$Microbes[sv_shapes$Microbes == "Mix(WW+WS)"] <- "Mix"
area_convert <- 13.2*3.7/46856
tail(sv_shapes)

sv_shapes$cam_angle <- as.numeric(unlist(lapply(strsplit(sv_shapes$imgname,"_"),function(i) i[3])))
df <- aggregate(data=sv_shapes,area~Barcodes+Microbes+Drought+DAP,FUN = "mean")
ggplot(df[df$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C"),],aes(DAP,area))+
  facet_grid(~Drought)+
  geom_smooth(aes(color=Microbes),method = "loess")+
  theme_light()
head(df)
#write.csv(sv_shapes,"bart_pheno3_raw.csv",row.names = F,quote = F)

#tail(aggregate(data=sv_shapes,oof~Microbes+Drought+DAP,"sum"))
sv_shapes$Barcodes[sv_shapes$Drought == "ABB" & sv_shapes$Microbes == "SynCom A" & sv_shapes$DAP == 20]


#*************************************************************************************************
# Watering data
#*************************************************************************************************
img_to_barcode <- read.csv("SnapshotInfo.csv",header = T,stringsAsFactors = F)
img_to_barcode <- img_to_barcode[img_to_barcode$tiles == "",]
img_to_barcode$timestamp <- as.POSIXct(strptime(img_to_barcode$timestamp,format = "%Y-%m-%d %H:%M:%S"))
colnames(img_to_barcode)[3] <- "Barcodes"
img_to_barcode <- join(img_to_barcode,assoc,by="Barcodes")
beg <- min(img_to_barcode$timestamp)
img_to_barcode$DAP <- floor(as.numeric((img_to_barcode$timestamp - beg)/60/60/24))+2
head(img_to_barcode)

ggplot(img_to_barcode[img_to_barcode$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C") ,],aes(DAP,weight.before))+
  facet_grid(~Drought)+
  scale_y_continuous(limits = c(750,1500))+
  geom_point(aes(color=Microbes))

sub_bars <- unique(img_to_barcode$Barcodes[img_to_barcode$DAP >20 & img_to_barcode$weight.before>900 & img_to_barcode$Drought =="ABB"])
sv_shapes <- sv_shapes[!(sv_shapes$Barcodes %in% sub_bars),]

#*************************************************************************************************
# Outlier Detection - Done above, don't rerun
#*************************************************************************************************
library(gputools)
chooseGpu(1)
cooksd <- cooks.distance(gpuGlm(data=sv_shapes,area~Microbes:Drought:as.factor(DAP)))
sv_shapes <- sv_shapes[cooksd < 3*mean(cooksd),]
plot(cooksd)
abline(h=3*mean(cooksd),col="red")
head(cooksd)

write.csv(cooksd < 3*mean(cooksd),"pheno3_outliers_bool.csv",row.names = F,quote = F)


#*************************************************************************************************
# Looking at leaf data
#*************************************************************************************************
leaves <- setNames(read.table("pheno3_leaves.txt",stringsAsFactors = F,header=F,sep=" "),c("meta","leaf_num","area","eucD","path_length","tort"))
leaves$id <- substring(as.character(sapply(leaves$meta,function(i) strsplit(i,"/")[[1]][2])),9)
leaves$imgname <- as.character(sapply(leaves$meta,function(i) strsplit(strsplit(i,"/")[[1]][3],"[.]")[[1]][1]))
leaves <- join(leaves,img_to_barcode[,c("id","Barcodes","timestamp")],by="id")
leaves <- join(leaves,assoc,by="Barcodes")
leaves$timestamp <- as.POSIXct(strptime(leaves$timestamp,format = "%Y-%m-%d %H:%M:%S"))
leaves$DAP <- floor(as.numeric(difftime(leaves$timestamp,beg,units = "days")))+2
leaves <- leaves[!is.na(leaves$Drought),]
leaves <- leaves[leaves$imgname %in% sv_shapes$imgname,]
#leaves <- leaves[leaves$imgname %in% sv_shapes$imgname[sv_shapes$oof==0],]


for(i in 5:25){
  which_dap <- i
  back <- data.frame(with(leaves[leaves$DAP==which_dap,],table(leaf_num,Microbes,Drought))[,,"ABB"])
  back_list <- split(back,back$Microbes)
  df <- data.frame(do.call("rbind",lapply(back_list,function(i) {i$prob <- i$Freq/i$Freq[1];i})))
  
  p <- ggplot(df[df$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C"),],aes(as.numeric(leaf_num),prob))+
    geom_line(aes(color=Microbes,group=Microbes))+
    xlab("Number of leaves")+
    ylab("Proportion of Plants")+
    scale_x_continuous(limits = c(0,18),breaks = seq(0,18,2))+
    theme_light()+
    theme(axis.text = element_text(size = 14),
      axis.title= element_text(size = 18))+ 
    theme(strip.background=element_rect(fill="gray50"),
      strip.text.x=element_text(size=14,color="white"),
      strip.text.y=element_text(size=14,color="white"))
  p
  ggsave(paste0("pheno3_WS_dap",which_dap,"_numLeavesProp_syncoms.png"),width=6.01,height=4.82,plot = p, dpi = 300)
}
which_dap <- 18
back <- data.frame(with(leaves[leaves$DAP==which_dap,],table(leaf_num,Microbes,Drought))[,,"ABB"])
back_list <- split(back,back$Microbes)
df <- data.frame(do.call("rbind",lapply(back_list,function(i) {i$prob <- i$Freq/i$Freq[1];i})))

p <- ggplot(df[df$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C"),],aes(as.numeric(leaf_num),prob))+
  geom_line(aes(color=Microbes,group=Microbes))+
  xlab("Number of leaves")+
  ylab("Proportion of Plants")+
  scale_x_continuous(limits = c(0,22),breaks = seq(0,22,2))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p
ggsave(paste0("pheno3_WS_dap",which_dap,"_numLeavesProp_syncoms.png"),width=6.01,height=4.82,plot = p, dpi = 300)


leaves1 <- leaves[leaves$imgname %in% sv_shapes$imgname[sv_shapes$oof == 0],]
df <- aggregate(data=leaves1, tort~Microbes+Drought+DAP+Barcodes,FUN = "mean")
head(df)
p <- ggplot(df[df$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C") & df$DAP == 16,],aes(Microbes,tort))+
  facet_wrap(~Drought)+
  geom_boxplot()+
  xlab("Microbes")+
  ylab("Tortuosity")+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p
wilcox.test(data=df[df$Microbes %in% c("SynCom B","SynCom A") & df$DAP == 20,],log10(tort)~Microbes)


#*************************************************************************************************
# Emergence Survival
#*************************************************************************************************
des <- sort(colnames(assoc)[!(colnames(assoc) %in% "Barcodes")])
dat <- do.call("rbind",lapply(split(sv_shapes,sv_shapes$Barcodes),function(i) if(any(i$area >= 10)){
  sub <- i[i$area >= 10,]
  sub[order(sub$DAP,decreasing=F),][1,]
}else{
  sub <- i[i$DAP == max(i$DAP),][1,]
  sub[,"DAP"] <- sub[,"DAP"]+1
  sub[1,]
}))
head(dat)
dat$srv <- with(dat,Surv(time=DAP,event=(!DAP==(max(DAP)+1))))
#bad <- dat$Barcodes[dat$DAP > 5]
fmla <- as.formula(paste0("srv~",paste(des,collapse="+")))
mod1 <- summary(survfit(fmla, data = dat, conf.type = "log-log"),time=min(sv_shapes$DAP):(max(sv_shapes$DAP)))
mod_df <- data.frame("DAP"=mod1$time,"strata"=as.character(mod1$strata),"surv"=mod1$surv,"low"=mod1$lower,"high"=mod1$upper,stringsAsFactors = F)
mod_df <- cbind(mod_df,setNames(data.frame(sapply(des,function(m){unlist(lapply(str_split(mod_df$strata,", "),function(i) trimws(str_split(i[str_detect(i,m)],"=")[[1]][2])))}),stringsAsFactors = F),des))

p <- ggplot(mod_df,aes(DAP,surv))+
  facet_grid(~Drought)+
  geom_line(aes(color=Microbes))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,.2))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p


df <- na.omit(join(statefile[as.character(statefile$timestamp) == "2018-05-26",],dat,by="Barcodes"))
df <- df[,!duplicated(colnames(df))]
df$Lane1 <- df$Lane+(as.numeric(as.character(df$GH))-1)*8
df <- df[df$Drought == "AAA",]

trait <- "DAP"
sub <- df
design_effects <- aggregate(data=sub,as.formula(paste0(trait,"~Drought+Microbes")),function(i) mean(as.numeric(unlist(i))))
design_sd <- aggregate(data=sub,as.formula(paste0(trait,"~Drought+Microbes")),function(i) sd(as.numeric(unlist(i))))
sub$trait_n <- apply(sub[,c("Drought","Microbes",trait)],1,function(i){(as.numeric(i[3])-design_effects[design_effects$Drought == i[1] & design_effects$Microbes ==i[2],trait])/design_sd[design_sd$Drought == i[1] & design_sd$Microbes ==i[2],trait]})
sub$Lane1 <- sub$Lane+(as.numeric(as.character(sub$GH))-1)*8
coordinates(sub) <- ~Lane1+Position
vgm1 <- variogram(trait_n~1,sub)
fit <- fit.variogram(vgm1,model=vgm(0,"Sph"))
sub.grid <- data.frame("x"=rep(1:38,each=30),"y"=rep(1:30,38))
coordinates(sub.grid) <- ~x+y
kriged <- krige(trait_n~1,sub,sub.grid,fit,maxdist=3)
out <- setNames(data.frame(kriged)[,1:3],c("x","y",paste0("pred_",5)))
out$gh <- sapply(out$x,function(i) if(i %in% 1:8){"GH1"}else if(i %in% 9:16){"GH2"}else if(i %in% 17:24){"GH3"}else{"GH4"})
out$gh <- ordered(out$gh,levels=c("GH4","GH3","GH2","GH1"))
head(out)
  
p <- ggplot(out,aes(x,y))+
  facet_grid(~gh,space = "free_x",scales = "free_x")+
  geom_tile(aes(fill=pred_5))+
  scale_x_reverse(expand = c(0, 0),breaks=c(seq(2,38,2)))+
  scale_y_continuous(expand = c(0, 0),breaks=c(seq(2,30,2)))+
  scale_fill_gradient2(mid = "gray95",high="darkgreen",limits=c(-2.5,2.5))+
  theme_light()+
  theme(panel.spacing = unit(0.03, "lines"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=12,color="white"),
    strip.text.y=element_text(size=12,color="white"))
p  
ggsave("pheno3_aaaOnly_dap3_emerge_spatial.png",width=6.9,height=4.9,plot = p, dpi = 300)


#*************************************************************************************************
# OOF Survival
#*************************************************************************************************
head(sv_shapes)
table(sv_shapes$oof,sv_shapes$DAP)
dat <- do.call("rbind",lapply(split(sv_shapes,sv_shapes$Barcodes),function(i) if(any(i$oof == 1)){
    sub <- i[i$oof == 1,]
    sub[order(sub$DAP),][1,]
  }else{
    i[i$DAP == max(i$DAP),][1,]
  }))
head(dat)

dat$srv <- with(dat,Surv(time=DAP,event=oof))
mod1 <- summary(survfit(srv ~ Drought+Microbes, data = dat, conf.type = "log-log"),time=min(sv_shapes$DAP):max(sv_shapes$DAP))
mod_df <- data.frame("DAP"=mod1$time,"strata"=as.character(mod1$strata),"surv"=mod1$surv,"low"=mod1$lower,"high"=mod1$upper,stringsAsFactors = F)
mod_df$Drought <- unlist(lapply(strsplit(mod_df$strata,","),function(i) strsplit(trimws(i[1]),"Drought=")[[1]][2]))
mod_df$Microbes <- unlist(lapply(strsplit(mod_df$strata,","),function(i) strsplit(trimws(i[2]),"Microbes=")[[1]][2]))
tail(mod_df)

p <- ggplot(mod_df,aes(DAP,surv))+
  facet_grid(~Microbes)+
  #geom_ribbon(aes(ymin=low,ymax=high),fill="gray60")+
  geom_line(aes(color=Drought))+
  ylab("Out Of Frame Risk")+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,.2))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
ggsave("pheno4_ww_oof_risk.png",width=9.95,height=3.9,plot = p, dpi = 300)



#*************************************************************************************************
# Lighting data
#*************************************************************************************************
p <- ggplot(sv_shapes,aes(hour,det))+
  geom_jitter(width = 0.5)+
  scale_x_continuous(breaks = seq(from=0,to=24,by=2))+
  scale_y_continuous(limits=c(-200,2))+
  ylab("Correction Strength")+
  xlab("Time of Day (hr)")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 14),
        axis.title= element_text(size = 18))+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(legend.position='none')
p
ggsave("pheno3_image_variability.png",width=4.89,height=4.69,plot = p, dpi = 300)


#*************************************************************************************************
# Quick plots
#*************************************************************************************************
sub <- sv_shapes
p <- ggplot(sub[sub$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C") ,],aes(DAP,area*area_convert))+
  facet_wrap(~Drought)+
  geom_smooth(aes(color=Microbes),method="loess")+
  ylab(~~Area~(cm^2))+
  scale_color_manual(values=c("gray20",muted("red",60,100),muted("green",60,100),muted("cyan",60,100)))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
        axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno3_area_syncom_trends.png",width=8.03,height=3.9,plot = p, dpi = 300)

p <- ggplot(sv_shapes[sv_shapes$DAP == 20,],aes(Microbes,area*area_convert))+
  facet_wrap(~Drought)+
  ggtitle("DAP = 20")+
  geom_boxplot(aes(color=Microbes))+
  ylab(~~Area~(cm^2))+
  xlab("")+
  #scale_color_manual(values=c("gray20",muted("red",60,100),muted("green",60,100),muted("cyan",60,100)))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("pheno3_area_syncom_dap20_boxplots.png",width=7.5,height=3.9,plot = p, dpi = 300)
test <- aggregate(data=sv_shapes,area~Drought+Microbes+DAP,FUN="mean")
head(test)
lapply(split(test[test$DAP == 25,],test$Drought),function(i) paste(i[order(i$area),"Microbes"],collapse = " < "))
b[order(b$`(area * area_convert)`,decreasing = T),]

sub <- sv_shapes[sv_shapes$Barcodes== "Fm001AB056588" & str_detect(sv_shapes$imgname,"VIS_SV_0"),]
sub$area_wave <- wavethresh::wd(sub$area,2)
p <- ggplot(sub,aes(DAP,area*area_convert))+
  geom_line()+
  geom_point()+
  scale_y_continuous(limits = c(0,150))+
  ylab(~~Area~(cm^2))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p


#*************************************************************************************************
# Logistic Example
#*************************************************************************************************
df <- data.frame("DAP"=0:25,"Area"=sapply(0:25,function(i) 104.23/(1+(1-0.001)/0.001*exp(-0.41*i))))
p <- ggplot(df,aes(DAP,Area))+
  geom_line()+
  geom_jitter(height = 6,size=2)+
  scale_y_continuous(limits = c(0,120))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p
ggsave("logistic_example.png",width=5.17,height=4.62,plot = p, dpi = 300)


#*************************************************************************************************
# R^2 and Parital Correlations
#*************************************************************************************************
dat <- sv_shapes[sv_shapes$DAP==23,]
tail(dat)
shapes <- c("area","hull_area","solidity","perimeter","width","height","cmx","cmy","hull_verticies","ex","ey","emajor","eminor","angle","eccen","circ","round","ar","fd", "det")
H2 <- c()
for(e in shapes){
  model <- lmer(eval(parse(text=e))~(1|Drought)+(1|Microbes)+(1|Drought:Microbes),data = dat)
  re<-as.numeric(VarCorr(model))
  res<-attr(VarCorr(model), "sc")^2
  interaction.var <- re[1]
  microbe.var<-re[2]
  drought.var<-re[3]
  tot.var<-sum(re,res)
  unexp <- 1-sum(re)/sum(re,res)
  
  h2 <- c((microbe.var/tot.var),
          (drought.var/tot.var),
          (interaction.var/tot.var),
          unexp)
  H2 <- rbind(H2,h2)
}
H2 <- data.frame(H2,row.names = shapes)
H2$Shape <- rownames(H2)
rownames(H2) <- NULL
colnames(H2) <- c("Microbe","Drought","Interaction","Unexplained","Shape")
H2$Shape <-  ordered(H2$Shape,levels=H2$Shape[order(H2$Unexplained)])
H2_melt <- melt(H2,id=c("Shape"))
H2_melt$variable <- ordered(H2_melt$variable,levels=c("Unexplained","Drought","Microbe","Interaction"))
head(H2_melt)
p <- ggplot(data=H2_melt,aes(Shape,value*100))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.text = element_text(size = 14),
        axis.title.y= element_text(size = 18),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill=NA, size=1,linetype = 1))+
  theme(legend.position = "top")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("pheno3_dap23_anova.png",width=6.17,height=4.54,plot = p, dpi = 300)


#*************************************************************************************************
# First Derivative Plot
#*************************************************************************************************
head(sv_shapes)
deriv <- do.call("rbind",lapply(split(sv_shapes,sv_shapes$Barcodes),function(i) aggregate(data=i,area~DAP+Barcodes+Drought+Microbes,FUN = "mean")))
deriv$diff <- unlist(lapply(split(deriv,deriv$Barcodes),function(i) c(0,diff(i$area))))
rownames(deriv) <- NULL

p <- ggplot(deriv,aes(DAP,diff))+
  facet_wrap(~Drought)+
  geom_smooth(aes(color=Microbes),method = "loess",se=F)+
  ylab("dArea / dDAP")+
  theme_light()+
  theme(axis.text = element_text(size = 14),
        axis.title= element_text(size = 18))+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno3_area_first_derivative.png",width=8.03,height=3.9,plot = p, dpi = 300)


#*************************************************************************************************
# Positional effects
#*************************************************************************************************
statefile <- NULL
for(my_dir in list.dirs("Statefiles")[-1]){
  for(i in list.files(my_dir)){
    temp <- read.table(paste0(substr(my_dir,0,20),"/",i),header = F,stringsAsFactors = F,sep = ";",skip = 1)
    num <- as.numeric(read.table(paste0(substr(my_dir,0,20),"/",i),header = F,stringsAsFactors = F,sep = ";",nrows = 1))
    meta <- substr(strsplit(i,"[.]")[[1]][3],5,7)
    statefile <- rbind(statefile,data.frame("Date"=substr(my_dir,12,20),"GH"=strsplit(meta,"_")[[1]][1],"Lane"=as.numeric(strsplit(meta,"_")[[1]][2]),"Position"=1:num,"Barcodes"=temp$V3,stringsAsFactors = F))
  }
}
statefile$GH <- ordered(statefile$GH,levels=5:1)
statefile$timestamp <- as.POSIXlt(strptime(statefile$Date,format = "%m%d%Y"))
statefile$Date <- paste0(lubridate::month(statefile$timestamp),lubridate::day(statefile$timestamp),lubridate::year(statefile$timestamp))

dat <- aggregate(data=sv_shapes,area~Barcodes+timestamp+Drought+Microbes+DAP,"mean")
dat$timestamp <- as.POSIXct(strptime(dat$timestamp,format = "%Y-%m-%d %H:%M:%S"))
#dat <- aggregate(data=nir,intensityAVG~Barcodes+as.character(timestamp)+Drought+Microbes+DAP,"mean")
#dat$timestamp <- as.POSIXct(strptime(dat$`as.character(timestamp)`,format = "%Y-%m-%d %H:%M:%S"))
dat$Date <- paste0(lubridate::month(dat$timestamp),lubridate::day(dat$timestamp),lubridate::year(dat$timestamp))
dat <- join(statefile,dat,by=c("Barcodes","Date"))
dat <- na.omit(dat[,-c(6,7)])

trait <- "area"
all_spat <- data.frame("x"=rep(1:38,each=30),"y"=rep(1:30,38))
for(day in 9:22){
  sub <- dat[dat$DAP == day & dat$Drought %in% c("ABB","AAA"),]
  design_effects <- aggregate(data=sub,as.formula(paste0(trait,"~Drought+Microbes")),function(i) mean(as.numeric(unlist(i))))
  design_sd <- aggregate(data=sub,as.formula(paste0(trait,"~Drought+Microbes")),function(i) sd(as.numeric(unlist(i))))
  sub$trait_n <- apply(sub[,c("Drought","Microbes",trait)],1,function(i){(as.numeric(i[3])-design_effects[design_effects$Drought == i[1] & design_effects$Microbes ==i[2],trait])/design_sd[design_sd$Drought == i[1] & design_sd$Microbes ==i[2],trait]})
  sub$Lane1 <- sub$Lane+(as.numeric(as.character(sub$GH))-1)*8
  coordinates(sub) <- ~Lane1+Position
  vgm1 <- variogram(trait_n~1,sub)
  fit <- fit.variogram(vgm1,model=vgm(0,"Sph"))
  sub.grid <- data.frame("x"=rep(1:38,each=30),"y"=rep(1:30,38))
  coordinates(sub.grid) <- ~x+y
  kriged <- krige(trait_n~1,sub,sub.grid,fit,maxdist=3)
  out <- setNames(data.frame(kriged)[,1:3],c("x","y",paste0("pred_",day)))
  all_spat <- join(all_spat,out,by=c("x","y"))
}
all_spat$avg <- rowMeans(all_spat[,3:(ncol(all_spat))],na.rm = T)
all_spat$var <- apply(all_spat[,3:(ncol(all_spat)-1)],1,function(i) var(i,na.rm = T))
all_spat$cv <- all_spat$var/all_spat$avg
all_spat$gh <- sapply(all_spat$x,function(i) if(i %in% 1:8){"GH1"}else if(i %in% 9:16){"GH2"}else if(i %in% 17:24){"GH3"}else{"GH4"})
all_spat$gh <- ordered(all_spat$gh,levels=c("GH4","GH3","GH2","GH1"))
head(all_spat)

p <- ggplot(all_spat,aes(x,y))+
  facet_grid(~gh,space = "free_x",scales = "free_x")+
  geom_tile(aes(fill=avg))+
  scale_x_reverse(expand = c(0, 0),breaks=c(seq(2,38,2)))+
  scale_y_continuous(expand = c(0, 0),breaks=c(seq(2,30,2)))+
  scale_fill_gradient2(mid = "gray95",high="darkgreen",limits=c(-2.5,2.5))+
  theme_light()+
  theme(panel.spacing = unit(0.03, "lines"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=12,color="white"),
    strip.text.y=element_text(size=12,color="white"))
p  
ggsave("pheno3_bothabb-aaa_allDAP_areaNorm_spatial.png",width=6.9,height=4.9,plot = p, dpi = 300)

sub <- dat[dat$DAP == 20 & dat$Drought == "ABB",]
sub$Lane1 <- sub$Lane+(as.numeric(as.character(sub$GH))-1)*8
colnames(sub)[4] <- "y"
colnames(sub)[10] <- "x"

df <- join(sub,all_spat[,c("x","y","avg")],by=c("x","y"))
df$status <- "Uncalibrated"
df1 <- df
df1$area <- df1$area-df1$area*df1$avg*0.5
df1$status <- "Calibrated"
df <- rbind(df,df1)
df$status <- ordered(df$status, levels=c("Uncalibrated","Calibrated"))

ggplot(df,aes(Microbes,area*area_convert))+
  facet_grid(~status)+
  geom_boxplot()

design_effects <- aggregate(data=df,as.formula(paste0(trait,"~Drought+Microbes")),function(i) mean(as.numeric(unlist(i))))
design_sd <- aggregate(data=df,as.formula(paste0(trait,"~Drought+Microbes")),function(i) sd(as.numeric(unlist(i))))
df$trait_n <- apply(df[,c("Drought","Microbes",trait)],1,function(i){(as.numeric(i[3])-design_effects[design_effects$Drought == i[1] & design_effects$Microbes ==i[2],trait])/design_sd[design_sd$Drought == i[1] & design_sd$Microbes ==i[2],trait]})
coordinates(df) <- ~x+y
vgm1 <- variogram(trait_n~1,df)
fit <- fit.variogram(vgm1,model=vgm(1,"Sph"))
df1.grid <- data.frame("x"=rep(1:38,each=30),"y"=rep(1:30,38))
coordinates(df1.grid) <- ~x+y
kriged <- krige(trait_n~1,df,df1.grid,fit,maxdist=3)
out <- setNames(data.frame(kriged)[,1:3],c("x","y","pred"))

p <- ggplot(out,aes(x,y))+
#  facet_grid(~GH,space = "free_x",scales = "free_x")+
  geom_tile(aes(fill=pred))+
  scale_x_reverse(expand = c(0, 0),breaks=c(seq(2,38,2)))+
  scale_y_continuous(expand = c(0, 0),breaks=c(seq(2,30,2)))+
  scale_fill_gradient2(mid = "gray95",high="darkgreen",limits=c(-2.5,2.5))+
  theme_light()
p


with(df[df$status == "Uncalibrated" & df$Microbes %in% c("SynCom A","SynCom B"),],t.test(area~Microbes))


my_list <- do.call(rbind,lapply(unique(dat$DAP),function(i){
  sub <- dat[dat$DAP == i & dat$Drought == "ABB",]
  sub$Lane1 <- sub$Lane+(as.numeric(as.character(sub$GH))-1)*8
  colnames(sub)[4] <- "y"
  colnames(sub)[10] <- "x"
  df <- join(sub,all_spat[,c("x","y","avg")],by=c("x","y"))
  df$area_c <- df$area-df$area*df$avg*0.5
  df
}))

ggplot(my_list,aes(area*area_convert,area_c*area_convert))+
  geom_point(aes(color=avg),size=3)+
  scale_color_gradient2(mid = "gray95",high="darkgreen")+
  geom_abline(slope=1, intercept=0,color="gray20")+
  theme_light()

ggplot(my_list,aes(DAP,area*area_convert))+
  geom_smooth(aes(color=Microbes),method = "loess")


#*************************************************************************************************
# Which model is best
#*************************************************************************************************
source("analysis.R")
source("modeling.models.R")

head(sv_shapes)
control.growth1 <- control.plant.growth.modeling(input=sv_shapes,
                                                 col.plant="Barcodes",
                                                 col.genotype="Microbes",
                                                 col.treatment = "Drought",
                                                 treatment = "AAA",
                                                 genotype.id = "Control",
                                                 col.day="DAP",
                                                 col.time="DAP",
                                                 proxy.trait="area")
control.growth1 <- control.plant.growth.modeling(input=sv_shapes[sv_shapes$DAP %in% 1:18,],
                                                 col.plant="Barcodes",
                                                 col.genotype="Microbes",
                                                 col.treatment = "Drought",
                                                 treatment = "ABA",
                                                 plant.id = "Fm001BB056694",
                                                 col.day="DAP",
                                                 col.time="DAP",
                                                 proxy.trait="area")
stress.growth2<-stressed.plant.growth.modeling(input=sv_shapes,
                                               col.plant="Barcodes",
                                               col.genotype="Microbes",
                                               col.treatment = "Drought",
                                               treatment = "ABA",
                                               genotype.id = "Control",
                                               col.day="DAP",
                                               col.time="DAP",
                                               proxy.trait="area",
                                               first.stress.day = 4,
                                               last.stress.day = 18)
stress.growth2<-stressed.plant.growth.modeling(input=sv_shapes,
                                               col.plant="Barcodes",
                                               col.genotype="Microbes",
                                               col.treatment = "Drought",
                                               treatment = "ABA",
                                               plant.id = "Fm001BB056694",
                                               col.day="DAP",
                                               col.time="DAP",
                                               first.stress.day = 4,
                                               last.stress.day = 18,
                                               proxy.trait="area",
                                               plot = T)


#*************************************************************************************************
# Logistic Growth Modeling for AAA
#*************************************************************************************************
source("analysis.R")
source("modeling.models.R")

back <- NULL
for(i in unique(sv_shapes$Barcodes[sv_shapes$Drought %in% c("AAA")])){
  control.growth1 <- control.plant.growth.modeling(input=sv_shapes,
                                                   col.plant="Barcodes",
                                                   col.genotype="Microbes",
                                                   col.treatment = "Drought",
                                                   treatment = sv_shapes$Drought[sv_shapes$Barcodes == i][1],
                                                   plant.id = i,
                                                   col.day="DAP",
                                                   col.time="DAP",
                                                   proxy.trait="area",
                                                   plot = F)
  if(!is.null(control.growth1$Inflection)){
    back <- rbind(back,data.frame("Barcodes" = i,
                                  "Inflection_point" = control.growth1$Inflection,
                                  "Inflection_rate" = control.growth1$InflectionGR,
                                  "Amax" = control.growth1$Kmax,
                                  "R2" = control.growth1$logistic$rsquared,
                                  stringsAsFactors = F))
  }
}
back <- join(back,assoc,by="Barcodes")
back$Microbes[back$Microbes == "Mix(WW+WS)"] <- "Mix"
back$Microbes <- ordered(back$Microbes, levels=c("Control","WW","WS","Mix","SynCom A","SynCom B","SynCom C"))
back <- back[back$R2 > 0.9,]


#*************************************************************************************************
# Two-stage modeling for ABA
#*************************************************************************************************
rec <- NULL
for(i in unique(sv_shapes$Barcodes[sv_shapes$Drought %in% c("ABA")])){
  if(i %in% c("Fm001BB057240","Fm001GB057392","Fm001BB056736","Fm001BB056799","Fm001GB056615","Fm001AB056378","Fm001GB056447")){next}
  stress.growth1<-control.plant.growth.modeling(input=sv_shapes[sv_shapes$DAP < 18,],
                                                 col.plant="Barcodes",
                                                 col.genotype="Microbes",
                                                 col.treatment = "Drought",
                                                 treatment = sv_shapes$Drought[sv_shapes$Barcodes == i][1],
                                                 plant.id = i,
                                                 col.day="DAP",
                                                 col.time="DAP",
                                                 proxy.trait="area",
                                                 plot = F)
  stress.growth2<-control.plant.growth.modeling(input=sv_shapes[sv_shapes$DAP >= 18,],
                                                 col.plant="Barcodes",
                                                 col.genotype="Microbes",
                                                 col.treatment = "Drought",
                                                 treatment = sv_shapes$Drought[sv_shapes$Barcodes == i][1],
                                                 plant.id = i,
                                                 col.day="DAP",
                                                 col.time="DAP",
                                                 proxy.trait="area",
                                                 plot = F)
  if(!(is.null(stress.growth1$Inflection)|is.null(stress.growth2$Inflection))){
    rec <- rbind(rec,data.frame("Barcodes" = i,
                                  "Inflection_point_d" = stress.growth1$Inflection,
                                  "Inflection_rate_d" = stress.growth1$InflectionGR,
                                  "Amax_d" = stress.growth1$Kmax,
                                  "R2_d" = stress.growth1$logistic$rsquared,
                                  "Inflection_point_r" = stress.growth2$Inflection,
                                  "Inflection_rate_r" = stress.growth2$InflectionGR,
                                  "Amax_r" = stress.growth2$Kmax,
                                  "R2_r" = stress.growth2$logistic$rsquared,
                                  stringsAsFactors = F))
  }
}
rec <- join(rec,assoc,by="Barcodes")
rec$Microbes[rec$Microbes == "Mix(WW+WS)"] <- "Mix"
rec$Microbes <- ordered(rec$Microbes, levels=c("Control","WW","WS","Mix","SynCom A","SynCom B","SynCom C"))
rec <- rec[rec$R2_r > 0.9,]
rec$IRR <- rec$Inflection_rate_r/rec$Inflection_rate_d
anova(lm(data=rec,IRR~Microbes))
p <- ggplot(rec,aes(Microbes,IRR))+
  facet_grid(~Drought)+
  geom_boxplot()+
  ylab("Recovery / Stress")+
  theme_light()+
  theme(axis.text = element_text(size = 12),
    axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
    strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("pheno3_ABA_rate_ratio.png",width=3.58,height=5.47,plot = p, dpi = 300)

tail(rec)
aggregate(data=back,Barcodes~Microbes,FUN=function(i)length(i))

test <- rec[rec$R2_r > 0.9 & rec$R2_d > 0.9,]


#*************************************************************************************************
# Stress elasticity and stress susceptibility index
#*************************************************************************************************
head(back)
tail(rec)

#combinding back and rec
b <- data.frame("Microbes"=c(as.character(back$Microbes),as.character(rec$Microbes)),"Drought"=c(back$Drought,rec$Drought),"Amax"=c(back$Amax,rec$Amax_r),"Infl"=c(back$Inflection_rate,rec$Inflection_rate_r), stringsAsFactors = F)
b$Microbes <- ordered(b$Microbes, levels=c("Control","WW","WS","Mix","SynCom A","SynCom B","SynCom C"))
p <- ggplot(b,aes(Drought,Amax*area_convert))+
  facet_grid(~Microbes)+
  geom_boxplot()+
#  scale_y_continuous(limits=c(7500,18000),oob=rescale_none)+ #for infl
  scale_y_continuous(limits=c(100,275),oob=rescale_none)+ #forAmax
  ylab("Amax")+
  theme_light()+
  theme(axis.text = element_text(size = 12),
    axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
    strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("pheno3_logistic_amax.png",width=8.4,height=4.42,plot = p, dpi = 300)

#stress elasticity 
test <- aggregate(data=b[b$Drought %in% c("AAA","ABA"),],Infl~Microbes+Drought,FUN = function(i){mean(i,na.rm=T)})
elas <- data.frame("Microbes"=test$Microbes[test$Drought == "ABA"],"Elas"=test$Infl[test$Drought == "ABA"]/test$Infl[test$Drought == "AAA"],stringsAsFactors = F)

p <- ggplot(elas,aes(Microbes,Elas))+
  geom_segment(aes(x=Microbes,y=0,xend=Microbes,yend=Elas))+
  geom_point(size=4)+
  ylab("Stress Elasticity")+
  scale_y_continuous(limits=c(0.85,1.15),breaks = seq(0.85,1.15,0.1),oob=rescale_none)+
  theme_light()+
  theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
        strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("pheno3_stress_elasticity.png",width=3.31,height=4.39,plot = p, dpi = 300)

#SI
head(sv_shapes)
test <- aggregate(data=sv_shapes[sv_shapes$DAP == 24,],area~Microbes+Drought,FUN = "mean")
si <- data.frame("Microbes"=test$Microbes[test$Drought == "ABB"],"SI"=1-test$area[test$Drought=="ABB"]/test$area[test$Drought=="AAA"],stringsAsFactors = F)
si$Microbes <- ordered(si$Microbes, levels=c("Control","WW","WS","Mix","SynCom A","SynCom B","SynCom C"))

p <- ggplot(si,aes(Microbes,SI))+
  geom_segment(aes(x=Microbes,y=0,xend=Microbes,yend=SI))+
  geom_point(size=4)+
  ylab("Stress Intensity")+
  scale_y_continuous(limits=c(0.2,0.6),breaks = seq(0.2,0.6,0.1),oob=rescale_none)+
  theme_light()+
  theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
        strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("pheno3_stress_intensity.png",width=3.31,height=4.39,plot = p, dpi = 300)


#*************************************************************************************************
# Helper functions
#*************************************************************************************************
get_color <- function(file_name,start,stop){
  color_data <- read.table(file_name,header = F,stringsAsFactors = F,sep = " ")[,-257]
  color_data$id <- as.character(sapply(color_data$V1,function(i) strsplit(strsplit(i,"/")[[1]][2],"snapshot")[[1]][2]))
  color_data$imgname <- as.character(sapply(color_data$V1,function(i) strsplit(strsplit(i,"/")[[1]][3],"[.]")[[1]][1]))
  color_data <- join(color_data,img_to_barcode[,c("id","Barcodes","timestamp")],by="id")
  color_data <- join(color_data,assoc,by="Barcodes")
  color_data$timestamp <- strptime(color_data$timestamp,format = "%Y-%m-%d %H:%M:%S")
  color_data$DAP <- floor(as.numeric((color_data$timestamp - beg)/60/60/24))+2
  color_data[,start:stop] <- t(apply(color_data[,start:stop],1,function(i){i/(sum(i,na.rm = T)+1)}))*100
  color_data$hr <- as.POSIXlt(color_data$timestamp)$hour
  #color_data$status <- sapply(color_data$hr,function(i) if(i %in% 7:21){"Daytime"}else{"Nighttime"})
  return(color_data)
}

hist_avg <- function(data,start,stop){
  sub <- data
  test <- data.frame(do.call("rbind",lapply(split(sub,sub$Drought),function(t){
    data.frame(do.call("rbind",lapply(split(t,t$Microbes),function(g){
      data.frame(do.call("rbind",lapply(split(g,g$DAP),function(m){
        colMeans(m[,start:stop],na.rm = T)
      }
      )))
    }
    )))
  })))
  return(test)
}

hist_sd <- function(data,start,stop){
  sub <- data
  test <- data.frame(do.call("rbind",lapply(split(sub,sub$Drought),function(t){
    data.frame(do.call("rbind",lapply(split(t,t$Microbes),function(g){
      data.frame(do.call("rbind",lapply(split(g,g$DAP),function(m){
        apply(m[,start:stop],2,function(i){sd(i,na.rm = T)})
      }
      )))
    }
    )))
  })))
  return(test)
}


#*************************************************************************************************
# NIR Data
#*************************************************************************************************
nir <- get_color("nir_color.txt",2,256)
nir$intensityAVG <- apply(nir[,3:255],1,function(i){sum((i/100)*(2:254),na.rm = T)})
nir$Microbes[nir$Microbes == "Mix(WW+WS)"] <- "Mix"
nir$Microbes <- ordered(nir$Microbes, levels=c("Control","WW","WS","Mix","SynCom A","SynCom B","SynCom C"))
nir <- nir[!is.na(nir$Microbes),]
nir <- nir[!(nir$Barcodes %in% empties),]
outliers <- read.csv("pheno3_area_outliers.csv",header = T,stringsAsFactors = F)
outliers$camera_angle <- unlist(lapply(strsplit(outliers$meta,"_"),function(i) i[3]))
outliers$unique_id <- paste(outliers$Barcodes,outliers$DAP,outliers$camera_angle,sep="_")
nir$camera_angle <- unlist(lapply(strsplit(as.character(nir$V1),"_"),function(i) i[3]))
nir$unique_id <- paste(nir$Barcodes,nir$DAP,nir$camera_angle,sep="_")
nir <- nir[!(nir$unique_id %in% outliers$unique_id),]
cooksd <- cooks.distance(lm(data=nir,intensityAVG~Microbes:Drought:as.factor(DAP)))
nir <- nir[cooksd < 3*mean(cooksd),]


#*************************************************************************************************
# Joyplot
#*************************************************************************************************
joyplot_response <- function(data,start,stop){
  sub <- data
  test <- hist_avg(sub,start,stop)
  test_sd <- hist_sd(sub,start,stop)
  
  test_avg <- do.call("rbind",lapply(split(test,rownames(test)),function(j){colMeans(j,na.rm = T)}))
  test_sd <- do.call("rbind",lapply(split(test_sd,rownames(test_sd)),function(j){colMeans(j,na.rm = T)}))
  test_avg <- data.frame(melt(test_avg))
  test_avg$sd <- data.frame(melt(test_sd))[,3]
  test_avg$Drought <- unlist(lapply(strsplit(as.character(test_avg$Var1),"[.]"),function(i)i[1]))
  test_avg$Microbes <- unlist(lapply(strsplit(as.character(test_avg$Var1),"[.]"),function(i)i[2]))
  test_avg$DAP <- unlist(lapply(strsplit(as.character(test_avg$Var1),"[.]"),function(i)i[3]))
  test_avg$Microbes <- ordered(test_avg$Microbes, levels=c("Control","WW","WS","Mix(WW+WS)","SynCom A","SynCom B","SynCom C"))
  test_avg <- test_avg[as.numeric(str_sub(test_avg$Var2,2,4)) <= 256,]
  test_avg$bin <- as.numeric(str_sub(test_avg$Var2,2,4))-2
  limits <- aes(ymax = value + 1.96*sd, ymin=value - 1.96*sd)
  sub <- test_avg[as.numeric(test_avg$DAP) > 0,]

  return(sub)
}
back <- joyplot_response(nir,2,256)
back$Drought <- ordered(back$Drought,levels=rev(c("ABA","AAA","ABB")))
back <- back[as.numeric(back$DAP) > 4,]
head(back)

p <- ggplot(data=back,aes(x=bin,y=factor(DAP,levels = 25:1), height=value))+
  facet_grid(~Microbes)+
  #geom_ribbon(limits,aes(fill=Drought))+
  geom_density_ridges(stat = "identity", aes(fill=Drought),alpha=0.5)+
  scale_fill_manual(values = rev(c(muted("green",c=100,l=60),muted("red",c = 100,l=60),"transparent")))+
  ylab("Days after planting")+
  xlab("")+
  xlab("NIR Grayscale Intensity")+
  theme_ridges(grid=F,center_axis_labels = T)+
  #theme(legend.position='none')+
  theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
        strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno3_nir_joy.png",width=11,height=6,plot = p, dpi = 300)
ggsave("pheno3_nir_histograms_dap25.png",width=7.59,height=7.14,plot = p, dpi = 300)


#*************************************************************************************************
# NIR Average Intensity
#*************************************************************************************************
#nir$intensityAVG <- apply(nir[,3:255],1,function(i){sum((i/100)*(2:254),na.rm = T)})

#trends
p <- ggplot(nir[nir$intensityAVG != 0 & nir$DAP >8,],aes(DAP,intensityAVG))+
  facet_grid(~Microbes)+
  geom_smooth(aes(color=Drought),method = "loess")+
#  geom_point(aes(color=Drought))+
#  scale_y_continuous(limits = c(75,80))+
#  geom_vline(xintercept = 4,linetype="dashed")+
  geom_vline(xintercept = 18,linetype="dashed")+
  theme_light()+
  theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
        strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno3_nir_intesityAVG_trends_dapg8.png",width=12.3,height=4.18,plot = p, dpi = 300)

#heatmap just of drought
test <- aggregate(data=nir[nir$intensityAVG != 0 & nir$DAP >8,],intensityAVG~Drought+DAP,FUN = function(i)mean(i,na.rm=T))
p <- ggplot(test,aes(DAP,Drought))+
  geom_tile(aes(fill=intensityAVG))+
  scale_fill_gradient2(limits=c(75,95),midpoint = mean(nir$intensityAVG[nir$Drought == "AAA" & nir$intensityAVG != 0 & nir$DAP == 15]),high ="gray10",low= "#56B1F7",mid = "#d7e4ef")+
  theme_light()+
  theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
        strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno3_nir_intesityAVG_heatmap_dapg8.png",width = 8.91,height = 3.93,plot = p, dpi = 300)

#heatmap broken out by microbes
test <- aggregate(data=nir[nir$intensityAVG != 0 & nir$DAP >8,],intensityAVG~Drought+Microbes+DAP,FUN = function(i)mean(i,na.rm=T))
b <- aggregate(data=sv_shapes[sv_shapes$DAP==25 & sv_shapes$Drought == "ABB",],area~Microbes,FUN="mean")
b$Microbes[b$Microbes == "Mix(WW+WS)"] <- "Mix"
test$Microbes <- ordered(test$Microbes,levels=b$Microbes[order(b$area,decreasing = F)])
test <- test[test$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C"),]
test$Microbes <- ordered(test$Microbes,levels=rev(c("Control","SynCom A","SynCom B","SynCom C")))
p <- ggplot(test,aes(DAP,Microbes))+
  facet_grid(~Drought)+
  geom_tile(aes(fill=intensityAVG))+
  scale_fill_gradient2(limits=c(75,95),midpoint = 0.3+mean(nir$intensityAVG[nir$Drought == "AAA" & nir$intensityAVG != 0 & nir$DAP == 15]),high ="gray10",low= "#56B1F7",mid = "#d7e4ef")+
  theme_light()+
  theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
        strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno3_nir_intesityAVG_heatmap_drought_dap25_syncoms_alphabetOrdered.png",width = 8.91,height = 3.93,plot = p, dpi = 300)

#NIR ~ Biomass
test <- aggregate(data=nir[nir$intensityAVG != 0 & nir$DAP >8,],intensityAVG~Microbes+Drought+DAP,FUN = function(i)mean(i,na.rm=T))
b <- aggregate(data=sv_shapes,area~Drought+Microbes+DAP,FUN="mean")
b$Microbes[b$Microbes == "Mix(WW+WS)"] <- "Mix"
df <- join(b,test)
p <- ggplot(df[df$DAP > 8,],aes(area*area_convert,intensityAVG))+
  facet_wrap(~DAP,scales = "free")+
  geom_smooth(method = "lm",se=F,color="gray20")+
  geom_point(aes(color=Microbes,shape=Drought),size=3)+
  xlab(~~Area~(cm^2))+
  theme_light()+
  theme(axis.text = element_text(size = 12),
    axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
    strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno3_nir_intesityAVG_biomass_allDAP_allDrought.png",width = 10,height = 8,plot = p, dpi = 300)
head(df)

#boxplot
p <- ggplot(nir[nir$DAP==25 & nir$intensityAVG != 0,],aes(Drought,intensityAVG))+
  facet_grid(~Microbes)+
  geom_boxplot()+
  theme_light()+
  #scale_y_continuous(limits=c(65,105))+
  scale_fill_gradient(limits=c(35,251),low = "white",high = "black")+
  theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
        strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("pheno3_nir_intesityAVG_dap25_boxplots.png",width=10.3,height=4.18,plot = p, dpi = 300)


#*************************************************************************************************
# Associating NIR intensity to stress measures 
#*************************************************************************************************
head(si)
head(elas)
nir_i <- aggregate(data=nir[nir$DAP == 22,],intensityAVG~Microbes+Drought,FUN = "mean")

df <- rbind(join(join(nir_i[nir_i$Drought == "AAA",],elas),si),
            join(join(nir_i[nir_i$Drought == "ABA",],elas),si),
            join(join(nir_i[nir_i$Drought == "ABB",],elas),si))

ggplot(df,aes(SI,intensityAVG))+
  geom_point(aes(color=Drought),size=3)+
  geom_smooth(aes(color=Drought),method = "loess",span=2)+
  theme_light()

summary(lm(data=df[df$Drought == "AAA",],intensityAVG~rescale(Elas)))



nir$id <- substring(as.character(sapply(nir$V1,function(i) strsplit(i,"/")[[1]][2])),9)

df <- join(sv_shapes,nir,by="id")
plot(df$area,df$intensityAVG)


#*************************************************************************************************
# VIS Data
#*************************************************************************************************
vis <- get_color("color.txt",2,182)
outliers <- read.csv("pheno3_area_outliers.csv",header = T,stringsAsFactors = F)
outliers$camera_angle <- unlist(lapply(strsplit(outliers$meta,"_"),function(i) i[3]))
outliers$unique_id <- paste(outliers$Barcodes,outliers$DAP,outliers$camera_angle,sep="_")
vis$camera_angle <- unlist(lapply(strsplit(as.character(vis$V1),"_"),function(i) i[3]))
vis$unique_id <- paste(vis$Barcodes,vis$DAP,vis$camera_angle,sep="_")
vis <- vis[!(vis$unique_id %in% outliers$unique_id),]
vis <- vis[vis$Microbes != "Blank",]
head(vis)

ks <- vis[vis$DAP == 15,]
ks[,2:181] <- t(apply(ks[,2:181],1,function(i)rescale(cumsum(i))))
ks <- ks[!is.na(ks$Microbes),]
ks_melt <- melt(ks[,c(2:181,which(colnames(ks) %in% c("Barcodes","Microbes","Drought","DAP")))],id=c("Barcodes","Microbes","Drought","DAP"))
ks_melt$bin <- as.numeric(str_sub(ks_melt$variable,2,4))-2
ks_melt <- aggregate(data=ks_melt,value~bin+Barcodes+Microbes+Drought,"mean")
head(ks_melt)
ggplot(data=ks_melt[ks_melt$value != 0.5 & ks_melt$Microbes %in% c("Control","WW"),], aes(bin,value))+
  facet_wrap(~Microbes)+
  geom_smooth(aes(color=Drought),method = "loess",span=0.1,se=F)+
#  scale_x_continuous(limits=c(0,75))+
  theme_light()



back <- joyplot_response(vis,2,182)
back$Drought <- ordered(back$Drought,levels=rev(c("ABA","AAA","ABB")))
back <- back[as.numeric(back$DAP) > 4,]
head(back)
p <- ggplot(data=back,aes(x=bin,y=factor(DAP,levels = 25:1), height=value))+
  facet_grid(~Microbes)+
  #geom_ribbon(limits,aes(fill=Drought))+
  geom_density_ridges(stat = "identity", aes(fill=Drought),alpha=0.5)+
  scale_fill_manual(values = rev(c(muted("green",c=100,l=60),muted("red",c = 100,l=60),"transparent")))+
  ylab("Days after planting")+
  xlab("")+
  xlab("NIR Grayscale Intensity")+
  theme_ridges(grid=F,center_axis_labels = T)+
  #theme(legend.position='none')+
  theme(axis.text = element_text(size = 12),
    axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
    strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p

head(vis)


makePCA <- function(data,day,start,stop){
  sub <- data[data$DAP==day,]
  channel.pca <- PCA(sub[,start:stop],graph = F)
  pca_df <- data.frame("Treatment"=sub$Drought,
    "PC1"=channel.pca$ind$coord[,1],
    "PC2"=channel.pca$ind$coord[,2])
  varexp <- signif(c(channel.pca$eig[1,2],channel.pca$eig[2,2]),4)
  
  p <- ggplot(data=pca_df, aes(PC1,PC2))+
#    geom_point(aes(color=Treatment),alpha=0.2)+
    geom_point(data=aggregate(cbind(PC1,PC2)~Treatment,pca_df,mean),aes(color=Treatment),size=5)+
    stat_ellipse(aes(fill=Treatment,color=Treatment),geom = "polygon",alpha=0.25)+
    #ylim(c(-10,10))+
    #xlim(c(-10.2,10))+
    xlab(paste("PC1 (",varexp[1],"%)",sep = ""))+
    ylab(paste("PC2 (",varexp[2],"%)",sep = ""))+
    #scale_color_manual(values = rev(c(muted("orange",l=30,c=100),muted("green",l=30,c=100),muted("purple",l=40,c=100))))+
    #scale_fill_manual(values = rev(c("orange",muted("green",l=30,c=100),"purple")))+
    geom_vline(xintercept = 0,linetype="dashed")+
    geom_hline(yintercept = 0,linetype="dashed")+
    theme_minimal()+
    theme(axis.text = element_text(size = 22),
      axis.title= element_text(size = 24))+
    theme(panel.border = element_rect(colour = "gray60", fill=NA, size=1,linetype = 1))
  p
}
p <- makePCA(vis,20,2,181)
p
p <- makePCA(nir,20,2,256)
p


hist_avg(vis,start = 2,stop = 181)

plot_histo_response <- function(data,day,which_m,mode=1){
  sub <- data[data$DAP==day & data$Microbes == which_m,]
  test_avg <- hist_avg(sub,start = 2,stop = 181)[2:71]
  test_sd <- hist_sd(sub,start = 2,stop = 181)[2:71]
  
  test_avg <- data.frame(melt(t(test_avg)))
  test_avg$sd <- data.frame(melt(t(test_sd)))[,3]
  test_avg$bin <- (2*(as.numeric(str_sub(test_avg$Var1,2,4))))
  limits <- aes(ymax = value + 1.96*sd, ymin=value - 1.96*sd)
  if(mode==1){
    p <- ggplot(data=test_avg,aes(bin,value))+
      ggtitle(paste(which_m," (Day ",day,")",sep="",collapse = ""))+
      facet_wrap(~Var2)+
      geom_ribbon(limits,fill="gray80")+
      geom_line(aes(colour=bin),size=2)+
      scale_color_gradientn(colors=hue_pal(l=65)(180)[1:70])+
      #    scale_x_continuous(breaks=seq(from=0,to=360,by=90))+
      scale_y_continuous(limits = c(-0.5,15))+
      ylab("Percentage of Mask Explained")+
      xlab("")+
      xlab("Hue Channel")+
      theme_light()+
      theme(legend.position='none')+
      theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 18))+
      theme(plot.title = element_text(hjust = 0.5),
        strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))
    return(p)
  }else{
    return(test_avg)
  }
  
}
p <- plot_histo_response(vis,day = 20,which_m = "SynCom A")
p


sub <- vis[vis$DAP==20,]
test_avg <- hist_avg(sub,start = 2,stop = 181)[2:71]
test_sd <- hist_sd(sub,start = 2,stop = 181)[2:71]
test_avg <- data.frame(melt(t(test_avg)))
test_avg$sd <- data.frame(melt(t(test_sd)))[,3]
test_avg$bin <- (2*(as.numeric(str_sub(test_avg$Var1,2,4))))
test_avg$Drought <- unlist(lapply(strsplit(as.character(test_avg$Var2),"[.]"),function(i)i[1]))
test_avg$Microbes <- unlist(lapply(strsplit(as.character(test_avg$Var2),"[.]"),function(i)i[2]))
test_avg$Microbes[test_avg$Microbes == "Mix(WW+WS)"] <- "Mix"
test_avg$Microbes <- ordered(test_avg$Microbes, levels=c("Control","WW","WS","Mix","SynCom A","SynCom B","SynCom C"))

p <- ggplot(data=test_avg,aes(x=bin,y=Microbes, height=value))+
  facet_grid(~Drought)+
  #geom_ribbon(limits,aes(fill=Drought))+
  geom_density_ridges(stat = "identity", aes(colour=Microbes),alpha=0.5,
                      position=position_points_jitter(width=0.05,height=0),
                      point_shape = "|",point_size=3)+
#  geom_rug(aes(x=1:nrow(test_avg),color=hue_pal(l=65)(180)[1:70]))+
#  scale_color_gradientn(colors=hue_pal(l=65)(180)[1:70])+
  #scale_fill_manual(values = rev(c(muted("green",c=100,l=60),muted("red",c = 100,l=60),"transparent")))+
#  ylab("Days after planting")+
  xlab("")+
  xlab("Hue Channel")+
  theme_ridges(grid=T,center_axis_labels = T)+
  theme(legend.position='none')+
  theme(axis.text = element_text(size = 12),
    axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
    strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p
