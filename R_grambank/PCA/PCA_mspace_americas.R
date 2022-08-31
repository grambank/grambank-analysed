source("requirements.R")
h_load(pkg = c("ape", "phytools", "geiger", "phangorn", "data.table", "BAMMtools", "viridis", "pals", "magick"))

#script written by Angela Chira.

##########################
########## tree & data ###
##########################

edget_fn <- "output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree"
if(!(file.exists(edget_fn))){
  source("spatiophylogenetic_modelling/processing/pruning_EDGE_tree.R")}
edget = read.tree(edget_fn)

data<-read_tsv("output/PCA/PCA_language_values.tsv", show_col_types = F)
allgb<-data

#subset data
cnames<-data[data$Language_ID%in%edget$tip.label,]$Language_ID
data<-data[data$Language_ID%in%cnames,]
x <- all(edget$tip.label%in%data$Language_ID)
data<-data[match(edget$tip.label, data$Language_ID),]
Language_meta_data <-  read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, Family_name, Name, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  mutate(Family_name = ifelse(is.na(Family_name), "Isolate", Family_name))
alllang<-Language_meta_data
Language_meta_data<-Language_meta_data[Language_meta_data$Language_ID%in%data$Language_ID,]
Language_meta_data<-Language_meta_data[match(data$Language_ID, Language_meta_data$Language_ID),]
data$Family_name<-Language_meta_data$Family_name
data$Name<-Language_meta_data$Name
data$Macroarea<-Language_meta_data$Macroarea

# all grambank
alllang<-alllang[alllang$Language_ID%in%allgb$Language_ID,]
alllang<-alllang[match(alllang$Language_ID, alllang$Language_ID),]
allgb$Macroarea<-alllang$Macroarea


############################
########## simple BM asr ###
############################

valuei<-data$PC1 ; names(valuei)<-data$Language_ID
fasr_pc1<-fastAnc(edget,valuei,vars=TRUE,CI=TRUE) # ML - can also use Bayesian

valuei<-data$PC2 ; names(valuei)<-data$Language_ID
fasr_pc2<-fastAnc(edget,valuei,vars=TRUE,CI=TRUE) # ML - can also use Bayesian

valuei<-data$PC3 ; names(valuei)<-data$Language_ID
fasr_pc3<-fastAnc(edget,valuei,vars=TRUE,CI=TRUE) # ML - can also use Bayesian

dfasr<-data.frame("Node"=names(fasr_pc1$ace),"PC1"= fasr_pc1$ace, "PC2"= fasr_pc2$ace, "PC3"= fasr_pc3$ace)
amtips<-data[data$Macroarea%in%c("South America", "North America"),]$Language_ID
is.monophyletic(edget, tips=amtips) # T
nodeam<-getMRCA(edget, tip=amtips)

dfasr$Nodeam<-0
dfasr[dfasr$Node%in%getDescendants(edget,nodeam),]$Nodeam<-1
table(dfasr$Nodeam)


###############################
########## pc1 + pc2 ~ time ###
###############################

source("PCA_mspace_americas_fxns.R")

ages <- nodeAges(edget, min.age = 0)
ages.nodes <- ages$nodeage
ages.tips <- ages$tip.age
age.root <- max(ages.nodes) # age of root

x <- identical(names(ages.nodes), dfasr$Node)
dfasr$NodeAge<-ages.nodes


# logged time intervals
range(dfasr$NodeAge)
age.intervals_linear<-seq(from = min(dfasr$NodeAge), to =  max(dfasr$NodeAg)+1e-12, by=1)
age.intervals_log <- assignColorBreaks(age.intervals_linear, NCOLORS=length(age.intervals_linear)-1, logcolor=T, method="linear")

# check
hist(dfasr$NodeAge, xlab="", main="", col="gray90", border = "gray90")
abline(v=age.intervals_log, col=rev(vector_col(age.intervals_linear)$cols))
# plot(age.intervals_linear~ age.intervals_linear, col=rev(vector_col(age.intervals_linear)$cols))

# PLOT -- 1 side: all, 1 side: Americas, same time scale
x <- dev.off()

# par(mfrow=c(1,2))

# pdf("output/PCA/mspace_americas.pdf", width=10, height=8)
png("output/PCA/mspace_americas.png",  width=10, height=8, res=480, units='in')

matlay<-matrix(c(1,3,2,3),nrow=2)
nf<-layout(matlay, widths = c(2,2,2,2), heights =c(2,1,1,1),TRUE)
# layout.show(nf)
par(mar=c(5,5,2,3))

plot(0, bty = 'n', pch = '', ylab = "PC2", xlab = "PC1", xlim= range(allgb$PC1), ylim = range(allgb$PC2), xaxt='n', yaxt='n',font.lab=2, cex.lab=1.3)
mtext(text="a", side=3,adj=0, line=+0.3,cex=1.3)
axis(1, lwd=2, tck=-0.01, cex.axis=1.1)
axis(2, lwd=2, tck=-0.01, cex.axis=1.1) #lwd.tick=0
par(xpd=F)
#abline(v=0, lty=2, col="gray90")
#abline(h=0, lty=2, col="gray90")

col_plot<-vector_col(age.intervals_linear)$cols
col_plot<-rev(col_plot)
age.intervals<-age.intervals_log
if(max(age.intervals)< max(dfasr$NodeAge)) age.intervals[age.intervals%in%max(age.intervals)]<-max(dfasr$NodeAge)+1e-12

# tips
# points(data$PC2 ~ data$PC1, col=scales::alpha(col_plot[1],1), pch=20)
# points(data[!data$Macroarea%in%c("South America", "North America"),]$PC2 ~ data[!data$Macroarea%in%c("South America", "North America"),]$PC1, 
#        col=scales::alpha(col_plot[1],1), pch=20)

points(allgb[!allgb$Macroarea%in%c("South America", "North America"),]$PC2 ~ allgb[!allgb$Macroarea%in%c("South America", "North America"),]$PC1, 
       col=scales::alpha(col_plot[1],1), pch=20) # all languages in GB for present slice

for(i in 1: (length(age.intervals)-1)) {
  
  cmin<-age.intervals[i]
  cmax<-age.intervals[i+1]
  
  # each set of points once - does without rev(age.intervals)
  pPC2<-dfasr[dfasr$NodeAge>=cmin & dfasr$NodeAge<cmax & dfasr$Nodeam==0,]$PC2
  pPC1<-dfasr[dfasr$NodeAge>=cmin & dfasr$NodeAge<cmax & dfasr$Nodeam==0,]$PC1

  colsp<-col_plot[i]

  points( pPC2 ~ pPC1,pch=20, col=scales::alpha(colsp,0.9),cex=1.2)
  paste0("n = ",length(pPC1))
  
  print(i)  
}


# shine on Americas
plot(0, bty = 'n', pch = '', ylab = "PC2", xlab = "PC1", xlim= range(allgb$PC1), ylim = range(allgb$PC2), xaxt='n', yaxt='n',font.lab=2, cex.lab=1.3)
mtext(text="b", side=3,adj=0, line=+0.3, cex=1.3)
axis(1, lwd=2, tck=-0.01, cex.axis=1.1)
axis(2, lwd=2, tck=-0.01, cex.axis=1.1) #lwd.tick=0
par(xpd=F)
#abline(v=0, lty=2, col="gray90")
#abline(h=0, lty=2, col="gray90")

col_plot<-vector_col(age.intervals_linear)$cols
col_plot<-rev(col_plot)
age.intervals<-age.intervals_log
if(max(age.intervals)< max(dfasr$NodeAge)) age.intervals[age.intervals%in%max(age.intervals)]<-max(dfasr$NodeAge)+1e-12

# everything else
# yp<-c(dfasr[dfasr$Nodeam==0,]$PC2, data[!data$Macroarea%in%c("South America", "North America"),]$PC2)
# xp<-c(dfasr[dfasr$Nodeam==0,]$PC1, data[!data$Macroarea%in%c("South America", "North America"),]$PC1)
yp<-c(dfasr[dfasr$Nodeam==0,]$PC2, allgb[!allgb$Macroarea%in%c("South America", "North America"),]$PC2)
xp<-c(dfasr[dfasr$Nodeam==0,]$PC1, allgb[!allgb$Macroarea%in%c("South America", "North America"),]$PC1)
points(yp ~ xp, col=scales::alpha("gray94",0.7), pch=20) 

# tips in Americas
# points(data[data$Macroarea%in%c("South America", "North America"),]$PC2 ~ data[data$Macroarea%in%c("South America", "North America"),]$PC1, 
#        col=scales::alpha(col_plot[1],1), pch=20)
points(allgb[allgb$Macroarea%in%c("South America", "North America"),]$PC2 ~ allgb[allgb$Macroarea%in%c("South America", "North America"),]$PC1, 
       col=scales::alpha(col_plot[1],1), pch=20) # all gb

for(i in 1: (length(age.intervals)-1)) {
  
  cmin<-age.intervals[i]
  cmax<-age.intervals[i+1]
  
  # each set of points once - does without rev(age.intervals)
  if(any(dfasr$NodeAge>=cmin & dfasr$NodeAge<cmax &dfasr$Nodeam==1)){
    pPC2<-dfasr[dfasr$NodeAge>=cmin & dfasr$NodeAge<cmax & dfasr$Nodeam==1,]$PC2
    pPC1<-dfasr[dfasr$NodeAge>=cmin & dfasr$NodeAge<cmax & dfasr$Nodeam==1,]$PC1
    
    colsp<-col_plot[i]
    # colsp<-dfasr[dfasr$NodeAge>=cmin & dfasr$NodeAge<cmax,]$NodeCol
    
    points( pPC2 ~ pPC1,pch=20, col=scales::alpha(colsp,0.9),cex=1.2)
    paste0("n = ",length(pPC1))
  }
  print(i)  
}


# add legend
par(xpd=T)
hist(dfasr$NodeAge, col="gray50", border = "gray50", xaxt='n', yaxt='n', xlim=c(0,120), # sorry for the manual xlim fix!
     xlab="", main="", ylab="", breaks=200)
clip(0, max(age.intervals_log), -50, -10)
abline(v=age.intervals_log, col=rev(vector_col(age.intervals_linear)$cols), lwd=2.5)
# abline(v=age.intervals_log, col=vector_col(age.intervals_linear)$cols, lwd=2.5)
mtext(text="c", side=3,adj=0, line=+0.3,cex=1.3)
axis(1, lwd=2, tck=-0.01, line=+1.5)
text(x=-8,y=-110, labels="Time (kya)", cex=1.2)

# add Americas arrow
xarrow<-dfasr[dfasr$Node%in%getMRCA(edget, tip=amtips),]$NodeAge
arrows(x0 = xarrow, y0=-170, y1 = -90, lwd=2,length = 0.10)
text(x=xarrow+12, y=-155, labels="root of the Americas", cex=1.2) # sorry for the manual placement!
dev.off()


##################
# ## GIF IT  ##### only what it is in the tree
##################

col_plot<-vector_col(age.intervals_linear)$cols
col_plot<-rev(col_plot)
age.intervals<-age.intervals_log
if(max(age.intervals)< max(dfasr$NodeAge)) age.intervals[age.intervals%in%max(age.intervals)]<-max(dfasr$NodeAge)+1e-12

# dfasr$NodeCol<-vector_col((-1)*dfasr$NodeAge)$cols
if(!"mspaceamericasgif"%in%list.files("output/PCA/")) dir.create("output/PCA/mspaceamericasgif")

i=length(age.intervals)-1
j=1
for(i in (length(age.intervals)-1):1 ) {
  
  cmin<-age.intervals[i]
  
  pPC2<-dfasr[dfasr$NodeAge>=cmin,]$PC2
  pPC1<-dfasr[dfasr$NodeAge>=cmin,]$PC1
  
  #colsp<-col_plot[i]
  colsp<-dfasr[dfasr$NodeAge>=cmin,]$Nodeam
  colsp[colsp==0]<-"gray90"
  colsp[colsp==1]<-"darkmagenta"
  
  png(paste0("output/PCA/mspaceamericasgif/image_",j, ".png"), width=8, height=6, res=300, units='in')
  # pdf(paste0("output/PCA/mspaceamericasgif/image_",j, ".pdf"), height=7, width=9)
  
  par(mar=c(10,8,3,3))
  plot(0, bty = 'n', pch = '', ylab = "PC2", xlab = "PC1", xlim= range(data$PC1), ylim = range(data$PC2), xaxt='n', yaxt='n',font.lab=2, cex.lab=1.3)
  axis(1, lwd=2, tck=-0.01, cex.axis=1.1)
  axis(2, lwd=2, tck=-0.01, cex.axis=1.1) #lwd.tick=0
  par(xpd=F)
  
  if(!"darkmagenta"%in%colsp) points( pPC2 ~ pPC1,pch=20, col=scales::alpha(colsp,0.5),cex=1.2)
  if("darkmagenta"%in%colsp){
    wam<-which(colsp=="darkmagenta") 
    points( pPC2[-wam] ~ pPC1[-wam],pch=20, col=scales::alpha(colsp,0.5)[-wam],cex=1.2)
    points( pPC2[wam] ~ pPC1[wam],pch=20, col=scales::alpha(colsp,0.5)[wam],cex=1.2)
  }
  
  paste0("n = ",length(pPC1))
  par(xpd=T)
  text(x=max(data$PC1)-2, y=max(data$PC2), labels=paste0("n = ",length(pPC1)), font = 2)
  
  # time bar & arrow
  par(new = TRUE)
  par(mar=c(10,8,3,3))
  par(xpd=T)
  xp<-(-1)*age.intervals
  plot(age.intervals~xp, bty = 'n', pch = '', ylab = "", xlab = "", xaxt='n', yaxt='n', xlim=c(-120,0))
  axis(1, lwd=2, tck=-0.01, line=+6)
  # arrows(x0 = (-1)*ceiling(max(age.intervals)), x1=(-1)*cmin, y0 =-40, lwd=2,length = 0.10)
  text(x=-134,y=-59, "Time (kya)")
  
  clip(-120, 0, -45, -40)
  xplog<-(-1)*age.intervals_log
  
  cs<-rev(vector_col(age.intervals_linear)$cols)
  abline(v=xplog[xplog> (-1)*cmin], lwd=2.5,col= "white") #col="#39897e"
  # abline(v=xplog[xplog<= (-1)*cmin], lwd=2.5, col=cs[xplog<= (-1)*cmin] )  #col="#E7242D"
  abline(v=xplog[xplog<= (-1)*cmin], lwd=2.5, col="#E7242D" )  #col="#E7242D"
  
  # root arrow
  # arrows(x0 = (-1)*cmin, y0=-55, y1 = -40, lwd=2,length = 0.10, col="red")
  dev.off()
  
  print(i)
  j<-j+1
  
}

# # tips at last i
png(paste0("output/PCA/mspaceamericasgif/image_",j, ".png"), width=8, height=6, res=300, units='in')
# pdf(paste0("output/PCA/mspaceamericasgif/image_",j, ".pdf"), height=7, width=9)

par(mar=c(10,8,3,3))
plot(0, bty = 'n', pch = '', ylab = "PC2", xlab = "PC1", xlim= range(data$PC1), ylim = range(data$PC2), xaxt='n', yaxt='n',font.lab=2, cex.lab=1.3)
axis(1, lwd=2, tck=-0.01, cex.axis=1.1)
axis(2, lwd=2, tck=-0.01, cex.axis=1.1) #lwd.tick=0
par(xpd=F)
yp<-c(dfasr$PC2, data$PC2)
xp<-c(dfasr$PC1, data$PC1)
colsp1<-dfasr$Nodeam
colsp1[colsp1==0]<-"gray90"
colsp1[colsp1==1]<-"darkmagenta"
colsp<-data$Macroarea
colsp[!colsp%in%c("South America", "North America")]<-"gray90"
colsp[colsp%in%c("South America", "North America")]<-"darkmagenta"
colsp<-c(colsp1,colsp)

if(!"darkmagenta"%in%colsp) points( yp ~ xp,pch=20, col=scales::alpha(colsp,0.5),cex=1.2)
if("darkmagenta"%in%colsp){
  wam<-which(colsp=="darkmagenta")
  points( yp[-wam] ~ xp[-wam],pch=20, col=scales::alpha(colsp,0.5)[-wam],cex=1.2)
  points( yp[wam] ~ xp[wam],pch=20, col=scales::alpha(colsp,0.5)[wam],cex=1.2)
}

paste0("n = ",length(pPC1))
par(xpd=T)
text(x=max(data$PC1)-2, y=max(data$PC2), labels=paste0("n = ",length(yp)), font = 2)
# time bar & arrow
par(new = TRUE)
par(mar=c(10,8,3,3))
par(xpd=T)
xp<-(-1)*age.intervals
plot(age.intervals~xp, bty = 'n', pch = '', ylab = "", xlab = "", xaxt='n', yaxt='n', xlim=c(-120,0))
axis(1, lwd=2, tck=-0.01, line=+6)
# arrows(x0 = (-1)*ceiling(max(age.intervals)), x1=(-1)*cmin, y0 =-40, lwd=2,length = 0.10)
text(x=-134,y=-59, "Time (kya)")
clip(-120, 0, -45, -40)
xplog<-(-1)*age.intervals_log
cs<-rev(vector_col(age.intervals_linear)$cols)
abline(v=xplog[xplog> (-1)*cmin], lwd=2.5,col="white" ) #col="#39897e"
# abline(v=xplog[xplog<= (-1)*cmin], lwd=2.5, col=cs[xplog<= (-1)*cmin] )  #col="#E7242D"
abline(v=xplog[xplog<= (-1)*cmin], lwd=2.5, col="#E7242D" )  #col="#E7242D"
dev.off()

# create gif
# https://www.nagraj.net/notes/gifs-in-r/
imgs <- list.files("output/PCA/mspaceamericasgif", full.names = T) ## list file names and read in
# imgs<-imgs[grep(".pdf", imgs)]
ordern<-seq(from =1 , to=length(imgs)) # order this by number not alphabetically
newimgs<-character()
for(i in 1:length(imgs)){
  imgsi<-strsplit(imgs[i], split="_")[[1]][2]
  imgsi<-strsplit(imgsi, split=".png")[[1]][1]
  #imgsi<-strsplit(imgsi, split=".pdf")[[1]][1]
  imgsi<-as.numeric(imgsi)
  newimgs[imgsi]<-imgs[i]
}

imgs<-newimgs
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list) ## join the images together
img_animated <- image_animate(img_joined, fps = 2) ## animate at 2 frames per second
# img_animated # view animated image
image_write(image = img_animated, path = "output/PCA/mspace_americas.gif") ## save to disk

