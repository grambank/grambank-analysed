source("spatiophylogenetic_modelling/analysis/featurewise/INLA_results_most_prep.R")

#aestehtics choices for plotting
dat.offset = 1.02
ln.offset = 1.06
lab.offset = 1.11
fsize = 3
cex = 1.8

cols <- colour_ramp(viridis(500))

#we're reading in INLA asr objects in qs-files that were outputted by INLA_results_ASR_get.R. We'll loop over each and extract the necessary information for plotting
INLA_ASR_objects_fns <- list.files(path = "output/spatiophylogenetic_modelling/INLA_spec_results_summaries/", 
                                   pattern = "INLA_ASR_OBJ_.*qs", full.names = T)


for(fn in INLA_ASR_objects_fns){
  #fn <- INLA_ASR_objects_fns[1]
  
  INLA_ASR_object <- qs::qread(fn)
  
  
  dual_model <- INLA_ASR_object$dual_model
  pred_df <- INLA_ASR_object$pred_df
  filename <- INLA_ASR_object$filename
  filename_png <- INLA_ASR_object$filename_png
  plot_title <- INLA_ASR_object$plot_title
  feature_df <- INLA_ASR_object$feature_df
  feature <- colnames(feature_df)[2]

root <- plogis(dual_model$summary.fixed$mean[1])

tips <- pred_df$pred[1:length(tree$tip.label)]
names(tips) <- pred_df$Language_ID[1:length(tree$tip.label)]
nodes <- c(root, pred_df$pred[-1:-length(tree$tip.label)])

tiff(file = filename, width = 15.27, height = 15.69, units = "in", res = 400)

par(mar = c(1,6,1,1),
    xpd = NA) ## xpd=NA prevents clipping in the margins (just turns off clipping entirely)

p <- contMap(ladderize(tree), tips,
             anc.states = nodes,
             method = "user",
             type = "fan", ftype = "off",
             plot = FALSE)

p <- setMap(p, viridis(500))
plot(p, type = "fan", ftype = "off",
     mar = c(4, 4, 6, 4), legend = FALSE)

add.color.bar(75, p$cols, title = "Predicted Probability\n",
              prompt = FALSE, x = 0, y = 10, lwd = 12,
              fsize = 1.5)

lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

tip_points <- cbind(lastPP$xx[1:Ntip(tree)], lastPP$yy[1:Ntip(tree)],
                    feature_df[match(tree$tip.label, feature_df$Language_ID), feature])
## plot 'yes' values, blank is assumed 'no'
points(tip_points[tip_points[ , 3] == 1, 1:2] * dat.offset, col = "black", pch = 19)
## plot missing data
points(tip_points[is.na(tip_points[ , 3]), 1:2] * dat.offset, col = "grey", pch = 19)

arc.cladelabels(text="Austronesian",node =   getMRCA(tree, tip = c("kana1286", "samo1305")), mark.node=FALSE,
                ln.offset = ln.offset , lab.offset = lab.offset,fsize=fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Otomanguean",node =   getMRCA(tree, tip = c("mali1285", "yatz1235")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Uto-Aztecan",node =   getMRCA(tree, tip = c("hopi1249", "isth1240")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Nuclear Trans New Guinea",node =   getMRCA(tree, tip = c("gira1247", "domm1246")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Afro-Asiatic",node =   getMRCA(tree, tip = c("glav1244", "xamt1239")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Indo-European",node =   getMRCA(tree, tip = c("port1283", "mode1248")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Atlantic-Congo",node =   getMRCA(tree, tip = c("tswa1255", "noon1242")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Tibeto-Burman",node =   getMRCA(tree, tip = c("koir1240", "cent2004")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Austroasiatic",node =   getMRCA(tree, tip = c("seme1247", "aheu1239")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Uralic",node =   getMRCA(tree, tip = c("livv1243", "lule1254")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

title(plot_title, cex.main = cex)

dev.off()

png(file = filename_png, width = 15.27, height = 15.69, units = "in", res = 400)


par(mar = c(1,6,1,1),
    xpd = NA) ## xpd=NA prevents clipping in the margins (just turns off clipping entirely)

p <- contMap(ladderize(tree), tips,
             anc.states = nodes,
             method = "user",
             type = "fan", ftype = "off",
             plot = FALSE)

p <- setMap(p, viridis(500))
plot(p, type = "fan", ftype = "off",
     mar = c(4, 4, 6, 4),
     legend = FALSE)

add.color.bar(75, p$cols, title = "Predicted Probability\n",
              prompt = FALSE, x = 0, y = 10, lwd = 12,
              fsize = 1.5)

lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

tip_points <- cbind(lastPP$xx[1:Ntip(tree)], lastPP$yy[1:Ntip(tree)],
                    feature_df[match(tree$tip.label, feature_df$Language_ID), feature])
## plot 'yes' values, blank is assumed 'no'
points(tip_points[tip_points[ , 3] == 1, 1:2] * dat.offset, col = "black", pch = 19)
## plot missing data
points(tip_points[is.na(tip_points[ , 3]), 1:2] * dat.offset, col = "grey", pch = 19)

arc.cladelabels(text="Austronesian",node =   getMRCA(tree, tip = c("kana1286", "samo1305")), mark.node=FALSE,
                ln.offset = ln.offset , lab.offset = lab.offset,fsize=fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Otomanguean",node =   getMRCA(tree, tip = c("mali1285", "yatz1235")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Uto-Aztecan",node =   getMRCA(tree, tip = c("hopi1249", "isth1240")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Nuclear Trans New Guinea",node =   getMRCA(tree, tip = c("gira1247", "domm1246")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Afro-Asiatic",node =   getMRCA(tree, tip = c("glav1244", "xamt1239")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Indo-European",node =   getMRCA(tree, tip = c("port1283", "mode1248")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Atlantic-Congo",node =   getMRCA(tree, tip = c("tswa1255", "noon1242")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Tibeto-Burman",node =   getMRCA(tree, tip = c("koir1240", "cent2004")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Austroasiatic",node =   getMRCA(tree, tip = c("seme1247", "aheu1239")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

arc.cladelabels(text="Uralic",node =   getMRCA(tree, tip = c("livv1243", "lule1254")), mark.node=FALSE,
                ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

title(plot_title, cex.main = cex)

dev.off()
}