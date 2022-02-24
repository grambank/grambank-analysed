# grounding distance
#

#script written by Sam Passmore
source('requirements.R')

# load variational covariance matrix taken from geoR::varcov_spatial
source('spatiophylogenetic_modelling/analysis/varcov_spatial.R')

#### Functions #### 
#### # Calculate Haversine distance in kilometers between two points (as the crow flies taking into account curvature of the earth).
earth.dist <- function (long1, lat1, long2, lat2)
{
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}

plotCircle <- function(LonDec, LatDec, Km, col = "black") {#Corrected function
  #LatDec = latitude in decimal degrees of the center of the circle
  #LonDec = longitude in decimal degrees
  #Km = radius of the circle in kilometers
  ER <- 6371 #Mean Earth radius in kilometers. Change this to 3959 and you will have your function working in miles.
  AngDeg <- seq(1:360) #angles in degrees 
  Lat1Rad <- LatDec*(pi/180)#Latitude of the center of the circle in radians
  Lon1Rad <- LonDec*(pi/180)#Longitude of the center of the circle in radians
  AngRad <- AngDeg*(pi/180)#angles in radians
  Lat2Rad <-asin(sin(Lat1Rad)*cos(Km/ER)+cos(Lat1Rad)*sin(Km/ER)*cos(AngRad)) #Latitude of each point of the circle rearding to angle in radians
  Lon2Rad <- Lon1Rad+atan2(sin(AngRad)*sin(Km/ER)*cos(Lat1Rad),cos(Km/ER)-sin(Lat1Rad)*sin(Lat2Rad))#Longitude of each point of the circle rearding to angle in radians
  Lat2Deg <- Lat2Rad*(180/pi)#Latitude of each point of the circle rearding to angle in degrees (conversion of radians to degrees deg = rad*(180/pi) )
  Lon2Deg <- Lon2Rad*(180/pi)#Longitude of each point of the circle rearding to angle in degrees (conversion of radians to degrees deg = rad*(180/pi) )
  polygon(Lon2Deg,Lat2Deg,lty=1, lwd = 3, border = col)
}

#### Data ####
languages <- read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, Family_name, Name, Longitude, Latitude, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) 

suppressWarnings(rownames(languages) <- languages$Language_ID)

phylogeny = read.tree("spatiophylogenetic_modelling/processed_data/jaeger_pruned.tree")

#jitter points that are at exactly the same coordinates
duplicate_coords = languages[duplicated(languages[,c("Longitude", "Latitude")]) | duplicated(languages[,c("Longitude", "Latitude")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = languages$Language_ID %in% duplicate_coords
languages$Latitude[duplicate_rowid] = jitter(languages$Latitude[duplicate_rowid], factor = 2)
languages$Longitude[duplicate_rowid] = jitter(languages$Longitude[duplicate_rowid], factor = 2)

# reorder to match phylogeny
languages = languages[match(phylogeny$tip.label, languages$Language_ID),]

#### Distance matrices ####
## Phylogenetic
phylo_dist = vcv.phylo(phylogeny)

euclidean_dist = distm(languages[,c("Longitude", "Latitude")],
                       fun = distHaversine)
dimnames(euclidean_dist) = list(languages$Language_ID, languages$Language_ID)
# scale 
euclidean_dist = scales::rescale(euclidean_dist)

#### Make spatial matrix #####
parameters = data.frame(kappa = kappa_vec,
                        sigma =  sigma_vec)

spatial_parameters = map2(parameters$kappa, parameters$sigma,function(k, s){
  spatial_covar_mat = varcov.spatial(languages[,c("Longitude", "Latitude")], cov.pars = c(1, s), kappa = k)$varcov
  dimnames(spatial_covar_mat) = list(languages$Language_ID, languages$Language_ID)
  spatial_covar_mat
})

#### Make plotting matrix ####
# all_dist = c(1 - euclidean_dist[lower.tri(euclidean_dist)])
paired_names = expand.grid(rownames(euclidean_dist), colnames(euclidean_dist))

labs = paired_names[as.vector(upper.tri(euclidean_dist,diag=F)),]
plot_df = cbind(labs, 1 - euclidean_dist[upper.tri(euclidean_dist,diag=F)])
colnames(plot_df) = c("L1", "L2", "haversine_dist")
plot_df$exp_dist = scales::rescale(exp(1 - euclidean_dist[upper.tri(euclidean_dist,diag=F)]))
plot_df$phylo_dist = phylo_dist[upper.tri(phylo_dist,diag=F)]

## subset to a plottable set
plot_n = 1000

sample_idx = ceiling(seq(1, nrow(plot_df)-1, length.out = plot_n))
plot_ss = plot_df[sample_idx,]
plot_ss$index = sample_idx
plot_ss = plot_ss[order(plot_ss$haversine_dist),]

spatialkappa_lines = lapply(spatial_parameters, function(x) {
  d = sort(c(x[lower.tri(x)]), decreasing = TRUE)
  sample_idx = seq(1, length(d), length.out = plot_n)  
  d[sample_idx]
})

distpair_names = data.frame(first = 
                              c("gugu1253", "yami1254", "awar1249", "auuu1241"), 
                            second = 
                              c("juan1238", "kwin1241", "taul1251", "urim1252"))

## Distances values
distances = apply(distpair_names, 1, function(x){
  longlat1 = languages %>% dplyr::filter(Language_ID == x[1])
  longlat2 = languages %>% dplyr::filter(Language_ID == x[2])

  earth.dist(longlat1$Longitude, longlat1$Latitude,
             longlat2$Longitude, longlat2$Latitude)
  })

distances = round(distances, 2)

cols = c("black", brewer.pal(7, "Set1"))

legend_text = c("1 - Haversine",
                "Phylogenetic distance",
                paste0("Spatial: k = ", parameters[1,1],"; s = ", parameters[1,2]),
                paste0("Spatial: k = ", parameters[2,1],"; s = ", parameters[2,2]),
                paste0("Spatial: k = ", parameters[3,1],"; s = ", parameters[3,2]),
                paste0("Spatial: k = ", parameters[4,1],"; s = ", parameters[4,2]))


tiff("spatiophylogenetic_modelling/figures/spatial_varyingparameters.tiff", width = 8, height = 8, res = 400, units = "in")
plot(x = plot_ss$haversine_dist, y = plot_ss$haversine_dist, 
     type = "l", main = "Distance between all languages", 
     ylim = c(0, 1),
     xlab = "Haversine distance (effectively linear distance)",
     ylab = "Distance metrics (Scaled 0 - 1)"
     )
for(i in seq_along(spatialkappa_lines)){
    lines(x = plot_ss$haversine_dist, y = rev(spatialkappa_lines[[i]]), col = cols[i+2], lwd = 2)
}

lines(x = plot_ss$haversine_dist, y =  sort(plot_ss$phylo_dist), col = cols[2], lwd = 2)
abline(v = c(0.6, 0.85, 0.9, 0.95), lty = "dashed")
text(x =  c(0.6, 0.85, 0.9, 0.95) + 0.02, y = 0.5, labels = paste0(distances, " km"), srt = 90, font = 2)
legend(0.05, 1.0, 
       legend=legend_text,
       col=cols, lty=1, cex=0.8, lwd = 3)
x <- dev.off()

png("spatiophylogenetic_modelling/figures/spatial_varyingparameters.png", width = 8, height = 8, res = 400, units = "in")
plot(x = plot_ss$haversine_dist, y = plot_ss$haversine_dist, 
     type = "l", main = "Distance between all languages", 
     ylim = c(0, 1),
     xlab = "Haversine distance (effectively linear distance)",
     ylab = "Distance metrics (Scaled 0 - 1)"
)
for(i in seq_along(spatialkappa_lines)){
  lines(x = plot_ss$haversine_dist, y = rev(spatialkappa_lines[[i]]), col = cols[i+2], lwd = 2)
}

lines(x = plot_ss$haversine_dist, y =  sort(plot_ss$phylo_dist), col = cols[2], lwd = 2)
abline(v = c(0.6, 0.85, 0.9, 0.95), lty = "dashed")
text(x =  c(0.6, 0.85, 0.9, 0.95) + 0.02, y = 0.5, labels = paste0(distances, " km"), srt = 90, font = 2)
legend(0.05, 1.0, 
       legend=legend_text,
       col=cols, lty=1, cex=0.8, lwd = 3)
x <- dev.off()



