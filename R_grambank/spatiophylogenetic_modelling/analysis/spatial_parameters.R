# grounding covariance in distances
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

#### Data ####
languages <- read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, 
                Family_name, 
                Name, 
                Longitude, 
                Latitude, 
                Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) 

suppressWarnings(rownames(languages) <- languages$Language_ID)

phylogeny = read.tree("output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree")

# jitter points that are at exactly the same coordinates
duplicate_coords = languages[
  duplicated(languages[,c("Longitude", "Latitude")]) | 
    duplicated(languages[,c("Longitude", "Latitude")], 
               fromLast = TRUE), 
  "Language_ID"]

duplicate_rowid = languages$Language_ID %in% duplicate_coords

languages$Latitude[duplicate_rowid] = 
  jitter(languages$Latitude[duplicate_rowid], 
         factor = 2)

languages$Longitude[duplicate_rowid] = 
  jitter(languages$Longitude[duplicate_rowid], 
         factor = 2)

# reorder to match phylogeny
languages = languages[match(phylogeny$tip.label, languages$Language_ID),]

#### Distance matrices ####

## Phylogenetic distance based on branch lengths
phylo_dist = vcv.phylo(phylogeny)

## Spatial distance based on Haversine distances
euclidean_dist = distm(languages[,c("Longitude", "Latitude")],
                       fun = distHaversine) / 1000

dimnames(euclidean_dist) = list(languages$Language_ID, languages$Language_ID)

#### Make spatial matrix #####
parameters = data.frame(kappa = kappa_vec,
                        sigma =  sigma_vec)

# Make the spatial covariance matrices from the set of parameters in parameters
spatial_parameters = map2(parameters$kappa, 
                          parameters$sigma,
                          function(k, s){
                            spatial_covar_mat = 
                              varcov.spatial(
                                languages[,c("Longitude", "Latitude")], 
                                     cov.pars = c(1, s), 
                                kappa = k)$varcov
                            
                            dimnames(spatial_covar_mat) = 
                              list(languages$Language_ID,
                                   languages$Language_ID)
                            
                            spatial_covar_mat
})

## What is the covariance value for a set of known pairs
# we are interested in these distances (km)
interest_distances = c(500, 1000, 2000, 4000)
error_dist = 20
# find pairs at approximately these distances
dist_pairs = 
  sapply(interest_distances,
       function(id){
         pair = which(euclidean_dist > id - error_dist & 
                        euclidean_dist < id + error_dist, arr.ind = TRUE)[1,]
         rownames(euclidean_dist)[pair] 
       }
)

dist_pairs = data.frame(t(dist_pairs))
colnames(dist_pairs) = c("first", "second")

## Calculate those distances in KM
distances = apply(dist_pairs, 1, function(x){
  longlat1 = languages %>% dplyr::filter(Language_ID == x[1])
  longlat2 = languages %>% dplyr::filter(Language_ID == x[2])
  
  earth.dist(longlat1$Longitude, longlat1$Latitude,
             longlat2$Longitude, longlat2$Latitude)
})
dist_pairs$distances = round(distances, 2)

# Calculate the covariance of pairs at varying distances
covariances = lapply(spatial_parameters, function(s) 
  diag(s[dist_pairs$first, dist_pairs$second]))
covariances = do.call(cbind, covariances)
covariances = round(covariances, 2)
colnames(covariances) = paste("kappa", kappa_vec, "sigma", sigma_vec, sep = "_")

dist_pairs = cbind(dist_pairs, covariances)


#### Make plotting matrix ####
# rescale euclidean distances to a 0 - 1 scale
euclidean_dist = scales::rescale(euclidean_dist)

paired_names = expand.grid(rownames(euclidean_dist), colnames(euclidean_dist))

labs = paired_names[as.vector(upper.tri(euclidean_dist,diag=F)),]

plot_df = cbind(labs, 1 - euclidean_dist[upper.tri(euclidean_dist,diag=F)])
colnames(plot_df) = c("L1", "L2", "haversine_dist")

plot_df$exp_dist = scales::rescale(
  exp(
    1 - euclidean_dist[upper.tri(euclidean_dist,diag=F)]
    )
  )
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


cols = c("black", brewer.pal(7, "Set1"))

legend_text = c("1 - Haversine",
                "Phylogenetic distance",
                paste0("Spatial: k = ", parameters[1,1],"; s = ", parameters[1,2]),
                paste0("Spatial: k = ", parameters[2,1],"; s = ", parameters[2,2]),
                paste0("Spatial: k = ", parameters[3,1],"; s = ", parameters[3,2]),
                paste0("Spatial: k = ", parameters[4,1],"; s = ", parameters[4,2]))


# tiff("output/spatiophylogenetic_modelling/figures/spatial_varyingparameters.tiff", width = 8, height = 8, res = 400, units = "in")
plot(x = plot_ss$haversine_dist, y = plot_ss$haversine_dist, 
     type = "l", main = "Distance between all languages", 
     ylim = c(0, 1),
     xlab = "Haversine distance (effectively linear distance)",
     ylab = "Distance metrics (Scaled 0 - 1)"
     )
for(i in seq_along(spatialkappa_lines)){
    lines(x = plot_ss$haversine_dist, y = rev(spatialkappa_lines[[i]]), col = cols[i+2], lwd = 2)
}

abline(v = c(0.6, 0.85, 0.9, 0.95), lty = "dashed")
text(x =  c(0.6, 0.85, 0.9, 0.95) + 0.02, y = 0.5, labels = paste0(distances, " km"), srt = 90, font = 2)
legend(0.05, 1.0, 
       legend=legend_text,
       col=cols, lty=1, cex=0.8, lwd = 3)
x <- dev.off()

png("output/spatiophylogenetic_modelling/figures/spatial_varyingparameters.png", width = 8, height = 8, res = 400, units = "in")
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
