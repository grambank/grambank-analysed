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
interest_distances = c(250, 500, 1000, 2000, 3000)
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
paired_names = expand.grid(rownames(euclidean_dist), colnames(euclidean_dist))

plot_df = data.frame(distance = euclidean_dist[1,])
covariance_dist = sapply(spatial_parameters, function(dd) dd[1,])
colnames(covariance_dist) = paste("kappa", kappa_vec, "sigma", sigma_vec, sep = "_")

plot_df = cbind(plot_df, covariance_dist)

plot_df = plot_df[order(plot_df$distance),]

plot_df$distance_std = plot_df$distance / max(plot_df$distance)

cols = RColorBrewer::brewer.pal(5, "Set1")
plot(x = plot_df$distance,
     y = 1 - plot_df$distance_std, 
     type = "l", 
     ylab = "Covariance",
     xlab = "Haversine Distance (km)",
     ylim = c(0, 1),
     xlim = c(0, 10000)) # focus on 10,000 km
for(i in 2:6){
  lines(x = plot_df$distance, 
        y = plot_df[,i], col = cols[i-1], lwd = 2)
}
abline(v = dist_pairs$distances, lty = "dashed")
text(x = dist_pairs$distances,
     y = 0.2 + c(0, 0.1, 0.2, 0.3, 0.4),
     labels = paste0(ceiling(dist_pairs$distances), "Km"), 
     srt = 90, 
     cex = 0.75,pos = 4)

# Make legend
legend_text = c("1 - Haversine",
                paste0("Spatial: k = ", parameters[1,1],"; s = ", parameters[1,2]),
                paste0("Spatial: k = ", parameters[2,1],"; s = ", parameters[2,2]),
                paste0("Spatial: k = ", parameters[3,1],"; s = ", parameters[3,2]),
                paste0("Spatial: k = ", parameters[4,1],"; s = ", parameters[4,2]))

leg_cols = c("black", cols)
legend(7500, 1.0, 
       legend=legend_text,
       col=leg_cols, lty=1, cex=0.8, lwd = 3)
