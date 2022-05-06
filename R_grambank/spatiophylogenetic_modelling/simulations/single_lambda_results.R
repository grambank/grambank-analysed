# single lambda results
library(ggplot2)
library(caper)

output = readRDS('lambda0.7_simulation.RDS')
tree = read.tree('data/jaeger_pruned.tree')

get_lambda_inla = function(fit, effect){
  hyper_summary = fit$summary.hyperpar
  eff_1 = hyper_summary[effect, "mean"]
  eff_2 = sum(hyper_summary[,"mean"])
  binomial_constant = pi^2/3
  (1 / eff_1) / 
    ((1/eff_1) + (1/eff_2) + (1/binomial_constant))
}

get_lambda_brms = function(fit, effect){
  fff = summary(fit)
  estimates = sapply(fff$random, function(x) x$Estimate)
  binomial_constant = (pi^3)/ 3
  estimates[effect]^2  / (sum(estimates^2) + binomial_constant)
}

pagels_lambda_results = lapply(output, function(x) x$pagels_lambda)
inla_results = lapply(output, function(x) x$inla_model)
brms_results = lapply(output, function(x) x$brms_model)

# Estimates 
pagels_lambda_estimate = sapply(pagels_lambda_results, function(x) x$opt$lambda)
inla_estimate = sapply(inla_results, get_lambda_inla, effect = "Precision for phy_id_int")
brms_estimate = sapply(brms_results, get_lambda_brms, effect = "glottocodes")

## Phylo.D test
ys = sapply(output, function(x) x$y)
taxa = rownames(ys)
sum(taxa %in% tree$tip.label)
ys = apply(ys, 2, as.numeric)
ys = as.data.frame(ys)
ys$taxa = taxa

d_stat = list()
for(i in 1:ncol(ys)){
  var = colnames(ys)[i]
  
  print(
    table(ys[,i])
  )
  
  d_stat[[i]] = eval(
            substitute(
              phylo.d(data = ys,
                      names.col = taxa, 
                      phy = tree,
                      binvar = this_var), 
              list(this_var=as.name(var)))) 
}


d_estimate = sapply(
  d_stat, function(x) x$DEstimate
)

plot_df = data.frame(pagels_lambda = pagels_lambda_estimate,
                     inla = inla_estimate,
                     brms = brms_estimate,
                     phylo_d = d_estimate
                     )
plot_df = plot_df[order(plot_df$pagels_lambda),]
present = sapply(output, function(x) table(x$y)["1"])

plot_long = reshape2::melt(plot_df)
plot_long$id = rep(1:20, times = 4)
plot_long$present = rep(present, each = 4)

ggplot(plot_long, aes(y = value, x = variable, fill = variable)) + 
  geom_boxplot() + 
  geom_jitter(height = 0, width = 0.05, aes(col = id, size = present)) + 
  geom_hline(yintercept = 0.7, lty = "dashed", col = "red") + 
  ylab("Lambda estimate") + 
  xlab("Estimator") + 
  theme(legend.position = "bottom",
        legend.title = element_blank())

ggplot(plot_df, aes(y = pagels_lambda, x = inla)) + 
  geom_point()

ggplot(plot_df, aes(y = pagels_lambda, x = brms)) + 
  geom_point()

ggplot(plot_df, aes(y = inla, x = brms, size = present)) + 
  geom_point()

ggplot(plot_df, aes(y = inla, x = present, size = present)) + 
  geom_point()
