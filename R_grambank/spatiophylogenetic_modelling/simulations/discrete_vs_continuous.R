## Continious trait test
library(ape)
library(geiger)
library(ggplot2)

tree = read.tree('./data/jaeger_pruned.tree')
lambda = 0.7
N = 50

rescaled_tree = geiger::rescale(tree,
                  lambda,
                  model = "lambda")

disc = c()
cont = c()
for(i in 1:N){
  print(i)
  y_disc = rTraitDisc(
    rescaled_tree,
    k = 2,
    freq = 0.5,
    states = 0:1
  )
  
  y_cont = rTraitCont(
    phy = rescaled_tree, 
    model = "BM", 
  )
  
  pagel_disc = fitDiscrete(tree, 
                           factor(y_disc), 
                           transform = "lambda")
  
  pagel_cont = fitContinuous(tree, 
                             y_cont, 
                             model = "lambda")
  
  disc[i] = pagel_disc$opt$lambda
  cont[i] = pagel_cont$opt$lambda
}

plot_df = data.frame(continuous = cont,
                     discrete = disc)
plot_long = reshape::melt(plot_df)
plot_long$id = rep(1:N, times = 2)

ggplot(plot_long, 
       aes(y = value, x = variable, fill = variable)) + 
  geom_boxplot() + 
  geom_jitter(height = 0, width = 0.05, aes(col = id)) + 
  geom_hline(yintercept = 0.7, lty = "dashed", col = "red") + 
  ylab("Lambda estimate") + 
  xlab("Estimator") + 
  theme(legend.position = "bottom",
        legend.title = element_blank())
