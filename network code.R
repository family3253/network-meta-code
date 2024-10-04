# Network Meta-Analysis Script  
# Author: Yechao Chen  
# Description: This R script performs a comprehensive network meta-analysis.  
#              It includes data preparation, network construction, visualization,  
#              and statistical analysis to evaluate and compare the effects of  
#              multiple treatments. The script also encompasses model diagnostics,  
#              bias assessment, and ranking of treatments based on their performance.  
#              The analysis aims to provide insights not only into the effectiveness  
#              but also the potential biases and consistency across studies within  
#              the dataset, ultimately contributing to evidence-based decision-making   
#              in clinical settings.  

# Load necessary libraries  
library(gemtc)  
library(ggplot2)  
library(dplyr)  
library(knitr)  
library(kableExtra)  
library(igraph)  
library(tidyr)  
library(netmeta)  

# Set working directory and read dataset  
setwd("")  
data <- read.csv("xx.csv")  
colnames(data) <- c("study", "treatment", "responders", "sampleSize")  

# Create network object  
network <- mtc.network(data)  

# Generate edges for network visualization  
edges <- data %>%  
  group_by(study) %>%  
  do(data.frame(t(combn(.$treatment, 2)))) %>%  
  setNames(c("from", "to")) %>%  
  group_by(from, to) %>%  
  summarise(weight = n(), .groups = 'drop')  

# Initialize nodes  
nodes <- data.frame(name = unique(data$treatment))  

# Create and configure igraph object  
g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)  
E(g)$width <- E(g)$weight * 0.5  

# Plot network  
plot(g,  
     edge.color = "darkgrey",  
     vertex.size = data$sampleSize / 4,  
     vertex.color = "darkorange",  
     vertex.label.cex = 0.8,  
     layout = layout_in_circle,  
     main = "Network Plot")  

# Conduct network meta-analysis  
pairwise_res <- pairwise(treat = treatment, event = responders, n = sampleSize, studlab = study, data = data)  
e.netmeta <- netmeta(pairwise_res, comb.fixed = FALSE, comb.random = TRUE)  

comparison_ef <- setdiff(e.netmeta$trts, "VOR")  
ord_ef <- c(comparison_ef, "VOR")  

# Perform bias test using Egger's method  
bias_test <- metabias(e.netmeta, order = ord_ef)  
print(bias_test)  

# Set random colors for plotting  
set.seed(123)  
colors <- sample(rainbow(length(e.netmeta$comparisons)))  

# Plot funnel with Egger's regression line  
funnel(e.netmeta, order = ord_ef, col = colors, linreg = TRUE, method.bias = "Egger")  

legend("topleft", legend = comparison_ef, col = colors, pch = 19, cex = 0.75)  

# Fit random effects models  
model_consistent <- mtc.model(network, type = "consistency", n.chain = 4, likelihood = "binom", link = "log", linearModel = "random", dic = TRUE)  
results_consistent <- mtc.run(model_consistent, n.adapt = 5000, n.iter = 20000)  

model_inconsistent <- mtc.model(network, type = "use", n.chain = 4, likelihood = "binom", link = "log", linearModel = "random")  
results_inconsistent <- mtc.run(model_inconsistent, n.adapt = 5000, n.iter = 20000)  

# Produce diagnostic plots  
gelman.diag(results_consistent)  
gelman.plot(results_consistent, ask = FALSE)  

# Create and display league table  
league <- relative.effect.table(results_consistent)  
league_df <- as.data.frame(round(exp(league), 2))  

kable(league_df, caption = "Pairwise Comparison League Table", format = "pipe", align = "c", digits = 2) %>%  
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), font_size = 12)  

# Ranking plots  
ranks <- rank.probability(results_consistent, preferredDirection = 1)  
plot(ranks, beside = TRUE, cex.names = 0.5)  
plot(ranks, cumulative = TRUE)  

# Calculate SUCRA values  
sucra_values <- sucra(ranks)  
print(sucra_values)
