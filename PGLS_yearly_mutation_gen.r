library(ggplot2)
library(ggpubr)
library(ape)
library(nlme)
library(geiger)
library(dplyr)
library(rr2)

#mutation rate per year vs. generation time with PGLS correction

#read in Newick file
tree = read.tree("data/timetree_list.nwk")

#check for and resolve polytomies
resolved_tree = multi2di(tree, random = T)
is.binary(resolved_tree)

#make sure tree is ultrametric
resolved_tree = chronos(resolved_tree)
is.ultrametric(resolved_tree)

#read in data
data = read.csv("data/mutation_rates_final.csv", h = T, row.names = 1)
data_u_gen = na.omit(data)

#filter data to only include species in the tree
tree_species = resolved_tree$tip.label
filtered_u_gen = data_u_gen[rownames(data_u_gen) %in% tree_species, ]

obj = name.check(resolved_tree, filtered_u_gen)
obj
new_tree = drop.tip(resolved_tree, obj$tree_not_data)

#define covariance structure
spp = rownames(filtered_u_gen)
corBM = corBrownian(phy = new_tree, form = ~spp)
corBM

#log transformations
filtered_u_gen$log_u = log10(filtered_u_gen$u_per_year)
filtered_u_gen$log_gen_time = log10(filtered_u_gen$gen_time_yr)

#fit GLS model
pgls = gls(log_u ~ log_gen_time, data = filtered_u_gen, correlation = corBM)
summary(pgls)

#get PGLS coefficient and p-value
model_summary = summary(pgls)
pgls_coef_slope = round(model_summary$tTable[2, "Value"], 2)
pgls_coef_intercept = round(model_summary$tTable[1, "Value"], 2)
pgls_p = round(model_summary$tTable[2, "p-value"], 22)

#create regression equation string
regression_eq = paste("y = ", pgls_coef_slope, "x", " + ", pgls_coef_intercept)

#extract observed and expected values
observed_values = log10(filtered_u_gen$u_per_year)
fitted_values = fitted(pgls)

#make dataframe for plotting
plot_data = data.frame(log_gen_time_yr = log10(filtered_u_gen$gen_time_yr), observed_values = observed_values, fitted_values = fitted_values, group = filtered_u_gen$group)

#plot
pdf("u_gen_yearly_pgls2.pdf", width = 9, height = 5)
fig1 = ggplot(plot_data, aes(x = log_gen_time_yr, y = observed_values)) +
	geom_point(aes(colour = group)) +
	geom_line(aes(y = fitted_values), colour = "black") +
	labs(x = "Log10 Generation Time (years)", y = "Log10 Mutation Rate (per year)") +
	guides(colour = guide_legend(title = NULL)) +
	theme_minimal() +
	annotate("text", x = max(plot_data$log_gen_time_yr), y = max(plot_data$observed_values), label = paste("PGLS Coef =", pgls_coef_slope, ", p =", pgls_p), hjust = 1, vjust = 2.3, size = 4) +
    annotate("text", x = min(plot_data$log_gen_time_yr), y = max(plot_data$observed_values - 0.2), label = regression_eq, hjust = -3.78, vjust = 2.3, size = 4)
print(fig1)
dev.off()

r_squared = R2(mod = pgls, phy = T)
print(r_squared)
