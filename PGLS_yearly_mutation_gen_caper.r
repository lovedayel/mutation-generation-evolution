library(ggplot2)
library(ggpubr)
library(ape)
library(caper)
library(dplyr)

#mutation rate per year vs. generation time with PGLS correction

#read in Newick file
tree = read.tree("data/timetree_list.nwk")

#check for and resolve polytomies
resolved_tree = multi2di(tree, random = TRUE)
is.binary(resolved_tree)

#make sure tree is ultrametric
resolved_tree = chronos(resolved_tree)
is.ultrametric(resolved_tree)

#read in data
data = read.csv("data/mutation_rates_final.csv", h = TRUE, row.names = 1)
data_u_gen = data[, c(1,5,6)]
data_u_gen = na.omit(data_u_gen)

#filter data to only include species in the tree
tree_species = resolved_tree$tip.label
filtered_u_gen = data_u_gen[rownames(data_u_gen) %in% tree_species, ]

obj = name.check(resolved_tree, filtered_u_gen)
obj
new_tree = drop.tip(resolved_tree, obj$tree_not_data)

#ensure species names are in a column
filtered_u_gen$species <- rownames(filtered_u_gen)

#create comparative data object for caper
comp_data = comparative.data(phy = new_tree, data = filtered_u_gen, names.col = "species", vcv = TRUE, warn.dropped = TRUE)

#log transformations
comp_data$data$log_u = log10(comp_data$data$u_per_year)
comp_data$data$log_gen_time = log10(comp_data$data$gen_time_yr)

#fit PGLS model
pgls_model = pgls(log_u ~ log_gen_time, data = comp_data)
summary(pgls_model)

#get PGLS coefficient and p-value
model_summary = summary(pgls_model)
pgls_coef_slope = round(model_summary$coefficients[2, "Estimate"], 2)
pgls_coef_intercept = round(model_summary$coefficients[1, "Estimate"], 2)
pgls_raw_p = model_summary$coefficients[2, "Pr(>|t|)"]
pgls_p = ifelse(pgls_raw_p < 2.2e-16, "p < 2.2e-16", formatC(pgls_raw_p, format = "e", digits = 2))
print(pgls_p)

#create regression equation string
regression_eq = paste("y = ", pgls_coef_slope, "x", " + ", pgls_coef_intercept)

#extract observed and fitted values
observed_values = log10(comp_data$data$u_per_year)
fitted_values = fitted(pgls_model)

#make dataframe for plotting
plot_data = data.frame(log_gen_time_yr = log10(comp_data$data$gen_time_yr), observed_values = observed_values, fitted_values = fitted_values, group = comp_data$data$group)

#plot
pdf("u_gen_yearly_pgls_caper.pdf", width = 9, height = 5)
fig1 = ggplot(plot_data, aes(x = log_gen_time_yr, y = observed_values)) +
    geom_point(aes(colour = group)) +
    geom_line(aes(y = fitted_values), colour = "black") +
    labs(x = "Log10 Generation Time (years)", y = "Log10 Mutation Rate (per year)") +
    guides(colour = guide_legend(title = NULL)) +
    theme_minimal() +
    annotate("text", x = max(plot_data$log_gen_time_yr), y = max(plot_data$observed_values), 
             label = paste("PGLS Coef =", pgls_coef_slope, ", ", pgls_p), 
             hjust = 0.89, vjust = 2.3, size = 4) +
    annotate("text", x = min(plot_data$log_gen_time_yr), y = max(plot_data$observed_values - 0.2), label = regression_eq, hjust = -3.68, vjust = 2.3, size = 4)
print(fig1)
dev.off()