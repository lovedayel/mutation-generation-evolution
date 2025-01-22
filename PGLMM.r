library(ape)
library(tidyverse)
library(geiger)
library(phyr)
library(tibble)

#run PGLMM of log(mutation rate per year) vs. log(generation time) with phylogenetic group as a fixed effect

#read in Newick file
tree = read.tree("data/timetree_list.nwk")

#check for and resolve polytomies
resolved_tree = multi2di(tree, random = T)
is.binary(resolved_tree)

#make sure tree is ultrametric
resolved_tree = chronos(resolved_tree)
is.ultrametric(resolved_tree)

#read in data and get rid of NAs
data = read.csv("data/mutation_rates_final.csv", h = T, row.names = 1)
data_u_gen = data[,c(1,2,5,6)]
data_u_gen = na.omit(data_u_gen)

#log transform variables
data_u_gen = data_u_gen %>%
	mutate(log_u_per_year = log10(u_per_year), log_gen_time = log10(gen_time_yr))

#filter data to only include species in the tree
tree_species = resolved_tree$tip.label
filtered_u_gen = data_u_gen[rownames(data_u_gen) %in% tree_species, ]

#filter tree to only include species in the data
obj = name.check(resolved_tree, filtered_u_gen)
obj
new_tree = drop.tip(resolved_tree, obj$tree_not_data)

#change rownames to species column
filtered_u_gen = tibble::rownames_to_column(filtered_u_gen, "species")

#fit PGLMM model with phylogenetic group as fixed effect
pglmm_model = pglmm(log_u_per_year ~ log_gen_time * group + (1 | species),
			data = filtered_u_gen,
			family = "gaussian",
			cov_ranef = list(species = new_tree),
			REML = F)
summary(pglmm_model)
