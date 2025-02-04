library(ggplot2)
library(gridExtra)
library(cowplot)
library(ape)
library(caper)
library(dplyr)
library(rr2)

#read in Newick file
tree = read.tree("data/timetree_list.nwk")

#resolve polytomies
resolved_tree = multi2di(tree, random = TRUE)

#ensure tree is ultrametric
resolved_tree = chronos(resolved_tree)

#read in data
data = read.csv("data/mutation_rates_final.csv", header = TRUE, row.names = 1)
data$group[data$group == "" | data$group == " "] <- NA
data = data[, c("u_per_year", "gen_time_yr", "group")]
data_u_gen = na.omit(data)

#filter data to only include species in the tree
tree_species = resolved_tree$tip.label
filtered_u_gen = data_u_gen[rownames(data_u_gen) %in% tree_species, ]

#make sure tree and data match
obj = name.check(resolved_tree, filtered_u_gen)
new_tree = drop.tip(resolved_tree, obj$tree_not_data)

#log transformations
filtered_u_gen$log_u = log10(filtered_u_gen$u_per_year)
filtered_u_gen$log_gen_time = log10(filtered_u_gen$gen_time_yr)

#list to store PGLS results for each group
pgls_results = list()

#iterate through each group
for (group in unique(filtered_u_gen$group)) {
  
  #subset data for the group
  group_data = filtered_u_gen %>% filter(group == !!group)
  
  #prune tree to match group data
  group_tree = drop.tip(new_tree, setdiff(new_tree$tip.label, rownames(group_data)))
  
  #make sure species names are in a column
  group_data$species = rownames(group_data)
  
  #convert data to comparative object for caper
  comp_data = comparative.data(phy = group_tree, data = group_data, names.col = "species", vcv = T, warn.dropped = T)
  
  #fit PGLS model
  pgls_model = pgls(log_u ~ log_gen_time, data = comp_data, lambda = "ML")
  model_summary = summary(pgls_model)
  print(model_summary)

  #extract fitted values and assign species names
  fitted_values = setNames(fitted(pgls_model), rownames(comp_data$data))
  fitted_values = fitted_values[rownames(group_data)]
  
  #store results
  pgls_results[[group]] = list(
    coef = round(model_summary$coefficients[2, "Estimate"], 2),
    intercept = round(model_summary$coefficients[1, "Estimate"], 2),
    p_value = round(model_summary$coefficients[2, "Pr(>|t|)"], 4),
    se = round(model_summary$coefficients[2, "Std. Error"], 3),
    fitted_values = fitted_values, 
    adj_r2 = round(model_summary$adj.r.squared, 2),
    group_data = group_data
  )
}

#table of PGLS results
table_data = bind_rows(lapply(names(pgls_results), function(group) {
  res = pgls_results[[group]]
  data.frame(
    group = group,
    PGLS_coef = res$coef,
    p_value = res$p_value,
    R2 = res$adj_r2,
    n_species = nrow(res$group_data),
    SE_slope = res$se
  )
}))

table_grob = tableGrob(table_data)

#combine results for plotting
plot_data = bind_rows(lapply(names(pgls_results), function(group) {
  res = pgls_results[[group]]
  data.frame(
    log_gen_time_yr = res$group_data$log_gen_time,
    observed_values = res$group_data$log_u,
    fitted_values = res$fitted_values,
    group = group
  )
}))
    
#plot figures only
pdf("PGLS_yearly_mutation_gen_groups_plots.pdf", width = 10, height = 10)
fig1 = ggplot(plot_data, aes(x = log_gen_time_yr, y = observed_values)) +
  geom_point(colour = "red") +
  geom_line(aes(y = fitted_values), colour = "blue") +
  labs(x = "Log10 Generation Time (years)", y = "Log10 Mutation Rate (per year)") +
  facet_wrap(~group, scales = "free", ncol = 3) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 16, margin = margin(t = 10)), axis.title.y = element_text(size = 16, margin = margin(r = 10)))
print(fig1)
dev.off()  
  
combined_plot = plot_grid(fig1, table_grob, ncol = 2, rel_widths = c(2, 1))

#plot figures & table
pdf("PGLS_yearly_mutation_gen_groups_comb.pdf", width = 16, height = 10)
print(combined_plot)
dev.off()



#test for statistically significant difference between slopes & intercepts of two largest groups (mammals and reptiles)

#extract slopes, intercepts and SEs for mammals and reptiles
mammals_slope = pgls_results[["mammals"]][["coef"]]
birds_slope = pgls_results[["birds"]][["coef"]]
mammals_se = pgls_results[["mammals"]][["se"]]
birds_se = pgls_results[["birds"]][["se"]]
mammals_intercept = pgls_results[["mammals"]][["intercept"]]
birds_intercept = pgls_results[["birds"]][["intercept"]]

#calculate the difference in slopes and its standard error
slope_difference = mammals_slope - birds_slope
se_difference = sqrt(mammals_se^2 + birds_se^2)
intercept_difference = mammals_intercept - birds_intercept

#perform t-test on the difference between slopes
t_value = slope_difference / se_difference
df = min(nrow(pgls_results[["mammals"]][["group_data"]]) - 2, 
         nrow(pgls_results[["birds"]][["group_data"]]) - 2)  

p_value = 2 * pt(-abs(t_value), df)

cat("Mammals slope:", mammals_slope, "SE:", mammals_se, "\n")
cat("Birds slope:", birds_slope, "SE:", birds_se, "\n")
cat("Slope difference:", slope_difference, "\n")
cat("t-value:", t_value, "\n")
cat("p-value:", p_value, "\n")

#perform t-test on the difference between intercepts
t_value2 = intercept_difference / se_difference
df2 = min(nrow(pgls_results[["mammals"]][["group_data"]]) - 2, 
         nrow(pgls_results[["birds"]][["group_data"]]) - 2)

p_value2 = 2 * pt(-abs(t_value2), df2)

cat("Mammals intercept:", mammals_intercept, "SE:", mammals_se, "\n")
cat("Birds intercept:", birds_intercept, "SE:", birds_se, "\n")
cat("Intercept difference:", intercept_difference, "\n")
cat("t-value:", t_value2, "\n")
cat("p-value:", p_value2, "\n")


