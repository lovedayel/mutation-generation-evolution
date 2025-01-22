set.seed(21)

data = read.csv("mutation_rates_final.csv", h = T)

#calculate observed variance for real data 
log_generation_times = log(data$gen_time_yr)
Vobs_LGT = var(log_generation_times, na.rm = T)
print(Vobs_LGT)

n = 147 
observed_slope = 0.31  
num_repetitions = 100

calculate_slope = function(x, y) {
  cov(x, y) / var(x)
}

simulate_GT_mutation = function(V_e, num_reps) {
  slopes <- numeric(num_reps)
  
  for (rep in 1:num_reps) {
    Vtrue_LGT = Vobs_LGT - V_e
    
    if (Vtrue_LGT <= 0) {
      return(NULL)
    }
    
    log_GT_true = rnorm(n, mean = 1, sd = sqrt(Vtrue_LGT))
    error = rnorm(n, mean = 0, sd = sqrt(V_e))
    log_GT_obs = log_GT_true + error
    
    constant = runif(1, -1, 1)
    log_u_gen = log_GT_true + constant
    
    slopes[rep] = calculate_slope(log_GT_obs, log_u_gen)
  }
  
  mean_slope = mean(slopes)
  return(list(slope = mean_slope, V_e = V_e, Vtrue_LGT = Vtrue_LGT))
}

#run simulation for different V_e values
V_e_values = seq(0, Vobs_LGT, length.out = 20)
results = lapply(V_e_values, function(v) simulate_GT_mutation(v, num_repetitions))
results = results[!sapply(results, is.null)]

slopes = sapply(results, function(x) x$slope)
V_e_ratio = sapply(results, function(x) x$V_e / x$Vtrue_LGT)

#find the V_e ratio that produces the observed slope
target_index = which.min(abs(slopes - observed_slope))
target_ratio = V_e_ratio[target_index]

pdf("gen_time_simulation2.pdf", width = 7, height = 5)
fig1 = plot(V_e_ratio, slopes, type = "l", 
     xlab = "V(e) / V(LGT)", ylab = "Mean Slope")
abline(h = observed_slope, col = "red", lty = 2)
abline(v = target_ratio, col = "blue", lty = 2)
print(fig1)
dev.off()

cat("To achieve the observed slope of", observed_slope, 
    "the ratio V(e) / Vtrue(LGT) needs to be approximately", round(target_ratio, 2), "\n")