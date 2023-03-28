# Example CTMC simulation of a simple SIS model, parallel version, with multiple R_0 values
# Here, we just want to know how many simulations see extinction, so we do minimal post-processing

# BEWARE !!!
#
# With 100 sims and the number of values used, this can either be very lengthy on a machine with low thread
# count or very expensive in RAM on a machine with high thread count..

library(GillespieSSA2)
library(parallel)
library(dplyr)
library(latex2exp)

# Source a file with a few helpful functions for plotting (nice axes labels, crop figure)
source("useful_functions.R")

# Make data frame from list of results
# Assumes all list entries are similar
make_df_from_list = function(in_list) {
  names_list = names(in_list[[1]])
  tmp = list()
  for (n in names_list) {
    tmp[[n]] = unlist(lapply(in_list, function(x) x[[n]]))
  }
  OUT = do.call(cbind.data.frame, tmp)
  return(OUT)
}

# Need a function that runs one simulation and returns a result. While we're at it,
# we also return an interpolated solution
run_one_sim = function(params_vary, params) {
  R_0 = params_vary$R_0
  I_0 = params_vary$I_0
  # Initial conditions
  IC <- c(S = (params$Pop-I_0), I = I_0)
  # R0=(beta/gamma)*S0, so beta=R0*gamma/S0
  beta = R_0*params$gamma/(params$Pop-params$I_0)
  params_local <- c(gamma = params$gamma, beta = beta)
  reactions <- list(
    # propensity function effects name for reaction
    reaction("beta*S*I", c(S=-1,I=+1), "new_infection"),
    reaction("gamma*I", c(S=+1,I=-1), "recovery")
  )
  set.seed(NULL)
  sol <- ssa(
    initial_state = IC,
    reactions = reactions,
    params = params_local,
    method = ssa_exact(),
    final_time = params$t_f
  )
  # If the final time is less than t_f, we've hit extinction
  if (sol$state[dim(sol$state)[1],"I"] == 0) {
    sol$extinct = TRUE
  } else {
    sol$extinct = FALSE
  }
  # To lighten the load, select only a few items to return
  OUT = list()
  OUT$R_0 = R_0
  OUT$I_0 = I_0
  OUT$final_time = sol$time[length(sol$time)]
  OUT$final_I = sol$state[length(sol$time),"I"]
  OUT$extinct = sol$extinct
  # Return result
  return(OUT)
}

# To run in parallel, it useful to put parameters in a list
params = list()
params$Pop = 100
params$gamma = 1/5
params$R_0 = 1.5
params$t_f = 100
params$I_0 = 1
# R0=(beta/gamma)*S0, so beta=R0*gamma/S0
params$beta = params$R_0*params$gamma/(params$Pop-params$I_0)
# Number of simulations. We may want to save all simulation parameters later,
# so we add it here
params$number_sims = 500

# Detect number of cores (often good to use all but 1, i.e. detectCores()-1)
# Here, we also test if we are on a 3990X or similar, as R doesn't do 128 threads
nb_cores <- detectCores()
if (nb_cores > 124) {
  nb_cores = 124
}
# Values of I_0 we consider (we do those sequentially)
values_I_0 = c(1,2,3,5,10)
# Run main computation loop (iterative part)
for (I_0 in values_I_0) {
  writeLines(paste0("I_0=",I_0))
  # To process efficiently in parallel, we make a list with the different parameter values
  # we want to change, which will fed by parLapply to run_one_sim
  params_vary = list()
  # Vary R_0 and I_0
  i = 1
  for (R_0 in c(seq(0.5, 3, by = 0.05))) {
    for (j in 1:params$number_sims) {
      params_vary[[i]] = list()
      params_vary[[i]]$R_0 = R_0
      params_vary[[i]]$I_0 = I_0
      i = i+1
    }
  }
  # Initiate cluster
  cl <- makeCluster(nb_cores)
  # Export needed library to cluster
  clusterEvalQ(cl,{
    library(GillespieSSA2)
  })
  # Export needed variable and function to cluster
  clusterExport(cl,
                c("params",
                  "run_one_sim",
                  "params_vary"),
                envir = .GlobalEnv)
  # Run main computation (parallel part)
  SIMS = parLapply(cl = cl, 
                 X = params_vary, 
                 fun =  function(x) run_one_sim(x, params))
  stopCluster(cl)
  saveRDS(SIMS, file = sprintf("SIMS_I0_%02d.Rds", I_0))
}


#SIMS = readRDS(file = "SIMS_I0_1.Rds")

# Use dplyr syntax: count the number of extinctions (TRUE and FALSE)
# after grouping by R_0 value
results = make_df_from_list(SIMS) %>%
  count(I_0, R_0, extinct)
# Take the result and prepare to negate those where n=params$number_sims
tmp = results %>%
  filter(n == params$number_sims)
tmp_insert = tmp
tmp_insert$extinct = ifelse(tmp$extinct == TRUE, FALSE, TRUE)
tmp_insert$n = rep(0, dim(tmp_insert)[1])
# Now inject this result and keep extinct==TRUE
results = rbind(results, tmp_insert) %>%
  filter(extinct == TRUE) %>%
  arrange(I_0, R_0)
results$pct_extinct = results$n/params$number_sims*100

# Values I_0
values_I_0 = unique(results$I_0)

# Are we plotting for a dark background
plot_blackBG = TRUE
if (plot_blackBG) {
  colour = "white"
} else {
  colour = "black"
}

# Prepare y-axis for human readable form
y_axis = make_y_axis(c(0, 100))

# Plot
png(file = "../FIGS/extinctions_fct_R0_I0.png",
    width = 1200, height = 800, res = 200)
if (plot_blackBG) {
  par(bg = 'black', fg = 'white') # set background to black, foreground white
}
tmp = results %>%
  filter(I_0 == values_I_0[1])
plot(tmp$R_0, tmp$pct_extinct,
     type = "b", lwd = 2,
     ylim = c(0, 100), yaxt = "n",
     col.axis = colour, cex.axis = 1.25,
     col.lab = colour, cex.lab = 1.1,
     xlab = TeX("$R_0$"), ylab = "Percentage extinctions")
for (i in 2:length(values_I_0)) {
  tmp = results %>%
    filter(I_0 == values_I_0[i])
  lines(tmp$R_0, tmp$pct_extinct,
        type = "l", lwd = 2, lty = i)
}
abline(v = 1, lty = "dotted")
# legend("topright",
#        legend = TeX(sprintf("$I_0$=%d", values_I_0)),
#        lty = 1:length(values_I_0), lwd = rep(2, length(values_I_0)))
axis(2, at = y_axis$ticks, labels = y_axis$labels, 
     las = 1,
     col.axis = colour,
     cex.axis = 0.75)
dev.off()
crop_figure(file = "../FIGS/extinctions_fct_R0_I0.png")

