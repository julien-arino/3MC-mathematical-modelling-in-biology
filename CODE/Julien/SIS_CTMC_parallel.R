# Example CTMC simulation of a simple SIS model, parallel version
library(parallel)
library(GillespieSSA2)
library(deSolve)

# Source a file with a few helpful functions for plotting (nice axes labels, crop figure)
source("useful_functions.R")

# Need a function that runs one simulation and returns a result. While we're at it,
# we also return an interpolated solution
run_one_sim = function(params) {
  IC <- c(S = (params$Pop-params$I_0), I = params$I_0)
  params_local <- c(gamma = params$gamma, beta = params$beta)
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
    final_time = params$t_f,
    log_firings = TRUE    # This way we keep track of events as well
  )
  # Interpolate result (just I will do)
  wanted_t = seq(from = 0, to = params$t_f, by = 0.01)
  sol$interp_I = approx(x = sol$time, y = sol$state[,"I"], xout = wanted_t)
  names(sol$interp_I) = c("time", "I")
  # Return result
  return(sol)
}

# The ODE, in order to compare with solutions to the CTMC
rhs_SIS_ODE = function(t, x, p) {
  with(as.list(x), {
    dS = p$gamma*I-p$beta*S*I
    dI = p$beta*S*I-p$gamma*I
    list(c(dS, dI))
  })
}

# To run in parallel, it useful to put parameters in a list
params = list()
params$Pop = 100
params$gamma = 1/5
params$R_0 = 1.5
params$t_f = 100
params$I_0 = 1
# R0 would be (beta/gamma)*S0, so beta=R0*gamma/S0
params$beta = params$R_0*params$gamma/(params$Pop-params$I_0)
# Number of simulations. We may want to save all simulation parameters later,
# so we add it here
params$number_sims = 50


# Detect number of cores (often good to use all but 1, i.e. detectCores()-1)
# Here, we also test if we are on a 3990X or similar, as R doesn't do 128 threads
nb_cores <- detectCores()
if (nb_cores > 124) {
  nb_cores = 124
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
                "run_one_sim"),
              envir = .GlobalEnv)
# Run computation
tictoc::tic()
SIMS = parLapply(cl = cl, 
                 X = 1:params$number_sims, 
                 fun =  function(x) run_one_sim(params))
tictoc::toc()
stopCluster(cl)

# The following is if running iteratively rather than in parallel
if (FALSE) {
  SIMS = lapply(X = 1:params$number_sims, 
                FUN =  function(x) run_one_sim(params))
}

# To compute means of trajectories, we have to work a little. There are two
# types of means: regular mean, using all trajectories (even those having reached
# zero) and mean conditioned on non extinction (mean of trajectories not absorbed
# at zero during the time interval under consideration). We prepare a list that 
# will contain both means.
mean_I = list(time = SIMS[[1]]$interp_I$time,
              I_all = rep(0, length(SIMS[[1]]$interp_I$time)),
              I_no_extinction = rep(0, length(SIMS[[1]]$interp_I$time)))
nb_sims_with_NA = 0
for (i in 1:params$number_sims) {
  if (any(is.na(SIMS[[i]]$interp_I$I))) {
    # If we do not condition on extinction, we need to set NAs to 0
    tmp = SIMS[[i]]$interp_I$I
    tmp[which(is.na(tmp))] = 0
    mean_I$I_all = mean_I$I_all+tmp
    nb_sims_with_NA = nb_sims_with_NA+1
  } else {
    mean_I$I_all = mean_I$I_all + SIMS[[i]]$interp_I$I
    mean_I$I_no_extinction = mean_I$I_no_extinction + SIMS[[i]]$interp_I$I
  }
}
mean_I$I_all = mean_I$I_all/params$number_sims
mean_I$I_no_extinction = mean_I$I_no_extinction/(params$number_sims-nb_sims_with_NA)

# Now simulate the ODE
IC <- c(S = (params$Pop-params$I_0), I = params$I_0)
sol_ODE = ode(y = IC,
              func = rhs_SIS_ODE,
              times = seq(from = 0, to = params$t_f, by = 0.1),
              parms = params)

# Determine maximum value of I for plot
I_max = max(unlist(lapply(SIMS, function(x) max(x$interp_I$I, na.rm = TRUE))))
# Prepare y-axis for human readable form
y_axis = make_y_axis(c(0, I_max))

# Are we plotting for a dark background
plot_blackBG = TRUE
if (plot_blackBG) {
  colour = "white"
} else {
  colour = "black"
}

# Plot
png(file = "../FIGS/many_CTMC_sims_with_means.png",
    width = 1200, height = 800, res = 200)
if (plot_blackBG) {
  par(bg = 'black', fg = 'white') # set background to black, foreground white
}
plot(mean_I$time, SIMS[[1]]$interp_I$I,
     type = "l", lwd = 0.2, 
     ylim = c(0, max_I), xaxs = "i",
     col.axis = colour, cex.axis = 1.25,
     col.lab = colour, cex.lab = 1.1,
     yaxt = "n",
     xlab = "Time (days)", ylab = "Prevalence")
for (i in 2:params$number_sims) {
  lines(mean_I$time, SIMS[[i]]$interp_I$I,
        type = "l", lwd = 0.2)
}
lines(mean_I$time, mean_I$I_all,
      type = "l",
      lwd = 5, col = "darkorange4")
lines(mean_I$time, mean_I$I_no_extinction,
      type = "l",
      lwd = 5, col = "red")
lines(sol_ODE[,"time"], sol_ODE[,"I"],
      type = "l",
      lwd = 5, col = "dodgerblue4")
legend("topleft",
       legend = c("Solutions", "Mean", 
                  "Mean (not extinct)", "ODE"),
       cex = 0.6,
       col = c(ifelse(plot_blackBG, "white", "black"), 
               "darkorange4", "red", "dodgerblue"),
       lty = c(1,1,1,1), lwd = c(0.5, 2.5, 2.5, 2.5))
axis(2, at = y_axis$ticks, labels = y_axis$labels, 
     las = 1,
     col.axis = colour,
     cex.axis = 0.75)
dev.off()
crop_figure(file = "../FIGS/many_CTMC_sims_with_means.png")
