library(deSolve)
library(latex2exp)

# Source a file with a few helpful functions for plotting (nice axes labels, crop figure)
source(sprintf("%s/../functions_useful.R", here::here()))

# Right hand side of the ODE
rhs_SIS_ODE = function(t, x, p) {
  with(as.list(x), {
    change = 
    dS = p$gamma*I-p$beta*S*I
    dI = p$beta*S*I-p$gamma*I
    list(c(dS, dI))
  })
}

# Compute EP for parameters
value_EP = function(params) {
  OUT = list()
  if (params$R_0<1) {
    OUT$S_EP = params$Pop
    OUT$E_EP = 0
    OUT$col = "dodgerblue4"
  } else {
    OUT$S_EP = 1/params$R_0*params$Pop
    OUT$E_EP = (1-1/params$R_0)*params$Pop
    OUT$col = "darkorange4"
  }
  return(OUT)
}

# It useful to put parameters in a list
params = list()
params$Pop = 100000
params$gamma = 1/5
params$R_0 = 1.5
params$t_f = 150
params$I_0 = 2
# R0 is (beta/gamma)*S0, so beta=R0*gamma/S0
params$beta = params$R_0*params$gamma/(params$Pop-params$I_0)

# Initial conditions. It is useful to have named variables here (S=.., I=..),
# so that the ODE algorithm knows what variables are passed to it (otherwise, 
# the RHS would need S=x[1], etc.)
IC <- c(S = (params$Pop-params$I_0), I = params$I_0)

# Set time span. Can be either end points or specified time points. We do the latter
# (when fitting, which we are not doing here, this is useful)
tspan = seq(from = 0, to = params$t_f, by = 0.1)

# Now simulate the ODE. We will loop on several values of R_0, to make the first
# figure in the slides
R_0 = c(0.8, 1.5, 2.5)
# Save results in a list
sol_ODE = list()
EP = list()
# Now loop on R_0
for (r in R_0) {
  params$R_0 = r
  # Name for list entry
  entry_name = sprintf("$R_0$=%1.1f",r)
  # Need to recompute beta for current R_0
  params$beta = params$R_0*params$gamma/(params$Pop-params$I_0)
  # Call numerical integrator
  sol_ODE[[entry_name]] = ode(y = IC,
                              func = rhs_SIS_ODE,
                              times = tspan,
                              parms = params)
  EP[[entry_name]] = value_EP(params)
  EP[[entry_name]]$lty = which(r == R_0)
}

# Get maximum value of I across all simulations for plot. Note the use of lapply.
max_I = max(unlist(lapply(sol_ODE, function(x) max(x[,"I"]))))

# Plot
png(file = sprintf("%s/FIGURES/ODE_SIS.png", here::here()),
    width = 1200, height = 800, res = 200)
y_axis = plot_hr_yaxis(sol_ODE[[1]][,"time"], sol_ODE[[1]][,"I"],
                       y_range = c(0, max_I),
                       type = "l", lwd = 5, col = EP[[1]]$col, lty = EP[[1]]$lty,
                       xlab = "Time (days)", ylab = "Number infectious")
points(x = params$t_f, y = EP[[1]]$E_EP*y_axis$factor, 
       col = EP[[1]]$col, pch = 19, cex = 2)
for (i in 2:length(sol_ODE)) {
  lines(sol_ODE[[i]][,"time"], sol_ODE[[i]][,"I"]*y_axis$factor,
        type = "l", lwd = 5, col = EP[[i]]$col, lty = EP[[i]]$lty)
  points(x = params$t_f, y = EP[[i]]$E_EP*y_axis$factor, 
         col = EP[[i]]$col, pch = 19, cex = 2)
}
legend("topleft", legend = TeX(names(EP)), cex = 0.8,
       col = unlist(lapply(EP, function(x) x$col)),
       lty = unlist(lapply(EP, function(x) x$lty)),
       lwd = c(3,3,3))
dev.off()
crop_figure(file = sprintf("%s/FIGURES/ODE_SIS.png", here::here()))

