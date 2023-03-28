# Running a birth-death process as a CTMC using classic Gillespie
# 
# This demos both the classic Gillespie algorithm and how things can go bad
# because of shrinking inter-event times.

b = 0.025   # Birth rate
d = 0.01   # Death rate
t_0 = 0    # Initial time
N_0 = 100  # Initial population

# Vectors to store time and state. Initialise with initial condition.
t = t_0
N = N_0

t_f = 1000  # Final time

# We'll track the current time and state (could also just check last entry in t
# and N, but will take more operations)
t_curr = t_0
N_curr = N_0

while (t_curr<=t_f) {
  xi_t = (b+d)*N_curr
  # The exponential number generator does not like a rate of 0 (when the 
  # population crashes), so we check if we need to quit
  if (N_curr == 0) {
    break
  }
  tau_t = rexp(1, rate = xi_t)
  t_curr = t_curr+tau_t
  v = c(b*N_curr, xi_t)/xi_t
  zeta_t = runif(n = 1)
  pos = findInterval(zeta_t, v)+1
  switch(pos,
         { 
           # Birth
           N_curr = N_curr+1
         },
         {
           # Death
           N_curr = N_curr-1
         })
  N = c(N, N_curr)
  t = c(t, t_curr)
}

# Make a title for all figures
main_title = paste0("b=",b,", d=",d)
fig_title = gsub("\\.", "_", paste0("b=",b,"__d=",d))

# Do we output the results?
OUTPUT_PLOT = TRUE

# Plot the result
if (OUTPUT_PLOT) {
  png(filename = paste0("../FIGS/CTMC_birth_death_sol_", fig_title, ".png"),
      width = 800, height = 600, res = 120)
}
plot(t, N, type = "l",
     main = main_title)
abline(h = N_0, lty = 2)
if (OUTPUT_PLOT) {
  dev.off()
  crop_figure(paste0("../FIGS/CTMC_birth_death_sol_", fig_title, ".png"))
}

# Plot the inter-event time. Diff removes one entry, beware..
plot(1:(length(t)-1), diff(t),
     xlab = "Jump #", ylab = "Inter-event time",
     main = main_title,
     type = "l")

# Plot the inter-event time. Diff removes one entry, beware..
if (OUTPUT_PLOT) {
  png(filename = paste0("../FIGS/CTMC_birth_death_ie_vs_t_", fig_title, ".png"),
      width = 800, height = 600, res = 120)
}
plot(t[2:length(t)], diff(t),
     xlab = "Time", ylab = "Inter-event time",
     main = main_title,
     type = "l")
if (OUTPUT_PLOT) {
  dev.off()
  crop_figure(paste0("../FIGS/CTMC_birth_death_ie_vs_t_", fig_title, ".png"))
}
