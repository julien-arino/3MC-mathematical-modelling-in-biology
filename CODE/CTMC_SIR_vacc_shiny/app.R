#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# CTMC SIR with vaccination

library(shiny)
library(deSolve)
library(parallel)
library(adaptivetau)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Continuous-time Markov chain SIR with vaccination"),
    
    # Sidebar with a slider input for the number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("inv_gamma",
                        "Average infectious period (days):",
                        min = 0,
                        max = 30,
                        value = 7),
            sliderInput("R0",
                        "Basic reproduction number R0:",
                        min = 0,
                        max = 20,
                        value = 2,
                        step = 0.1),
            sliderInput("p",
                        "Percentage vaccinated:",
                        min = 0,
                        max = 100,
                        value = 0,
                        step = 1),
            sliderInput("tf",
                        "Final time:",
                        min = 1,
                        max = 200,
                        value = 100,
                        step = 1),
            sliderInput("nb_sims",
                        "Number of simulations:",
                        min = 1,
                        max = 200,
                        value = 25,
                        step = 1)
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("a_distPlot",width="750px")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    SIS <- function (t, x, parms) {
        with(as.list(x), {
            dS <- params$gamma*I-params$beta*S*I
            dI <- params$beta*S*I-params$gamma*I
            list(c(dS, dI))
        })
    }
    # Encode potential transitions:
    transitions_nodemog = cbind(c(-1,1,0), # New infection (-1S,+1I,0R)
                                c(0,-1,1)) # Recovery (0S,-1I,+1R)

    # Function to calculate transition rates, given variables and parameters
    lvrates_nodemog <- function(x, params, t) {
        return(c(params$beta*x["S"]*x["I"], # Rate of new infection
                 params$gamma*x["I"])) # Rate of recovery
    }
    
    one_sim_ssa_tau_leap <- function(sim) {
        # Perform the stochastic simulation!
        r = ssa.adaptivetau(IC, 
                            transitions_nodemog, 
                            lvrates_nodemog, 
                            params, 
                            tf = params$tf,
                            tl.params = list(epsilon=0.01))
        ii_interp <- approx(r[,"time"], r[,"I"], 
                            params$times, ties = "ordered", rule = 2)
        results <- data.frame(ii_interp$y)
        colnames(results) <- c("I")
        return(results)
    }

    # Expression that generates a histogram. The expression is
    # wrapped in a call to renderPlot to indicate that:
    #
    #  1) It is "reactive" and therefore should re-execute automatically
    #     when inputs change
    #  2) Its output type is a plot
    
    output$a_distPlot <- renderPlot({
        params <- list()
        
        params$nb_sims = input$nb_sims

        params$S0 = round(999*(1-input$p/100))
        params$I0 = 1
        IC = c(S = params$S0,
               I = params$I0,
               R = 0)
        
        params$t0 = 0
        params$tf = input$tf
        params$times <- seq(from = params$t0, 
                            to = params$tf, 
                            by = 0.1)
        params$length_times <- length(params$times)
        
        RESULTS <- mat.or.vec(params$nb_sims,
                              params$length_times)

        if (input$inv_gamma>0)
            params$gamma = 1/input$inv_gamma
        else
            params$gamma = 0
        #params$beta <- input$R0*params$gamma/(params$S0+params$I0)
        params$beta <- input$R0*params$gamma/params$S0
        
        inputs = 1:params[["nb_sims"]]
        
        numCores <- detectCores()
        cl <- makeCluster(numCores)  
        clusterExport(cl,c("params","IC","ssa.adaptivetau"),
                      envir=environment())
        RESULTS_TMP = parLapply(cl, inputs, one_sim_ssa_tau_leap)  
        stopCluster(cl)  
        
        RESULTS <- matrix(unlist(RESULTS_TMP), ncol = params[["length_times"]], byrow = TRUE)
        
        y_max = max(max(RESULTS))
        nb_sims_non_extinct = 0
        results_raw <- params[["times"]]*0
        results_conditioned <- params[["times"]]*0
        # Now plot the results
        for (sim in 1:params[["nb_sims"]]) {
            if (sim == 1) {
                plot(params[["times"]],RESULTS[sim,],
                     xlab="Time (days)",
                     ylab="Number of infectious individuals",
                     type="l",lwd=0.5,ylim = c(0,y_max))
            } else {
                lines(params[["times"]],RESULTS[sim,],lwd=0.5,type="l")
            }
            # Prepare sum of results for average
            results_raw = results_raw+RESULTS[sim,]
            # Prepare sum of results conditioned on non-extinction for average
            if (RESULTS[sim,params[["length_times"]]]>0) {
                results_conditioned <- results_conditioned+RESULTS[sim,]
                nb_sims_non_extinct <- nb_sims_non_extinct+1
            }
        }
        
        # Plot non-conditioned average
        results_raw = results_raw / params[["nb_sims"]]
        lines(params[["times"]],results_raw,col="blue",type="l",lwd=4)
        # Plot average conditioned on non-extinction
        results_conditioned = results_conditioned / nb_sims_non_extinct
        lines(params[["times"]],results_conditioned,col="red",type="l",lwd=4)
        
        legend("topleft", legend = c("Trajectories","Average (all realizations)","Average (non-extinction)"), 
               col=c("black","blue","red"), 
               lwd = c(1,2,2), lty = c(1,1,1))
        
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
