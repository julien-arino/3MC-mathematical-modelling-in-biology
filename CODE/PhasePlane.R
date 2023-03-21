# Stephanie Portet (University of Manitoba)
# 2023

library(deSolve)
library(ggplot2)
library(ggquiver)

xi=yi=0
xs=ys=7

#create a grid of points
grid.by = 0.2
seqx = seq(xi,xs,grid.by)
seqy = seq(yi,ys,grid.by)
seq0 = 0*seqx
grid_dir = expand.grid(x=seqx,y=seqy)


###xy-plane
## dy/dx=f(x,y)/g(x,y)
f<-function(x,y,p){
  res = p$a*x-p$c*x*y
  return(res)
}

g<-function(x,y,p){
  res = -p$b*y+p$d*x*y
  return(res)
}


# Right hand side of the ODE
rhs_ODE = function(t, z, p) {
  with(as.list(z), {
     dx = p$a*x-p$c*x*y
     dy = -p$b*y+p$d*x*y
     # dx = y*(x-1)
     # dy = x*(2-y)
    list(c(dx, dy))
  })
}

##
params = list()
params$x0 = 1
params$y0 = 1
params$a = 2
params$b = 1
params$c = 0.5
params$d = 1
params$tf = 100

## assemble the data in a data frame
Direction = data.frame(x=grid_dir$x,
                       y=grid_dir$y,
                       u_=f(grid_dir$x,grid_dir$y,params),
                       v_=g(grid_dir$x,grid_dir$y,params),
                       u=f(grid_dir$x,grid_dir$y,params)/sqrt(f(grid_dir$x,grid_dir$y,params)^2+g(grid_dir$x,grid_dir$y,params)^2),
                       v=g(grid_dir$x,grid_dir$y,params)/sqrt(f(grid_dir$x,grid_dir$y,params)^2+g(grid_dir$x,grid_dir$y,params)^2))

## plot the direction field
directionfield_=ggplot()+
  geom_quiver(data=Direction,aes(x=x,y=y,u=u,v=v))+
  theme_classic()
print(directionfield_)
####################################333333


IC <- c(x = params$x0, y = params$y0)
# Set time span. 
tspan = seq(from = 0, to = params$tf, by = 0.1)



sol_ODE = list()
sol_ODE = ode(y = IC,
          func = rhs_ODE,
          times = tspan,
          parms = params)

## put the solution in dataframe
Solution = data.frame(x=sol_ODE[,2],
                      y=sol_ODE[,3])

## plot direction field and trajectories and nullclines
directionfield_all=ggplot()+
  geom_quiver(data=Direction,aes(x=x,y=y,u=u,v=v))+
  geom_path(data = Solution,aes(x=x,y=y),color="red")+
  annotate("line",x=seqx,y=seq0,lty=1,col="orange")+
  annotate("line",x=params$b/params$d,y=seqy,lty=1,col="orange")+
  annotate("line",x=0,y=seqy,lty=1,col="green")+
  annotate("line",x=seqx,y=params$a/params$c,lty=1,col="green")+
  theme_classic()
print(directionfield_all)


plot(sol_ODE[,2],sol_ODE[,3],type = "l",xlab = "x",ylab = "y")
plot(sol_ODE[,1],sol_ODE[,2],type = "l",col="red",ylim=c(xi,xs),xlim = c(0,params$tf),xlab = "time",ylab = "Dependent var.")
lines(sol_ODE[,1],sol_ODE[,3],type = "l",col="blue")
legend("topleft",legend= c("x(t)","y(t)"),col = c("red","blue"),lty = c(1,1), horiz=TRUE, cex = 0.4)

