#two objectives, unidimensional example MOP2
MOP2 <- function(x)
{
  xmod <- x*4 - 2
  if (is.null(nrow(x)))
  { 
    n <- length(xmod)
    y1 <- 1 - exp(-sum((xmod - 1/sqrt(n))^2) )
    y2 <- 1 - exp(-sum((xmod + 1/sqrt(n))^2) )
    Y <- matrix(c(y1,y2),1,2)
  } else
  {
    n <- ncol(xmod)
    y1 <- 1 - exp(-rowSums((xmod - 1/sqrt(n))^2) )
    y2 <- 1 - exp(-rowSums((xmod + 1/sqrt(n))^2) )
    Y <- cbind(y1,y2)
  }
  
  return(Y)
}

require(GPareto)
start_time <- Sys.time()
set.seed(25)
d <- 2
fname <- MOP2
plotParetoGrid(MOP2) 
end_time <- Sys.time()

a  <- (end_time - start_time)
a


start_time <- Sys.time()
require(mco)
r1 <- nsga2(MOP2, 2, 2,
            generations=150, popsize=200,
            cprob=0.7, cdist=20,
            mprob=0.2, mdist=20,
            lower.bounds=rep(0, 2),
            upper.bounds=rep(1, 2))
plot(r1,xlab ="f1", ylab = "f2")

end_time <- Sys.time()

b  <- (end_time - start_time)
b

#  two objectives, two dimensions example P1
P1 <-
  function(x){
    if(is.null(dim(x))){
      x <- matrix(x, nrow = 1) 
    }
    b1<-15*x[,1]-5
    b2<-15*x[,2]
    return(cbind((b2-5.1*(b1/(2*pi))^2+5/pi*b1-6)^2 +10*((1-1/(8*pi))*cos(b1)+1),
                 -sqrt((10.5-b1)*(b1+5.5)*(b2+0.5)) - 1/30*(b2 -5.1*(b1/(2*pi))^2-6)^2 - 1/3*((1-1/(8*pi))*cos(b1)+1)
    ) 
    )
  }

set.seed(25)
d <- 2
fname <- P1
plotParetoGrid(P1)

library(DiceDesign)
n_var <- 2
f_name <- "P1"
n.grid <- 26
test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
n_appr <- 15
design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 42)$design)$design, 1)
response.grid <- t(apply(design.grid, 1, f_name))
Front_Pareto <- t(nondominated_points(t(response.grid)))
mf1 <- km(~., design = design.grid, response = response.grid[,1])
mf2 <- km(~., design = design.grid, response = response.grid[,2])
EHI_grid <- crit_EHI(x = as.matrix(test.grid), model = list(mf1, mf2),
                     critcontrol = list(refPoint = c(300,0)))
filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
               matrix(EHI_grid, n.grid), main = "Expected Hypervolume Improvement",
               xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors,
               plot.axes = {axis(1); axis(2);
                 points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
               }
)
