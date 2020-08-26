require(GPareto)

pow <- function(x,n){
  return(x^n)
}

M3 <-
  function(x){
    if(is.null(dim(x))){
      x <- matrix(x, nrow = 1) 
    }
    b1<-x[,1]
    b2<-x[,2]
    
    
    c1 <- 1-pnorm(0.05, mean =0.394, sd= 0.241, lower.tail = TRUE, log.p = FALSE)
    c2 <- 1-pnorm(0.05, mean =0.456, sd= 0.403, lower.tail = TRUE, log.p = FALSE)
    c3 <- 1-pnorm(0.05, mean =0.875, sd= 0.488, lower.tail = TRUE, log.p = FALSE)
    return(cbind(-b1*pow(1.394,0.9)/0.9 - b2*pow(1.456,0.75)/0.75 - (1-b1-b2)*pow(1.875,0.65)/0.65,
                 b1*(1.62+c1)+b2*(0.92+c2)+(1-b1-b2)*(1.67+c3)
    ) 
    )
  }

set.seed(25)
n_var <- 2
fname <- M3
lower <- rep(0, n_var)
upper <- rep(0.5, n_var)
res <- easyGParetoptim(fn=fname, lower=lower, upper=upper, budget=40,
                       control=list(method="EHI", inneroptim="pso", maxit=20))
par(mfrow=c(1,2))
plotGPareto(res)
title("Pareto Front")
plot(res$history$X, main="Pareto set", col = "red", pch = 15)
points(res$par, col="blue", pch = 17)

xy <- res$par

a <- xy[,-1]
b <- xy[,1]
c <- (1-a-b)

data <- as.data.frame(xy)
data1 <- as.data.frame(c)
write.table(data,file = "data.csv",sep = ",",quote = FALSE, append = FALSE, na = "NA")

dat <- read.csv('D:\\data.csv',header = TRUE)




plot(dat[,1], ylab="weight",main="weight of different assets in pareto set", type="l",lty=1)
lines(dat[,2], ann = F,col="red",lty=1)
lines(dat[,3], ann = F,col="blue",lty=1)


#legend("topright",cex=0.4,c("Gold","NASDAQ","Apple"),col=c("black","red","blue"),lty=1)
legend("topright", cex=0.2,inset=.01, c("Gold","NASDAQ","Apple"),lty=1, col=c("black","red","blue"))
