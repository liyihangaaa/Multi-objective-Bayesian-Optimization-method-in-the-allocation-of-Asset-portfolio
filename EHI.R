## Analytical expression of the Expected Hypervolume Improvement with 2 objectives
## To avoid numerical instabilities, the new point is evaluated only if it is not too close to an existing observation.



EHI<- function(x, model, critcontrol=NULL, paretoFront = NULL){
  
  n.obj <- length(model)
  d <- model[[1]]@d
  if (!is.matrix(x)) x <- matrix(x, 1, d)  
  n.candidates <- nrow(x)
  
  if(is.null(paretoFront) || is.null(critcontrol$refPoint)){
    observations <- Reduce(cbind, lapply(model, slot, "y"))
    if(is.null(paretoFront))
      paretoFront <- t(nondominated_points(t(observations)))
  }
  
  if (is.unsorted(paretoFront[,1])){
    paretoFront <- paretoFront[order(paretoFront[,1]),]
  }
  
  refPoint <- critcontrol$refPoint
  if (is.null(refPoint)){
    if(is.null(critcontrol$extendper)) critcontrol$extendper <- 0.2
    # refPoint    <- matrix(apply(paretoFront, 2, max) + 1, 1, n.obj)
    PF_range <- apply(paretoFront, 2, range)
    refPoint <- matrix(PF_range[2,] + pmax(1, (PF_range[2,] - PF_range[1,]) * critcontrol$extendper), 1, n.obj)
    cat("No refPoint provided, ", signif(refPoint, 3), "used \n")
  }else {
    
    # mu    <- rep(NaN, n.obj)
    # sigma <- rep(NaN, n.obj)
    # for (i in 1:n.obj){    
    #   pred     <- predict(object=model[[i]], newdata=x.new, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
    #   mu[i]    <- pred$mean
    #   sigma[i] <- pred$sd
    # }
    pred <- predict_kms(model, newdata=x, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
    mu <- t(pred$mean)
    sigma <- t(pred$sd)
    
    ## A new x too close to the known observations could result in numerical problems
    # check <- apply(x, 1, checkPredict, model = model, type = type, distance = critcontrol$distance, threshold = critcontrol$threshold)
    check <- checkPredict(x, model, threshold = critcontrol$threshold, distance = critcontrol$distance, type = type)
    resu <- rep(0, n.candidates)
    resu[check] <- -1
    resu[!check] <- EHI_2d_wrap_Rcpp(paretoFront, refPoint, mu[!check,,drop=FALSE], sigma[!check,,drop=FALSE])
    return(resu)
  }
}
