#do some simulations! 

rm(list = ls())
setwd("//tsclient/C/Users/Mijke/Documents/__STUDIES__/CausalIndicators/MediationModelSimulation")

setwd("C:/Users/Mijke/Documents/__STUDIES__/CausalIndicators/MediationModelSimulation")
source("Functions.R")
windowsFonts(A=windowsFont("Verdana"))

######################## SIM 2.1 ###################################################
#replicate sim3 from Study 1, except set constant residual variance of M at .25
# reliability (M) = 1/1.25 = .8
# now both models are misspecified. 

sim2.1 <- function(i){
  a <- .6 #increase this to increase the range of possible values.  
  b <- .4
  res.M <- .25
  avec <- runif(4, 0, .8) 
  data <- generateCov(a = a, b = b, avec = avec, covm = NA, res.M = res.M)
  while (min(eigen(data[-c(5,8),-c(5,8)])$values) < 1e-10) {    
    avec <- runif(4, 0, .8) 
    data <- generateCov(a = a, b = b, avec = avec, covm = NA, res.M = res.M)
  }
  meana <- mean(avec)
  vara <- var(avec)
  rangea <- range(avec)[2] - range(avec)[1]
  amat <- avec %*% t(avec)
  gamma <- a/sum(avec)
  covm <- (1 - 4*gamma^2 - res.M)/(12*gamma^2)   #NEW so total var(M) includes residual variance
  
  fitComp <- sem(model = CompMod, sample.cov = data, sample.nobs = 100000, warn = FALSE)  
  fitLV <- try(sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE, warn = FALSE), silent = TRUE)  
    
  conv.LV <- try(inspect(fitLV, "converged"), silent = TRUE)
  SEs.LV <- TRUE  
  if (identical(class(try(vcov(fitLV), silent = TRUE)), "try-error") == TRUE) {SEs.LV <- FALSE}
  if (conv.LV == FALSE) {fit.LV <- NA} else (fit.LV <- inspect(fitLV, "fit")["fmin"])  
    
  return(c(meana, rangea, vara, covm, gamma,
           inspect(fitLV, "std")$lambda[1:4,1],
           inspect(fitComp, "std")$beta[1,3], #a path composite model standardized
           inspect(fitComp, "std")$beta[2,1], #b path composite model standardized 
           inspect(fitLV, "std")$beta[1,3], #a path LV model standardized
           inspect(fitLV, "std")$beta[2,1], #b path LV model standardized
           inspect(fitComp, "fit")["fmin"], #fit of composite model
           conv.LV, SEs.LV, fit.LV))    
}

X <- 1:10000
ptm <- proc.time()
out <- sapply(X, sim2.1)
proc.time() - ptm

save(out, file = "sim2.1.Rdata")
load(file = "sim2.1.Rdata")
results <- as.data.frame(t(out))

#new to reflect added columns: 
colnames(results) <- c("mean.a", "range.a", "var.a", "covm", "gamma", 
                       "lambda1.std", "lambda2.std", "lambda3.std", "lambda4.std",
                       "CompEst.a.std", "CompEst.b.std", "LVEst.a.std", "LVEst.b.std", 
                       "Comp.fmin", "conv.LV", "SEconv.LV", "LV.fmin")


#fix gamma and covm values: 
results$gamma <- .6 /(4*results$mean.a)
results$covm <- (1/(12*results$gamma^2)) - 1/3
head(results)

save(results, file = "sim2.1.fixed.Rdata")

#discard rows that did not converge or that generated non-positive-definite sigma: 
results2 <- results[is.na(results$LV.fmin) == FALSE,]
results2 <- results2[results2$SEconv.LV == 1,]
dim(results2) #4666

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
ryPal <- colorRampPalette(c('red3','yellow'))

results2$col.var.a <- rbPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]
results2$col.covm <- ryPal(10)[as.numeric(cut(results2$covm,breaks = 10))]

#compute indirect effects:
results2$LVEst.ab <- results2$LVEst.a * results2$LVEst.b
results2$CompEst.ab <- results2$CompEst.a * results2$CompEst.b
results2$LVEst.ab.std <- results2$LVEst.a.std * results2$LVEst.b.std
results2$CompEst.ab.std <- results2$CompEst.a.std * results2$CompEst.b.std

#compute fit: 
results2$pop.RMSEA <- sqrt(results2$LV.fmin/9)
results2$pop.RMSEA.comp <- sqrt(results2$Comp.fmin/1)

#plot a and b estimates:
jpeg("sim2.1.Plot1.jpg", width = 12, height = 4, units = "in", res = 800) 
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, 1/3, 0, 1))
  plot(results2$covm, results2$CompEst.a, type = "p", xlab = "covariance between causal indicators",
       ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 5), pch = ".", cex = 2, xlim = c(0,1))
  points(results2$covm, results2$LVEst.a, col = results2$col.var.a, pch = ".", cex = 2)
  
  par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
  plot(results2$covm, results2$CompEst.b, type = "p", xlab = "covariance between causal indicators",
       ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 2, xlim = c(0,1))
  points(results2$covm, results2$LVEst.b, col = results2$col.var.a, pch = ".", cex = 2) 
  legend(x = -.1, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "darkred"), bty = "n")
  
  par(fig = c(2/3, 1, 0, 1), new = TRUE)
  plot(results2$covm, results2$CompEst.ab, type = "p", 
       xlab = "covariance between causal indicators",
       ylab = "estimate of a*b path", col = "darkgreen", ylim = c(.22, .32), pch = ".", cex = 2, xlim = c(0,1))
  points(results2$covm, results2$LVEst.ab, col = results2$col.var.a, pch = ".", cex = 2)
dev.off()  


#plot a and b standardized estimates:
jpeg("sim2.1.Plot2.jpg", width = 12, height = 4, units = "in", res = 800) 
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, 1/3, 0, 1))
  plot(results2$covm, results2$CompEst.a.std, type = "p", xlab = "correlation between causal indicators",
       ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 1), pch = ".", xlim = c(0,1), cex = 2)
  abline (h = .6)  
  points(results2$covm, results2$LVEst.a.std, col = results2$col.var.a, pch = ".", cex = 2)

  par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
  plot(results2$covm, results2$CompEst.b.std, type = "p", xlab = "correlation between causal indicators",
       ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = ".", xlim = c(0,1), cex = 2)
  abline (h = .4)  
  points(results2$covm, results2$LVEst.b.std, col = results2$col.var.a, pch = ".", cex = 2) 
  legend(x = 0, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "darkred"), bty = "n")

  par(fig = c(2/3, 1, 0, 1), new = TRUE)
  plot(results2$covm, results2$CompEst.ab.std, type = "p", 
       xlab = "correlation between causal indicators", xlim = c(0,1),
       ylab = "estimate of a*b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 2)
  abline (h = .24)  
  points(results2$covm, results2$LVEst.ab.std, col = results2$col.var.a, pch = ".", cex = 2)
dev.off()  


#plot fit
jpeg("sim2.1.Plot3.jpg", width = 9, height = 4, units = "in", res = 800) 
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  plot(results2$covm, results2$pop.RMSEA, type = "p", xlab = "correlation between causal indicators",
       ylab = "population RMSEA", col = results2$col.var.a, ylim = c(0, 1), 
       pch = ".", cex = 2)
  points(results2$covm, results2$pop.RMSEA.Comp, col = "darkgreen", pch = ".", cex = 2)
dev.off()  


# Study 2.2. vary zeta and hold covm = .4 constant. #####

sim2.2 <- function(i){
  a <- .6 
  b <- .4
  covm <- .4
  meana <- sqrt((3*a^2/4)*(covm + (1/3))) #phew
  a4 <- 1
  while (a4 > .8 || a4 < 0){
    avec <- runif(3, 0, .8) #in sim2 this had max = .5
    a4 <- meana*4 - sum(avec)    
  }
  avec <- c(avec, a4)   
  vara <- var(avec)
  rangea <- range(avec)[2] - range(avec)[1]
  rho.M <- runif(1, 0, 1)
  res.M <- (1/rho.M) - 1  
  data <- generateCov(a = a, b = b, avec = avec, covm = NA, res.M = res.M)
  
  if (min(eigen(data[-5,-5])$values) > 0) { #ensure that sigma is positive definite! 
    gamma <- a/sum(avec)
    covm <- (1/(12*gamma^2)) - 1/3    
    
    fitComp <- sem(model = CompMod, sample.cov = data, sample.nobs = 100000)  
    fitLV <- sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  
    
    conv.LV <- try(inspect(fitLV, "converged"))
    SEs.LV <- TRUE  
    if (class(try(vcov(fitLV))) == "try-error") {SEs.LV <- FALSE}
    if (conv.LV == FALSE) {fit.LV <- NA} else (fit.LV <- inspect(fitLV, "fit")["fmin"])  
    
    return(c(meana, rangea, vara, covm, gamma, res.M,
             inspect(fitLV, "coef")$lambda[1:4,1],
             inspect(fitComp, "coef")$beta[1,3], #a path composite model
             inspect(fitComp, "coef")$beta[2,1], #b path composite model  
             inspect(fitLV, "coef")$beta[1,3], #a path LV model
             inspect(fitLV, "coef")$beta[2,1], #b path LV model
             inspect(fitLV, "std")$lambda[1:4,1],
             inspect(fitComp, "std")$beta[1,3], #a path composite model standardized
             inspect(fitComp, "std")$beta[2,1], #b path composite model standardized 
             inspect(fitLV, "std")$beta[1,3], #a path LV model standardized
             inspect(fitLV, "std")$beta[2,1], #b path LV model standardized
             inspect(fitComp, "fit")["fmin"], #fit of composite model
             conv.LV, SEs.LV, fit.LV))    
  } else {return (rep(NA, 26))}
}

X <- 1:10000
ptm <- proc.time()
out2.2 <- sapply(X, sim2.2)
proc.time() - ptm

save(out2.2, file = "sim2.2.10000.Rdata")
load(file = "sim2.2.10000.Rdata")
results <- as.data.frame(t(out2.2))

colnames(results) <- c("mean.a", "range.a", "var.a", "covm", "gamma", "res.M",
                       "lambda1", "lambda2", "lambda3", "lambda4",
                       "CompEst.a", "CompEst.b", "LVEst.a", "LVEst.b", 
                       "std.lam1", "std.lam2", "std.lam3", "std.lam4",
                       "std.CompEst.a", "std.CompEst.b", "std.LVEst.a", "std.LVEst.b", 
                       "Comp.fmin", "conv.LV", "SEconv.LV", "LV.fmin")


#discard rows that did not converge or that generated non-positive-definite sigma: 
results2 <- results[is.na(results$LV.fmin) == FALSE,]
results2 <- results2[results2$SEconv.LV == 1,]
dim(results2) #9636 rows left

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
ryPal <- colorRampPalette(c('red3','yellow'))

results2$col.var.a <- rbPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]
results2$col.covm <- ryPal(10)[as.numeric(cut(results2$covm,breaks = 10))]

#compute indirect effects:
results2$LVEst.ab <- results2$LVEst.a * results2$LVEst.b
results2$CompEst.ab <- results2$CompEst.a * results2$CompEst.b
results2$std.LVEst.ab <- results2$std.LVEst.a * results2$std.LVEst.b
results2$std.CompEst.ab <- results2$std.CompEst.a * results2$std.CompEst.b

#compute fit: 
results2$pop.RMSEA <- sqrt(results2$LV.fmin/9)
results2$pop.RMSEA.comp <- sqrt(results2$Comp.fmin/1)

results2$rho.M <- 1/(1+results2$res.M)

#plot a and b standardized estimates:
jpeg("sim2.2.Plot1.jpg", width = 12, height = 8, units = "in", res = 500) 
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, 1/3, .5, 1))
  plot(results2$rho.M, results2$std.CompEst.a, type = "p", 
       xlab = "prop. of var(M) accounted for by m1-m4",
       ylab = "estimate of a path", col = "darkgreen", ylim = c(0, 1), pch = ".", 
       xlim = c(0,1), cex = 2)
  abline (h = .6)  
  points(results2$rho.M, results2$std.LVEst.a, col = results2$col.var.a, pch = ".", cex = 2)
  
  par(fig = c(1/3, 2/3, .5, 1), new = TRUE)
  plot(results2$rho.M, results2$std.CompEst.b, type = "p", 
       xlab = "prop. of var(M) accounted for by m1-m4",
       ylab = "estimate of b path", col = "darkgreen", ylim = c(0, 1), pch = ".", 
       xlim = c(0,1), cex = 2)
  abline (h = .4)  
  points(results2$rho.M, results2$std.LVEst.b, col = results2$col.var.a, pch = ".", cex = 2) 
  legend(x = 0, y = 1, c("composite model", "reflective model"), 
         fill = c("darkgreen", "darkred"), bty = "n")
  
  par(fig = c(2/3, 1, .5, 1), new = TRUE)
  plot(results2$rho.M, results2$std.CompEst.ab, type = "p", 
       xlab = "prop. of var(M) accounted for by m1-m4", xlim = c(0,1),
       ylab = "estimate of a*b path", col = "darkgreen", ylim = c(0, 1), pch = ".", cex = 2)
  abline (h = .24)  
  points(results2$rho.M, results2$std.LVEst.ab, col = results2$col.var.a, pch = ".", cex = 2)

  #add fit
  par(fig = c(0, 1, 0, .5), new = TRUE)
  plot(results2$rho.M, results2$pop.RMSEA, type = "p", xlab = "prop. of var(M) accounted for by m1-m4",
       ylab = "population RMSEA", col = results2$col.var.a, ylim = c(0, 1), xlim = c(0, 1),
       pch = ".", cex = 2)
  points(results2$rho.M, results2$pop.RMSEA.Comp, col = "darkgreen", pch = ".", cex = 2)
dev.off()  

