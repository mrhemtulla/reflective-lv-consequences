#do some simulations! 

rm(list = ls())
setwd("C:/Users/Mijke/Documents/__STUDIES__/CausalIndicators/MediationModelSimulation")
source("Functions.R")
windowsFonts(A=windowsFont("Verdana"))

######################## SIM 3.1 ###################################################

sim3.1 <- function(i){
  a <- .6 #increase this to increase the range of possible values.  
  b <- .4
  
  data <- matrix(1, 7, 7)
  while(min(eigen(data[-7,-7])$values) < 0){
    beta <- runif(4, 0, .8) 
    data <- generateCov.net(a = a, b = b, beta = beta) 
  }
  
  meana <- mean(beta)
  vara <- var(beta)
  rangea <- range(beta)[2] - range(beta)[1]
  
  fitComp <- sem(model = CompMod, sample.cov = data, sample.nobs = 100000, warn = FALSE)  
  fitLV <- try(sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE, warn = FALSE), silent = TRUE)  
  fitNet <- sem(model = netMod, sample.cov = data, sample.nobs = 100000, warn = FALSE)
  
  conv.LV <- try(inspect(fitLV, "converged"), silent = TRUE)
  SEs.LV <- TRUE  
  if (identical(class(try(vcov(fitLV), silent = TRUE)), "try-error") == TRUE) {SEs.LV <- FALSE}
  if (conv.LV == FALSE) {fit.LV <- NA} else (fit.LV <- inspect(fitLV, "fit")["fmin"])
  
  return(c(meana, rangea, vara, 
           parameterEstimates(fitNet)$est[5:1], #5 coefficients from true model              
           inspect(fitLV, "std")$lambda[1:4,1],
           inspect(fitComp, "std")$beta[1,3], #a path composite model standardized
           inspect(fitComp, "std")$beta[2,1], #b path composite model standardized 
           inspect(fitLV, "std")$beta[1,3], #a path LV model standardized
           inspect(fitLV, "std")$beta[2,1], #b path LV model standardized
           inspect(fitComp, "fit")["fmin"], #fit of composite model
           fit.LV, conv.LV, SEs.LV))
}

X <- 1:100000
ptm <- proc.time()
out <- sapply(X, sim3.1)
proc.time() - ptm

save(out, file = "sim3.1.Rdata")
load(file = "sim3.1.Rdata")
results <- as.data.frame(t(out))

colnames(results) <- c("mean.a", "range.a", "var.a", "beta1", "beta2", "beta3", "beta4", "b",  
                       "lambda1", "lambda2", "lambda3", "lambda4",
                       "CompEst.a", "CompEst.b", "LVEst.a", "LVEst.b", 
                       "Comp.fmin", "LV.fmin", "conv.LV", "SEconv.LV")

#discard rows that did not converge: 
results2 <- results[is.na(results$LV.fmin) == FALSE,]
results2 <- results2[results2$SEconv.LV == 1,]
dim(results2) #99977

#remove rows with negative residual variances? 
results2 <- results2[results2$lambda1 < 1,]
results2 <- results2[results2$lambda2 < 1,]
results2 <- results2[results2$lambda3 < 1,]
results2 <- results2[results2$lambda4 < 1,]
dim(results2) #78265

#compute covariance among m variables: 
results2$a. <- .6/(results2$beta1*results2$beta3 + results2$beta2*results2$beta4)

results2$mean.cov <- apply(cbind(results2$beta1, results2$beta2, results2$beta3, results2$beta4, 
                    results2$beta1*results2$beta2, 
                    (results2$beta1*results2$beta3 + results2$beta2*results2$beta4)), 1, mean)
results2$var.cov <-  apply(cbind(results2$beta1, results2$beta2, results2$beta3, results2$beta4, 
                                results2$beta1*results2$beta2, 
                                (results2$beta1*results2$beta3 + results2$beta2*results2$beta4)), 1, var)
round(cor(results2[c("a.", "mean.cov", "var.cov", "mean.a", "var.a")]), 2)

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
rgPal <- colorRampPalette(c('darkgreen','yellow'))

results2$col.var.a <- rbPal(25)[as.numeric(cut(results2$var.a,breaks = 25))]
results2$col2.var.a <- rgPal(25)[as.numeric(cut(results2$var.a,breaks = 25))]
results2$col.a. <- rbPal(25)[as.numeric(cut(results2$a., breaks = 25))]
results2$col.mean.cov <- rbPal(25)[as.numeric(cut(results2$var.a,breaks = 25))]
results2$col.var.cov <- rgPal(25)[as.numeric(cut(results2$var.a,breaks = 25))]

#compute indirect effects:
results2$LVEst.ab <- results2$LVEst.a * results2$LVEst.b
results2$CompEst.ab <- results2$CompEst.a * results2$CompEst.b

#compute fit: 
results2$pop.RMSEA.LV <- sqrt(results2$LV.fmin/9)
results2$pop.RMSEA.comp <- sqrt(results2$Comp.fmin/1)

#plot a, b, ab estimates:
#NOTE for this study a is undefined. don't use this plot! 
jpeg("sim3.1.Plot1.jpg", width = 8, height = 8, units = "in", res = 800) 
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, .5, .5, 1))
  plot(results2$a., results2$CompEst.a, type = "p", xlab = "a' value",
       ylab = "estimate of a path", col = results2$col2.var.a, ylim = c(.5, 1), pch = ".", 
       cex = 1, xlim = c(.7, 1))
  points(results2$a., results2$LVEst.a, col = results2$col.var.a, pch = ".", cex = 1)
  legend(x = .7, y = .2, c("composite model", "reflective model"), 
         fill = c("green4", "darkred"), bty = "n")
  abline(h = .6, col = "black")
  
  par(fig = c(.5, 1, .5, 1), new = TRUE)
  plot(results2$a., results2$CompEst.b, type = "p", xlab = "a' value",
       ylab = "estimate of b path", col = results2$col2.var.a, ylim = c(0, .5), pch = ".", 
       cex = 1, xlim = c(.7, 1))
  points(results2$a., results2$LVEst.b, col = results2$col.var.a, pch = ".", cex = 1)
  abline(h = .4)

  par(fig = c(0, .5, 0, .5), new = TRUE)
  plot(results2$a., results2$LVEst.ab, type = "p", 
       xlab = "a' value", ylab = "estimate of a*b path", col = results2$col.var.a, 
       ylim = c(.2, .4), pch = ".", cex = 1, xlim = c(.7, 1))
  abline(h = .24)

  par(fig = c(.5, 1, 0, .5), new = TRUE)
  plot(results2$a., results2$CompEst.ab, type = "p", 
       xlab = "a' value", ylab = "estimate of a*b path", col = results2$col2.var.a, 
       ylim = c(.2, .4), pch = ".", cex = 1, xlim = c(.7, 1))
  abline(h = .24)

dev.off()  

#plot ab estimates only:
jpeg("sim3.1.Plot2.jpg", width = 8, height = 4, units = "in", res = 800) 
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, .5, 0, 1))
  plot(results2$mean.a, results2$CompEst.ab, type = "p", main = "Composite Model",
       xlab = "mean beta value", ylab = "estimate of a*b path", col = results2$col.a., #col = results2$col.a.,  
       ylim = c(.2, .4), pch = ".", cex = 1, xlim = c(.3, .7))
  abline(h = .24)

  par(fig = c(.5, 1, 0, 1), new = TRUE)
  plot(results2$mean.a, results2$LVEst.ab, type = "p", main = "Reflective Model", 
       xlab = "mean beta value", ylab = "estimate of a*b path", col = results2$col.a., #col = results2$col.a.,
       ylim = c(.2, .4), pch = ".", cex = 1, xlim = c(.3, .7))
  abline(h = .24)
dev.off()  


jpeg("sim3.1.Plot2b.jpg", width = 8, height = 4, units = "in", res = 800) 
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, .5, 0, 1))
plot(results2$a., results2$CompEst.ab, type = "p", main = "Composite Model",
     xlab = "a' value", ylab = "estimate of a*b path", col = results2$col.var.a, 
     ylim = c(.2, .4), pch = ".", cex = 1, xlim = c(.7, 1))
abline(h = .24)

par(fig = c(.5, 1, 0, 1), new = TRUE)
plot(results2$a., results2$LVEst.ab, type = "p", main = "Reflective Model", 
     xlab = "a' value", ylab = "estimate of a*b path", col = results2$col.var.a,
     ylim = c(.2, .4), pch = ".", cex = 1, xlim = c(.7, 1))
abline(h = .24)
dev.off()  


#plot fit
jpeg("sim3.1.Plot3.jpg", width = 9, height = 4, units = "in", res = 800) 
  par(mar = c(5, 5, 1,.5), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  plot(results2$a., results2$pop.RMSEA.LV, type = "p", xlab = "a' value",
       ylab = "population RMSEA", col = results2$col.var.a, ylim = c(0, .6), 
       pch = ".", cex = 1)
  points(results2$a., results2$pop.RMSEA.comp, col = results2$col2.var.a, pch = ".", cex = 1)
  legend(x = .4, y = .6, c("composite model", "reflective model"), 
       fill = c("green4", "darkred"), bty = "n")
dev.off()  


#correlate fit and bias: 
round(cor(results2[c("LVEst.ab", "CompEst.ab", "pop.RMSEA.LV", "pop.RMSEA.comp")]), 2)
