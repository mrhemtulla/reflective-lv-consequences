#Study 4: MIMIC model

rm(list = ls())
setwd("C:/Users/Mijke/Documents/__STUDIES__/CausalIndicators/MediationModelSimulation")
source("Functions.R")
windowsFonts(A=windowsFont("Verdana"))

######################## SIM 4.1 ###################################################

sim4.1 <- function(i){
  a <- .6 #increase this to increase the range of possible values.  
  b <- .4
  res.M <- .25
  
  data <- matrix(1, 7, 7)
  while(min(eigen(data[-c(5,8),-c(5,8)])$values) < 1e-10){
    covm <- -.1
    while(covm < 0){ #repeat until covm is positive
      avec <- runif(2, .3, .8) 
      gamma <- a/sum(avec)
      covm <- (1 - 2*gamma^2 - res.M)/(2*gamma^2)      
    }
    lamvec <- runif(2, .3, 1) 
    data <- generateCov.mimic(a = a, b = b, avec = avec, lamvec = lamvec, res.M = res.M, covm = NA) 
  }
  
  fitComp <- try(sem(model = CompMod, sample.cov = data, sample.nobs = 100000, warn = FALSE), silent = TRUE)
  fitLV <- try(sem(model = LVMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE, warn = FALSE), silent = TRUE)  
  #fitMimic <- try(sem(model = MimicMod, sample.cov = data, sample.nobs = 100000, std.lv = TRUE, warn = FALSE), silent = TRUE) 
  
  conv.LV <- try(inspect(fitLV, "converged"), silent = TRUE)
  conv.Comp <- try(inspect(fitComp, "converged"), silent = TRUE)
  SEs.LV <- TRUE  
  SEs.Comp <- TRUE
  
  if (identical(class(try(vcov(fitLV), silent = TRUE)), "try-error") == TRUE) {SEs.LV <- FALSE}
  if (identical(class(try(vcov(fitComp), silent = TRUE)), "try-error") == TRUE) {SEs.Comp <- FALSE}
  if (conv.LV == FALSE | conv.Comp == FALSE) {
    return(c(avec, lamvec, gamma, covm, rep(NA, 10), 0, NA, 0, NA))
  } else {
    return(c(avec, lamvec, gamma, covm,
             inspect(fitLV, "std")$lambda[1:4,1],
             inspect(fitComp, "std")$beta[1,3], #a path composite model standardized
             inspect(fitComp, "std")$beta[2,1], #b path composite model standardized 
             inspect(fitLV, "std")$beta[1,3], #a path LV model standardized
             inspect(fitLV, "std")$beta[2,1], #b path LV model standardized
             inspect(fitComp, "fit")["fmin"], #fit of composite model
             inspect(fitLV, "fit")["fmin"], conv.LV, SEs.LV, conv.Comp, SEs.Comp))
  }
}


X <- 1:10000
ptm <- proc.time()
out <- sapply(X, sim4.1)
proc.time() - ptm

save(out, file = "sim4.1.Rdata")
load(file = "sim4.1.Rdata")
results <- as.data.frame(t(out))

colnames(results) <- c("a1", "a2", "lam1", "lam2", "gamma", "covm",   
                       "lambda1", "lambda2", "lambda3", "lambda4",
                       "CompEst.a", "CompEst.b", "LVEst.a", "LVEst.b", 
                       "Comp.fmin", "LV.fmin", "conv.LV", "SEconv.LV",
                       "conv.Comp", "SEconv.Comp")

#discard rows that did not converge: 
results2 <- results[is.na(results$LV.fmin) == FALSE,]
results2 <- results2[results2$SEconv.LV == 1,]
dim(results2) #9997

#remove rows with negative residual variances? 
results2 <- results2[results2$lambda1 < 1,]
results2 <- results2[results2$lambda2 < 1,]
results2 <- results2[results2$lambda3 < 1,]
results2 <- results2[results2$lambda4 < 1,]
dim(results2) #9427

#compute covariance among m variables: 

results2$mean.cov <- apply(cbind(results2$covm, results2$gamma*results2$lam1, 
                                 results2$gamma*results2$lam1, results2$gamma*results2$lam2,
                                 results2$gamma*results2$lam2, results2$lam1*results2$lam2), 1, mean)
results2$var.cov <- apply(cbind(results2$covm, results2$gamma*results2$lam1, 
                                 results2$gamma*results2$lam1, results2$gamma*results2$lam2,
                                 results2$gamma*results2$lam2, results2$lam1*results2$lam2), 1, var)

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
rgPal <- colorRampPalette(c('darkgreen','yellow'))

results2$col.var.cov <- rbPal(25)[as.numeric(cut(results2$var.cov,breaks = 25))]
results2$col2.var.cov <- rgPal(25)[as.numeric(cut(results2$var.cov,breaks = 25))]
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
jpeg("sim4.1.Plot1.jpg", width = 8, height = 8, units = "in", res = 800) 
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, .5, .5, 1))
plot(results2$covm, results2$CompEst.a, type = "p", xlab = "cov between causal indicators",
     ylab = "estimate of a path", col = results2$col2.var.cov, ylim = c(.5, 1), pch = ".", 
     cex = 1, xlim = c(0, 1))
points(results2$covm, results2$LVEst.a, col = results2$col.var.cov, pch = ".", cex = 1)
legend(x = .6, y = 1, c("composite model", "reflective model"), 
       fill = c("green4", "darkred"), bty = "n")
abline(h = .6, col = "black")

par(fig = c(.5, 1, .5, 1), new = TRUE)
plot(results2$covm, results2$CompEst.b, type = "p", xlab = "cov between causal indicators",
     ylab = "estimate of b path", col = results2$col2.var.cov, ylim = c(0, .5), pch = ".", 
     cex = 1, xlim = c(0, 1))
points(results2$covm, results2$LVEst.b, col = results2$col.var.cov, pch = ".", cex = 1)
abline(h = .4)

par(fig = c(0, .5, 0, .5), new = TRUE)
plot(results2$covm, results2$LVEst.ab, type = "p", 
     xlab = "cov between causal indicators", ylab = "estimate of a*b path", col = results2$col.var.cov, 
     ylim = c(.2, .4), pch = ".", cex = 1, xlim = c(0, 1))
abline(h = .24)

par(fig = c(.5, 1, 0, .5), new = TRUE)
plot(results2$covm, results2$CompEst.ab, type = "p", 
     xlab = "cov between causal indicators", ylab = "estimate of a*b path", col = results2$col2.var.cov, 
     ylim = c(.2, .4), pch = ".", cex = 1, xlim = c(0, 1))
abline(h = .24)

dev.off()  

#plot ab estimates only:
jpeg("sim4.1.Plot2.jpg", width = 8, height = 4, units = "in", res = 800) 
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, .5, 0, 1))
plot(results2$covm, results2$CompEst.ab, type = "p", main = "Composite Model",
     xlab = "cov between causal indicators", ylab = "estimate of a*b path", 
     col = results2$col.var.cov, 
     ylim = c(.2, .4), pch = ".", cex = 1, xlim = c(0, 1))
abline(h = .24)

par(fig = c(.5, 1, 0, 1), new = TRUE)
plot(results2$covm, results2$LVEst.ab, type = "p", main = "Reflective Model", 
     xlab = "cov between causal indicators", ylab = "estimate of a*b path", 
     col = results2$col.var.cov, #col = results2$col.a.,
     ylim = c(.2, .4), pch = ".", cex = 1, xlim = c(0, 1))
abline(h = .24)
dev.off()  


jpeg("sim4.1.Plot2b.jpg", width = 8, height = 4, units = "in", res = 800) 
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
jpeg("sim4.1.Plot3.jpg", width = 9, height = 4, units = "in", res = 800) 
par(mar = c(5, 5, 1,.5), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
plot(results2$covm, results2$pop.RMSEA.LV, type = "p", xlab = "cov between causal indicators",
     ylab = "population RMSEA", col = results2$col.var.cov, ylim = c(0, .6), 
     pch = ".", cex = 1)
points(results2$covm, results2$pop.RMSEA.comp, col = results2$col2.var.cov, pch = ".", cex = 1)
legend(x = .4, y = .6, c("composite model", "reflective model"), 
       fill = c("green4", "darkred"), bty = "n")
dev.off()  


#correlate fit and bias: 
round(cor(results2[c("LVEst.ab", "CompEst.ab", "pop.RMSEA.LV", "pop.RMSEA.comp")]), 2)
