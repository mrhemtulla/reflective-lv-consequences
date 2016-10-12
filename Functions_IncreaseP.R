

generateCov <- function(a = .6, b = .4, m = 4, avec = NA, covm = NA, res.M = 0){
  #a: total influence of X on M (composite)
  #b: regression of Y on M
  #a1 - a4 are regression paths of causal indicators of M on X:
  #res.M is the residual variance of M (i.e., if it is more than just the sum of m1-m4)
  amat <- avec %*% t(avec)
  
  #if covm is not provided, then solve for gamma and covm
  if (is.na(covm) == TRUE) {
    gamma <- a/sum(avec)
    covm <- (1 - m*gamma^2) / (2*choose(m,2)*gamma^2)    
  } else {
    #if a is not provided, then solve for gamma and a
    gamma <- 1/(2*choose(m,2)*covm + m)
    a <- gamma*(sum(avec)) 
  }
  
  #variance of residuals for the causal indicators of M:
  evec <- 1 - avec^2
  
  #residual covariances among m1-m4: 
  emat <- covm - amat  
  
  Beta <- matrix(0, 3+m, 3+m)
  Beta[1:m,3+m] <- avec
  Beta[m+1,1:m] <- gamma
  Beta[m+2,m+1] <- b
  I <- diag(3+m)
  Psi <- matrix(0, 3+m, 3+m)
  Psi[1:m, 1:m] <- emat
  diag(Psi) <- c(evec, res.M, (1-b^2), 1)  
  
  Sigma <- solve(I - Beta) %*% Psi %*% t(solve(I - Beta))
  colnames(Sigma) <- c(paste("m", 1:m, sep = ""), "M", "Y", "X")
  return(Sigma)  
}

write.LV <- function(m) {  #creates lavaan syntax for the LV model given p
  LVmod <- "LM =~ "
  for (i in 1:(m-1)) {
    LVmod <- paste(LVmod, 'm', i, '+', sep = "")    
  }
  LVmod <- paste(LVmod, 'm', m, "\n LM ~ X \n Y ~ LM", sep = "")
  return(LVmod)
}


CompMod <- 'M ~ X
Y ~ M'

CompMod <- 'FM <~ m1 + m2 + m3 + m4
            FM ~ X
            Y ~ FM'

m <- 7
LVMod <- write.LV(m)

#avec <- rep(.1, m)
avec <- runif(m, -.8, .8)
sigma <- generateCov(a = .6, b = .4, m = m, avec = avec, covm = NA, res.M = 0)

fit.comp <- sem(model = CompMod, sample.cov = sigma, sample.nobs = 200)
summary(fit.comp, standardized = TRUE)

fit.LV <- sem(model = LVMod, sample.cov = sigma, sample.nobs = 200)
summary(fit.LV, standardized = TRUE)


row <- 1
results <- matrix(NA, 27, 16)
for (m in 3:30){
  LVMod <- write.LV(m)
  result <- matrix(NA, 1000, 7)
  for (rep in 1:1000){
    avec <- runif(m, -.8, .8)

    sigma <- matrix(0, m+3, m+3)
    while(min(eigen(sigma[-(m+3),-(m+3)])$values) <= 0){
      avec <- runif(m, 0, .8)
      sigma <- generateCov(a = .6, b = .4, m = m, avec = avec, covm + NA, res.M = 0) 
    }  

    fitComp <- sem(model = CompMod, sample.cov = sigma, sample.nobs = 200)
    fitLV <- try(sem(model = LVMod, sample.cov = sigma, sample.nobs = 200), silent = TRUE)    
    conv.LV <- try(inspect(fitLV, "converged"), silent = TRUE)
    SEs.LV <- TRUE  
    if (class(try(vcov(fitLV), silent = TRUE)) == "try-error") {SEs.LV <- FALSE}
    if (conv.LV == FALSE) {fit.LV <- NA} else (fit.LV <- inspect(fitLV, "fit")["fmin"])  

    comp.a <- inspect(fitComp, "coef")$beta[1,3] #a path composite model
    comp.b <- inspect(fitComp, "coef")$beta[2,1] #b path composite model  
    LV.a <- inspect(fitLV, "coef")$beta[1,3] #a path LV model
    LV.b <- inspect(fitLV, "coef")$beta[2,1] #b path LV model
    result[rep,] <- c(comp.a, comp.b, LV.a, LV.b, SEs.LV, conv.LV, fit.LV)  
  }  
  result <- as.data.frame(result)
  colnames(result) <- c("comp.a", "comp.b", "LV.a", "LV.b", "SEs.LV", "conv.LV", "fit.LV")
  result2 <- result[result$conv.LV == 1,]
  result2 <- result2[result2$SEs.LV == 1,]
  
  results[row,1] <- m
  results[row, 2] <- nrow(result2)   
  results[row, 3:9] <- apply(result2, 2, mean)
  results[row, 10:16] <- apply(result2, 2, var)
  
  row <- row+1
  print(row)
}

results <- data.frame(results)
colnames(results) <- c("no.indicators", "convergence", "mean.comp.a", "mean.comp.b", 
                       "mean.LV.a", "mean.LV.b", "mean.SEs.LV", "mean.conv.LV", 
                       "mean.fit.LV", "var.comp.a", "var.comp.b", 
                       "var.LV.a", "var.LV.b", "var.SEs.LV", "var.conv.LV", "var.fit.LV")

#save(results, file = "Results_NoIndicators.Rdata")


plot(results$no.indicators, results$mean.LV.a, type = "l", ylim = c(0, 1), col = "red4")
lines(results$no.indicators, results$mean.LV.a + 1.96*sqrt(results$var.LV.a), col = "red4")
lines(results$no.indicators, results$mean.LV.a - 1.96*sqrt(results$var.LV.a), col = "red4")
abline(h = .6, col = "darkgreen")


plot(results$no.indicators, results$mean.LV.b, type = "l", ylim = c(0, 1), col = "red4")
lines(results$no.indicators, results$mean.LV.b + 1.96*sqrt(results$var.LV.a), col = "red4")
lines(results$no.indicators, results$mean.LV.b - 1.96*sqrt(results$var.LV.a), col = "red4")
abline(h = .4, col = "blue3")

