match_label = function(Yhat,Y){
  KK = length(unique(Y))
  I = diag(KK)
  perms = permn(KK)
  Label_Yhat_rotate = matrix(0,length(Yhat)[1],factorial(KK))
  match_Yhat = rep(0,factorial(KK))
  Yhat_idx = matrix(0,length(Yhat)[1],KK)
  for(k in unique(Yhat)){
    Yhat_idx[which(Yhat==k),k]=1
  }
  for(i in 1:factorial(KK)){
    perm_vec = perms[[i]]
    Mperm = diag(length(perm_vec))[perm_vec,]
    Label_Yhat_rotate[,i] = max.col(Yhat_idx%*%Mperm,"first")
    match_Yhat[i] = mean(Label_Yhat_rotate[,i]-Y==0)
  }
  max_match = max(match_Yhat)
  rotate = which(match_Yhat==max_match)
  Yhat_match = Label_Yhat_rotate[,rotate[1]]
  perm_vec = perms[[rotate[1]]]
  Qperm = diag(length(perm_vec))[perm_vec,]
  return(list(Yhat_match = Yhat_match,prop_match = max_match,Qperm = Qperm))
}

# DR method for ATE

dr_manual <- function(data, propensity_formula, outcome_formula){
  e_model <- glm(propensity_formula, data = data, family = binomial)
  m0_model <- suppressWarnings(glm(outcome_formula, subset = (Z == 0), data = data, family = binomial))
  m1_model <- suppressWarnings(glm(outcome_formula, subset = (Z == 1), data = data, family = binomial))
  e <- predict(e_model, newdata = data,type = "response")
  m0 <- predict(m0_model,newdata = data, type="response")
  m1 <- predict(m1_model,newdata = data, type="response")
  
  mu1 <- with(data, mean((Z * Y - (Z - e) * m1) / e))
  mu0 <- with(data, mean(((1 - Z) * Y + (Z - e) * m0) / (1 - e)))
  dr.est <-  mu1-mu0
  
  
  I <- (with(data, ((Z * Y - (Z - e) * m1) / e -((1 - Z) * Y + (Z - e) * m0) / (1 - e)))-dr.est)
  dr.I.var <- sum(I^2)/nrow(data)^2
  
  I.mu1 <- (with(data, ((Z * Y - (Z - e) * m1) / e ))-mu1)
  I.mu0 <- (with(data, (((1 - Z) * Y + (Z - e) * m0) / (1 - e)))-mu0)
  cov.mu1.mu0 <- rbind(I.mu1,I.mu0)%*%cbind(I.mu1,I.mu0)/nrow(data)^2
  dr.rr <- mu1/mu0
  dr.var.rr = cov.mu1.mu0[1,1]/(mu0)^2-2*mu1*cov.mu1.mu0[1,2]/(mu0^3)+mu1^2*cov.mu1.mu0[2,2]/(mu0^4)
  
  return(list(estimate = dr.est, variance = dr.I.var, I = I, mu1=mu1, mu0=mu0, 
              cov.mu1.mu0 = cov.mu1.mu0,
              I.mu1=I.mu1, I.mu0=I.mu0, dr.rr = dr.rr, dr.var.rr = dr.var.rr))
}


# Estimation of ATE with population standardization

drTL_manual <- function(data, propensity_formula, density_ratio_formula, outcome_formula){
  Z <- data$Z
  S <- data$S
  Y <- data$Y
  e_model <- glm(propensity_formula, subset = (S == 1), data = data, family = binomial)
  dr_model <- glm(density_ratio_formula, data = data, family = binomial)
  m0_model <- suppressWarnings(glm(outcome_formula, subset = (Z == 0 & S == 0), data = data, family = binomial))
  m1_model <- suppressWarnings(glm(outcome_formula, subset = (Z == 1 & S == 0), data = data, family = binomial))
  
  m0 <- predict(m0_model,newdata = data, type="response")
  m1 <- predict(m1_model,newdata = data, type="response")
  e <- predict(e_model, type = "response", newdata = data)
  dr <- predict(dr_model, type = "response", newdata = data)
  
  r <- (1-dr)/dr
  r[r>10] <- 10
  mu1 <- sum(m1*(1-S)+Z*(Y-m1)*S*r/e)/sum(1-S)
  mu0 <- sum(m0*(1-S)+(1-Z)*(Y-m0)*S*r/(1 - e))/sum(1-S)
  dr.est <- sum(m1*(1-S)-m0*(1-S)+Z*(Y-m1)*S*r/e-(1-Z)*(Y-m0)*S*r/(1 - e))/sum(1-S)
  
  pS0 <- sum(1-S)/nrow(data)
  I <- m1*(1-S)-m0*(1-S)+Z*(Y-m1)*S*r/e-(1-Z)*(Y-m0)*S*r/(1 - e) - dr.est*(1-S)
  dr.I.var <- sum((I)^2)/pS0^2/nrow(data)^2
  
  I.mu1 <- m1*(1-S)+Z*(Y-m1)*S*r/e - mu1*(1-S)
  I.mu0 <- m0*(1-S)+(1-Z)*(Y-m0)*S*r/(1 - e) - mu0*(1-S)
  cov.mu1.mu0 <- rbind(I.mu1,I.mu0)%*%cbind(I.mu1,I.mu0)/nrow(data)^2
  dr.rr <- mu1/mu0
  dr.var.rr = cov.mu1.mu0[1,1]/(mu0)^2-2*mu1*cov.mu1.mu0[1,2]/(mu0^3)+mu1^2*cov.mu1.mu0[2,2]/(mu0^4)
  
  
  return(list(estimate = dr.est, variance = dr.I.var, I = I, mu1=mu1, mu0=mu0, 
              cov.mu1.mu0 = cov.mu1.mu0,
              I.mu1=I.mu1, I.mu0=I.mu0, dr.rr = dr.rr, dr.var.rr = dr.var.rr))
}

# Adaptive aggregation

drAgg <- function(data, propensity_formula, density_ratio_formula, outcome_formula, lambda.lower.bound = NULL){
  target = data[data$S==0,]
  dr.target = dr_manual(target, propensity_formula, outcome_formula)
  
  drTL_source = function(site.idx){
    source = data[data$S==site.idx,]
    source$S =1
    TLdata = rbind(target, source)
    
    dr = drTL_manual(TLdata, propensity_formula, density_ratio_formula, outcome_formula)
    I = rep(0,nrow(data))
    I[data$S==0] = dr$I[TLdata$S==0]
    I[data$S==site.idx] = dr$I[TLdata$S==1]
    
    I.mu1 = rep(0,nrow(data))
    I.mu1[data$S==0] = dr$I.mu1[TLdata$S==0]
    I.mu1[data$S==site.idx] = dr$I.mu1[TLdata$S==1]
    
    I.mu0 = rep(0,nrow(data))
    I.mu0[data$S==0] = dr$I.mu0[TLdata$S==0]
    I.mu0[data$S==site.idx] = dr$I.mu0[TLdata$S==1]
    
    return(list(estimate = dr$estimate, variance = dr$variance, I = I, 
                mu1=dr$mu1, mu0=dr$mu0, I.mu1=I.mu1, I.mu0=I.mu0, rr = dr$dr.rr, rr.variance = dr$dr.var.rr))
  }
  
  dr.sources = lapply(1:max(data$S),drTL_source)
  
  pS0 = sum(data$S==0)/nrow(data)
  I.T = rep(0,nrow(data))
  I.T[data$S==0] = dr.target$I/pS0
  I.y = I.T
  I.X = NULL
  dr.S.est = NULL
  dr.S.mu1 = NULL
  dr.S.mu0 = NULL
  dr.S.rr = NULL
  dr.S.rr.var = NULL
  for(i in 1:max(data$S)){
    I.x = I.T+dr.target$estimate-dr.sources[[i]]$I/pS0-dr.sources[[i]]$estimate
    I.X = cbind(I.X,I.x)
    dr.S.est = c(dr.S.est,dr.sources[[i]]$estimate)
    dr.S.mu1 = c(dr.S.mu1, dr.sources[[i]]$mu1)
    dr.S.mu0 = c(dr.S.mu0, dr.sources[[i]]$mu0)
    dr.S.rr = c(dr.S.rr, dr.sources[[i]]$rr)
    dr.S.rr.var = c(dr.S.rr.var, dr.sources[[i]]$rr.variance)
  }
  if(max(data$S)==1){
    eta.S <- coef(lm(I.y~I.X-1))
    eta = c(1-sum(eta.S),eta.S)
    lambda = NA
  }else{
    I.data = data.frame(I.y,I.X)
    #penalty.weight = -log(abs(dr.S.est-dr.target$estimate))
    #cv_model = cv.glmnet(I.X, I.y, alpha = 1, intercept = FALSE, penalty.factor = penalty.weight)
    penalty.weight = (abs(dr.S.est-dr.target$estimate))
    #cv_model = cv.glmnet(I.X, I.y, alpha = 1, intercept = FALSE, penalty.factor = penalty.weight, lambda = seq(0.05,0.1,by = 0.001))
    cv_model = cv.glmnet(I.X, I.y, alpha = 1, intercept = FALSE, penalty.factor = penalty.weight, type.measure="mse")
    lambda = max(cv_model$lambda.1se, lambda.lower.bound)
    fit_model = glmnet(I.X, I.y, alpha = 1, intercept = FALSE, lambda = lambda, penalty.factor = penalty.weight)
    eta.S = c(coefficients(fit_model)[1+1:max(data$S)])
    eta = c(1-sum(eta.S),eta.S)
    
  }
  
  dr.est = sum(c(dr.target$estimate,dr.S.est)*eta)
  dr.mu1  = sum(c(dr.target$mu1,dr.S.mu1)*eta)
  dr.mu0  = sum(c(dr.target$mu0,dr.S.mu0)*eta)
  
  
  
  cal.joint.Var = function(dr.target,res.TATE.Han.Agg, data){
    pS0 = sum(data$S==0)/nrow(data)
    I.T = rep(0,nrow(data))
    I.T[data$S==0] = dr.target$I/pS0
    
    I.S = NULL
    for(i in 1:length(dr.sources)){
      I.s = dr.sources[[i]]$I/pS0
      I.S = cbind(I.S,I.s)
    }
    I.all = cbind(I.T,I.S)
    return((t(I.all)%*%I.all)/(nrow(data))^2)
  }
  joint.var = cal.joint.Var(dr.target,dr.sources, data) 
  dr.var = eta%*%joint.var%*%eta
  
  
  cal.joint.Var.mu = function(dr.target,res.TATE.Han.Agg, data){
    pS0 = sum(data$S==0)/nrow(data)
    I.T.mu1 =I.T.mu0 = rep(0,nrow(data))
    I.T.mu1[data$S==0] = dr.target$I.mu1/pS0
    I.T.mu0[data$S==0] = dr.target$I.mu0/pS0
    
    
    I.S.mu1 = I.S.mu0 = NULL
    for(i in 1:length(dr.sources)){
      
      I.s1 = dr.sources[[i]]$I.mu1/pS0
      I.s0 = dr.sources[[i]]$I.mu0/pS0
      I.S.mu1 = cbind(I.S.mu1,I.s1)
      I.S.mu0 = cbind(I.S.mu0,I.s0)
    }
    I.all.mu1 = cbind(I.T.mu1,I.S.mu1)
    I.all.mu0 = cbind(I.T.mu0,I.S.mu0)
    var.mu1 = (t(I.all.mu1)%*%I.all.mu1)/(nrow(data))^2
    var.mu0 = (t(I.all.mu0)%*%I.all.mu0)/(nrow(data))^2
    cov.mu1.mu0 <- (t(I.all.mu1)%*%I.all.mu0)/(nrow(data))^2
    
    return(list(var.mu1,var.mu0, cov.mu1.mu0))
  }
  joint.var.mu = cal.joint.Var.mu(dr.target,dr.sources, data) 
  dr.var.mu1 = eta%*%joint.var.mu[[1]]%*%eta
  dr.var.mu0 = eta%*%joint.var.mu[[2]]%*%eta
  
  dr.cov.mu1.mu0 = eta%*%joint.var.mu[[3]]%*%eta
  dr.var.rr = dr.var.mu1/(dr.mu0)^2-2*dr.mu1*dr.cov.mu1.mu0/(dr.mu0^3)+dr.mu1^2*dr.var.mu0/(dr.mu0^4)
  #dr.rr = sum(eta*c(dr.target$dr.rr,dr.S.rr))
  #dr.var.rr = sum(eta^2*c(dr.target$dr.var.rr, dr.S.rr.var))
  
  return(list(estimate = dr.est, variance = dr.var, eta = eta, lambda = lambda, 
              mu1.est = dr.mu1, mu1.variance = dr.var.mu1, 
              mu0.est =dr.mu0, mu0.variance = dr.var.mu0,
              rr = dr.mu1/dr.mu0, rr.variance = dr.var.rr))
}

calATE = function(k, c, KSiteData){
  SiteData = KSiteData[KSiteData$S==k,]
  SiteData.c = SiteData[SiteData$class.est==c,]
  
  dr <- dr_manual(SiteData.c, propensity_formula , outcome_formula)
  return(data.frame(site = k,
                    estimate = c(dr$estimate, dr$dr.rr),
                    var = c(dr$variance, dr$dr.var.rr),
                    ll = c(dr$estimate-1.96*sqrt(dr$variance),
                           dr$dr.rr-1.96*sqrt(dr$dr.var.rr)),
                    ul = c(dr$estimate+1.96*sqrt(dr$variance),
                           dr$dr.rr+1.96*sqrt(dr$dr.var.rr)),
                    type = c("Risk difference","Risk Ratio")
  ))
}


calTATE.agg.CI <- function(k.target, source.site, c, KSiteData, lambda.lower.bound = NULL){
  CI.target <- calATE(k.target, c=c, KSiteData)
  
  TargetData = KSiteData[KSiteData$S==k.target,]
  TargetData.c = TargetData[TargetData$class.est==c,]
  TargetData.c$S = 0
  
  source = NULL
  for(i in 1:length(source.site)){
    SourceData = KSiteData[KSiteData$S==source.site[i],]
    SourceData.c = SourceData[SourceData$class.est==c,]
    SourceData.c$S = i
    source = rbind(source, SourceData.c)
  }
  dt = rbind(TargetData.c, source)
  res.agg <- drAgg(dt, propensity_formula, density_ratio_formula, outcome_formula, lambda.lower.bound)  
  
  weights.agg <- data.frame(site = c(k.target,source.site),weights = round(res.agg$eta,4))
  
  rd <- res.agg$estimate
  rd.var <- res.agg$variance
  rr <- res.agg$rr
  rr.var <- res.agg$rr.variance
  mu1 <- res.agg$mu1.est
  mu1.var <- res.agg$mu1.variance
  mu0 <- res.agg$mu0.est
  mu0.var <- res.agg$mu0.variance
  
  CI.target$site <- NULL
  CI.target$type <- c("Target-only risk difference", "Target-only risk ratio")
  
  CI.all <- rbind(CI.target,
                  data.frame(type = "Risk difference",
                             estimate = rd, var = rd.var,
                             ll = rd-1.96*sqrt(rd.var), 
                             ul = rd+1.96*sqrt(rd.var)),
                  data.frame(type = "Risk ratio",
                             estimate = rr, var = rr.var,
                             ll = rr-1.96*sqrt(rr.var), 
                             ul = rr+1.96*sqrt(rr.var)),
                  data.frame(type = "Positive mean",
                             estimate = mu1, var = mu1.var,
                             ll = mu1-1.96*sqrt(mu1.var), 
                             ul = mu1+1.96*sqrt(mu1.var)),
                  data.frame(type = "Negative mean",
                             estimate = mu0, var = mu0.var,
                             ll = mu0-1.96*sqrt(mu0.var), 
                             ul = mu0+1.96*sqrt(mu0.var))
  )
  
  
  
  return(list(CI=CI.all,weights = weights.agg, lambda = res.agg$lambda))
}
