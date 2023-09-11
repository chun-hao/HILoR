HILoR$set("private", "logistic", function(y, X, int_loc, lambda, 
                                          beta=NULL, lb=10, tol=10^-6, maxite=20){
    p <- ncol(X)
    n <- length(y)
    
    if(is.null(beta)){beta <- rep(0, p); beta[1] <- 1}
    
    for(i in 1:maxite){
        phat <- c(X%*%beta)
        phat <- ifelse(phat > lb, lb, ifelse(phat < -lb, -lb, phat)) 
        phat <- exp(phat)/(1+exp(phat))
        vhat <- phat*(1-phat)
        b_penalty <- beta; b_penalty[int_loc] <- 0
        beta_new <- beta + solve(t(X*vhat)%*%X+lambda*diag(p), t(X)%*%(y-phat)-lambda*beta)
        beta_new[is.na(beta_new)] <- beta[is.na(beta_new)] # add
        test <- max(abs(beta-beta_new))
        beta <- beta_new
        if(test < tol) break
    }
    return(beta)
})

HILoR$set("private", "p_vec", function(Data, beta){
    k <- length(beta)
    
    # event probability p_ij
    P <- sapply(1:k, function(i) exp(Data[[i+1]]%*%beta[[i]]))
    P <- P/(1+P)
    P <- ifelse(P < 1e-6, 1e-6, ifelse(P > 1-1e-6, 1-1e-6, P)) # add
    
    # bag probability pi_ij
    t1 <- apply(P, 1, function(x) prod(1-x))
    BP <- 1 - t1
    
    # bag probability pi^[2]_ij
    t2 <- apply(P, 1, function(x) sum(x/(1-x)))
    BP2 <- 1 - (1 + t2) * t1
    
    return(list(EP=P, BP=BP, BP2=BP2))
})

HILoR$set("private", "Var_MILReg", function(beta, Data, opt.type, lambda){
    k <- length(beta)
    oo <- private$p_vec(Data, beta)
    p_v <- oo$EP
    
    pi_1_v <- oo$BP
    pi_2_v <- oo$BP2
    
    if(opt.type == 1){
    X <- Data[[2]]*p_v[,1]*(1-pi_1_v)
    for(inner in 2:k) X <- cbind(X, Data[[inner+1]]*p_v[,inner]*(1-pi_1_v))
    V <- 1/(1-pi_1_v)/pi_1_v
    }
    if(opt.type == 2){
    X <- Data[[2]]*p_v[,1]*(1-pi_2_v-(1-pi_1_v)/(1-p_v[,1]))
    for(inner in 2:k) X <- cbind(X, Data[[inner+1]]*p_v[,inner]**(1-pi_2_v-(1-pi_1_v)/(1-p_v[,inner])))
    V <- 1/(1-pi_2_v)/pi_2_v
    }
    
    
    A <- solve(t(X*V)%*%X + diag(ncol(X))*lambda)
    return(A)
})

HILoR$set("private", "gamma_vec", function(Data, beta, type){
    
        temp <- private$p_vec(Data, beta)
        out <- temp$EP/temp$BP
        
        if(type == 2){
          omega <- apply(temp$EP, 1, function(x) sum(x/(1-x)))
          omega <- temp$EP/(1-temp$EP)/(1+omega)
          out <- (temp$EP-omega)/temp$BP2
          return(list(out=out, omega=omega))
        }
        
        return(out)
    })

HILoR$set("private", "MILReg_Inner", function(Data, lambda, beta=NULL, 
                                              tol, maxite=2000, type){
    
    k <- length(Data)
    Y <- Data[[1]]
    
    qq <- length(Data)
    int_loc <- NULL
    for(i in 2:length(Data)) int_loc[i-1] <- ncol(Data[[i]])
    int_loc <- c(0, cumsum(int_loc)) + 1
    int_loc <- int_loc[-length(int_loc)]
    
    # set initial values for beta
    if(is.null(beta)){
    beta <- NULL; beta <- as.list(beta)
    for(inner in 1:(k-1)) beta[[inner]] <- glm(Data[[1]]~-1+Data[[inner+1]], family=binomial) %>% coefficients
    }
    
    beta_new <- beta
    for(cc in 1:maxite){
    for(inner in 1:(k-1)){
      oo <- private$gamma_vec(Data, beta_new, type)
      if(type == 1) beta_new[[inner]] <- private$logistic(Data[[1]]*oo[,inner], Data[[inner+1]], int_loc, lambda, beta=NULL)
      if(type == 2) beta_new[[inner]] <- private$logistic(oo$omega[,inner] + Data[[1]]*oo$out[,inner], 
                                                  Data[[inner+1]], int_loc, lambda, beta=NULL)
    }
    
    test <- max(abs(unlist(beta)-unlist(beta_new)))
    beta <- beta_new
    #cat(cc, test, "\r")
    if(test < tol) break
    }
    
    PV <- private$p_vec(Data, beta)
    if(type == 1) pi_vec <- PV$BP
    if(type == 2) pi_vec <- PV$BP2
    logLik <- sum(Data[[1]]*log(pi_vec/(1-pi_vec))+log(1-pi_vec))
    
    return(list(beta=beta,ite=cc, logLik=logLik, EP=PV$EP, lambda=lambda))
})

HILoR$set("public", "fit", function(lambda = 0, maxite = 2000, tol = 10^-6, type = 1){
    stopifnot(!is.null(self$bag_label), !is.null(self$X))
    
    Data <- append(self$X, list(self$bag_label), 0)
    for(i in 1:self$n_I){
        Data[[i+1]] <- cbind(1,Data[[i+1]])
    }
    
    oo <- private$MILReg_Inner(Data, lambda, beta=NULL, maxite=maxite, tol=tol, type)
    V <- private$Var_MILReg(oo$beta, Data, opt.type=type, lambda)
    vv <- round(sqrt(diag(V)), digits=4)
    
    beta <- oo$beta
    names(beta) <- paste0("beta", 1:self$n_I)
    for(i in 1:self$n_I){
        beta[[i]] <- as.vector(beta[[i]])
        names(beta[[i]]) <- c("(Intercept)", self$var_names[[i]])
    }
    self$coefficients <- beta
    
    beta_vec <- round(unlist(oo$beta), digits=4)
    pv <- round(2*pnorm(-abs(beta_vec/vv)), digits=4)
    sig <- ifelse(pv<.01, "**", ifelse(pv<.05, "*", ""))
    OUT <- data.frame(beta_vec, se=vv, z=round(beta_vec/vv, digits=4), pvalue=pv, sig=sig)
    
    private$fitting_result <- list(beta=beta, se=vv, pv=pv, 
                                   test.table=OUT, EP=oo$EP, logLik= oo$logLik, ite=oo$ite)
})

HILoR$set("public", "summary", function(){
    stopifnot("Fit the model first!" = !is.null(self$coefficients))
    
    cat("Data Information:\n")
    cat("    Number of bags: ", self$n_B, "\n", sep = "")
    cat("    Number of instances per bag: ", self$n_I, "\n", sep = "")
    cat("    Number of covariates for each instance: ", 
        paste(sapply(self$X, ncol), collapse = ", "), "\n", sep = "")
    cat("    The data is of type ", self$type, " (A bag is positive if at least ", 
        self$type, " of the instances are positive)", "\n", sep = "")
    
    cat("\n")
    
    cat("Coefficients:\n\n")
    
    for(i in 1:n_I){
        cat("Instance ", i, ":\n")
        start <- ifelse(i == 1, 1, length(self$coefficients[[i-1]]) + 1)
        end <- start + length(self$coefficients[[i]]) - 1
        tmp <- private$fitting_result$test.table[start:end, ]
        row.names(tmp) <- names(self$coefficients[[i]])
        colnames(tmp) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Sig.")
        print(tmp)
        cat("\n")
    }
    
    cat("===============================\n")
    
    AIC <- -2*private$fitting_result$logLik + 2*length(unlist(self$coefficients))
    BIC <- -2*private$fitting_result$logLik + 
        log(length(self$n_B))*length(unlist(self$coefficients))
    
    cat("AIC: ", AIC, "\n")
    cat("BIC: ", BIC, "\n\n")
    
    cat("Number of iterations: ", private$fitting_result$ite, "\n")
})
