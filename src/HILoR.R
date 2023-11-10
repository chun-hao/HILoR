library(R6)
library(tidyverse)
library(data.table)

rm(list = ls())

HILoR <- R6Class("HILoR", list(
    n_B = NA,
    n_I = NA,
    bag_label = NULL,
    X = NULL,
    var_names = NULL,
    ins_label = NULL,
    type = NA,
    coefficients = NULL,
    pval = NULL,
    fitting_result = NULL,
    initialize = function(n_B, n_I) {
        self$n_B <- n_B
        self$n_I <- n_I
    },
    print = function(...) {
        cat("Number of bags: ", self$n_B, "\n", sep = "")
        cat("Number of instances per bag: ", self$n_I, "\n", sep = "")
        if(!is.null(self$X)){
            cat("Number of covariates for each instance: ", 
                paste(sapply(self$X, ncol), collapse = ", "), "\n", sep = "")
        }
        if(!is.na(self$type)){
            cat("The data is of type ", self$type, " (A bag is positive if at least ", 
                self$type, " of the instances are positive)", "\n", sep = "")
        }
        invisible(self)
    }
))


HILoR$set("public", "inputData", function(bag_label, X, type) {
    stopifnot(is.null(dim(bag_label)), 
              "lenght of bag_label must equal to n_B" = length(bag_label) == self$n_B)
    stopifnot("X must be a list of n_I matrices" = is.list(X), 
              length(X) == self$n_I, 
              all(sapply(X, nrow) == self$n_B),
              !(list(NULL) %in% lapply(X, colnames)))
    stopifnot("type must be 1 or 2" = type %in% c(1,2))
    
    self$bag_label <- bag_label
    self$X <- X
    self$var_names <- lapply(X, colnames)
    self$type <- type
    invisible(self)
})

HILoR$set("public", "toLR", function() {
    stopifnot(!is.null(self$X), !is.null(self$bag_label))
    
    out <- self$X[[1]]
    for(i in 2:n_I){
        out <- cbind(out, self$X[[i]])
    }
    out <- out[, !duplicated(colnames(out))]
    return(as.data.frame(cbind(Y = self$bag_label, out)))
})

HILoR$set("public", "toMILR", function(){
    dat <- self$X
    for(i in 1:n_I){
        dat[[i]] <- cbind(i, dat[[i]])
        colnames(dat[[i]]) <- c("ins_number", self$var_names[[i]])
    }
    out <- as_tibble(dat[[1]]) 
    for (i in 2:n_I){
        out <- full_join(out, as_tibble(dat[[i]]))
    }
    out <- out[,-1] %>% replace(is.na(.), 0)
    
    return(list(bag_label = rep(self$bag_label, times = self$n_I), 
                X = as.matrix(out), ID = rep(1:self$n_B, times = self$n_I)))
})

HILoR$set("public", "generate", function(beta, type = 1, seed = 1234){
    stopifnot(is.list(beta), length(beta) == self$n_I,
              !(list(NULL) %in% lapply(beta, names)))
    stopifnot("type must be 1 or 2" = type %in% c(1,2))
    
    self$var_names <- lapply(beta, function(x) names(x)[-1])
    self$type <- type
    
    
    set.seed(seed)
    X <- vector("list", length = self$n_I)
    all_var_name <- unique(unlist(self$var_names))
    p_total <- length(all_var_name)
    X_tmp <- matrix(rnorm(p_total*self$n_B), nrow = self$n_B, ncol = p_total)
    
    
    for(i in 1:self$n_I){
        X[[i]] <- X_tmp[, match(self$var_names[[i]], all_var_name)]
        colnames(X[[i]]) <- self$var_names[[i]]
    }
    self$X <- X
    
    ins_prob <- matrix(0, nrow = self$n_B, ncol = self$n_I) # prob of an instance being positive
    z <- matrix(0, nrow = self$n_B, ncol = self$n_I) # instance label
    for(i in 1:self$n_I){
        p <- exp(cbind(1,X[[i]]) %*% beta[[i]])
        ins_prob[,i] <- p/(1+p)
        z[,i] <- rbinom(self$n_B, 1, ins_prob[,i])
    }
    self$ins_label <- z
    if(type == 1){
        self$bag_label <- as.numeric(rowSums(z) > 0)
    } else if(type == 2){
        self$bag_label <- as.numeric(rowSums(z) > 1)
    } 

})

HILoR$set("public", "predict", function(newX, thres = 0.5){
    stopifnot("Fit the model first!" = !is.null(self$coefficients))
    
    n_new <- nrow(newX[[1]])
    
    ins_prob <- matrix(0, nrow = n_new, ncol = self$n_I) # prob of an instance being positive
    for(i in 1:self$n_I){
        p <- exp(cbind(1,newX[[i]]) %*% self$coefficients[[i]])
        ins_prob[,i] <- p/(1+p)
    }
    ins_pred <- (ins_prob > thres) * 1
    
    if(self$type == 1){
        bag_prob <- apply(ins_prob, 1, function(x) 1-prod(1-x))
    } else if(self$type == 2){
        bag_prob <- apply(ins_prob, 1, function(x) 1-(1+sum(x/(1-x)))*prod(1-x))
    } 
    bag_pred <- as.integer(bag_prob > thres)

    return(list(ins_prob = ins_prob, ins_pred = ins_pred, 
                bag_prob = bag_prob, bag_pred = bag_pred))
})

HILoR$set("public", "inputLR", function(bag_label, X, var_names, type){
    stopifnot(is.null(dim(bag_label)), 
              "lenght of bag_label must equal to n_B" = length(bag_label) == self$n_B)
    stopifnot("X must be a matrix with dim n_X x unique(unlist(var_names))" = is.matrix(X), 
              nrow(X) == self$n_B,
              ncol(X) == unique(unlist(var_names)),
              !is.null(colnames(X)))
    stopifnot("type must be 1 or 2" = type %in% c(1,2))
    stopifnot(length(var_names) == n_I)
    
    self$var_names <- var_names
    self$bag_label <- bag_label
    self$type <- type
    
    all_var_names <- unique(unlist(var_names))
    
    X_HIL <- vector("list", length = self$n_I)
    for(i in 1:self$n_I){
        X_HIL[[i]] <- X[, match(self$var_names[[i]], all_var_names)]
        colnames(X_HIL[[i]]) <- self$var_names[[i]]
    }
    self$X <- X_HIL
    invisible(self)
})

source("~/github/HILoR/src/HILoR_fit.R")



