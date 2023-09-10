source("~/github/HILoR/src/HILoR.R")


n_B <- 50
n_I <- 4

beta <- list(beta1 = c(-4,-2,0, 0, 1),
             beta2 = c(-5,0,1,0,-2), 
             beta3 = c(-4,-3,1,0), 
             beta4 = c(-6,-2,0,0))

names(beta[[1]]) <- c("intercept", paste0("X", c(1,2,7,8)))
names(beta[[2]]) <- c("intercept", paste0("X", c(3,4,7,8)))
names(beta[[3]]) <- c("intercept", paste0("X", c(5,7,8)))
names(beta[[4]]) <- c("intercept", paste0("X", c(6,7,8)))

type <- 1

data <- HILoR$new(n_B, n_I)
data$generate(beta, type = type)



# split the data into training set and testing set
prop <- c(0.7, 0.3) # training, test

n_train <- round(n_B * prop[1], digits = 0)
n_test <- n_B - n_train
ind <- sample(n_B, n_B, replace = FALSE)

train <- data$clone()
test <- data$clone()
train$bag_label<- data$bag_label[ind[1:n_train]]
train$X <- lapply(data$X, function(x) x[ind[1:n_train],])
train$n_B <- n_train
test$bag_label<- data$bag_label[ind[(n_train+1):n_B]]
test$X <- lapply(data$X, function(x) x[ind[(n_train+1):n_B],])
test$n_B <- n_test

if(!is.null(data$ins_label)){
    train$ins_label<- data$ins_label[ind[1:n_train],]
    test$ins_label<- data$ins_label[ind[(n_train+1):n_B],]
}

# fit an HILoR model tow the training set
train$fit()

# prediction for the testing set using the fitted model
train$predict(test$X)


