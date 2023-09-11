library(ModelMetrics)
source("src/HILoR.R")


n_B <- 1000
n_I <- 4

true_beta <- list(beta1 = c(-4,-2,0, 0, 1),
             beta2 = c(-5,0,1,0,-2), 
             beta3 = c(-4,-3,1,0), 
             beta4 = c(-6,-2,0,0))

names(true_beta[[1]]) <- c("intercept", paste0("X", c(1,2,7,8)))
names(true_beta[[2]]) <- c("intercept", paste0("X", c(3,4,7,8)))
names(true_beta[[3]]) <- c("intercept", paste0("X", c(5,7,8)))
names(true_beta[[4]]) <- c("intercept", paste0("X", c(6,7,8)))

type <- 1

# split the data into training set and testing set
prop <- c(0.7, 0.3) # training, test

n_train <- round(n_B * prop[1], digits = 0)
n_test <- n_B - n_train

train <- HILoR$new(n_train, n_I)
train$generate(true_beta, type = type, seed = 2023)

test <- HILoR$new(n_test, n_I)
test$generate(true_beta, type = type, seed = 2024)


# fit an HILoR model tow the training set
train$fit()

# Summary of the fitted model
train$summary()
cat("\n")

# prediction for the testing set using the fitted model
pred <- train$predict(test$X)
bag_auc <- auc(test$bag_label, pred$bag_prob)
cat("Bag prediction AUC: ", round(bag_auc, 4), "\n")
ins_auc <- auc(as.vector(test$ins_label), as.vector(pred$ins_prob))
cat("Instance prediction AUC: ", round(ins_auc, 4), "\n")



