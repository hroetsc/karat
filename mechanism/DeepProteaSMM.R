### karat projetc - PCPS mechanism ###
# description:  try inference of STS weights with single layer neural net
# input:        STS + DATA matrix
# output:       weights
# author:       HPR

library(dplyr)
library(stringr)
library(tensorflow)
library(tfdatasets)
library(keras)


### INPUT ###
load("data/ProteaSMM/PSP_STS/DATA.RData")


### MAIN PART ###
# ----- train/test split -----
pseudo = 1e-05

# get data matrices
X = DATA$X
t = log(DATA$t/100+pseudo)

ktest = which(DATA$substrateIDs %in% c("MM582","MM537"))

Xtest = X[ktest, ]
ttest = t[ktest, ]

Xtrain = X[-ktest, ]
ttrain = t[-ktest, ]

numParam = ncol(X)

# ----- construct simple model -----
build_model = function() {
  input <- layer_input(shape = numParam)
  
  output <- input %>%
    # layer_dense(units = numParam, activation = "linear", kernel_initializer = "he_uniform") %>%
    layer_dense(units = 1, activation = "linear", use_bias = F,
                kernel_initializer = initializer_random_uniform(minval = 0, maxval = 2, seed = 42),
                kernel_constraint = constraint_minmaxnorm(min_value = 0, max_value = 2))
  
  model <- keras_model(input, output)
  
  model %>% 
    compile(
      loss = "mse",
      optimizer = optimizer_adam(learning_rate = 1e-04),
      metrics = list("mean_squared_error")
    )
  
  return(model)
}

model <- build_model()

# ----- training -----
early_stop <- callback_early_stopping(monitor = "val_loss", patience = 10)

history <- model %>%
  fit(x = Xtrain,
      y = rowMeans(ttrain, na.rm = T),
      epochs = 2e02, validation_split = .2, verbose = 2,
      callbacks = list(early_stop))

parameters = get_weights(model)
parameters = parameters[[1]]

# ----- prediction -----
tsim <- model %>% predict(Xtest)

pcc = cor(rowMeans(ttest, na.rm = T), tsim)
plot(x = rowMeans(ttest, na.rm = T),
     y = tsim, pch = 16,
     xlab = "true", ylab = "predicted",
     main = paste0("PCC = ", round(pcc,4)))
abline(a = 0, b = 1, col = "red")





