##########################################
#      Convergent Cross Mapping       
#
# Created by: Sarah Pirikahu
# Creation date: 24 March 2023
##########################################



# Embedding dimension (i.e. the n umber of lags used to build up the shadow manifold)
# determined based on the prediction skill of the model.The rEDM vingette 
# https://ha0ye.github.io/rEDM/articles/rEDM.html uses simplex function. 
# EmbedDimension is a wrapper around the simplex function to get E only out 

# specify library set = how much data to fit too
lib_max <- dim(d1)[1]/2
lib <- paste0("1 ", lib_max)

# specify pred = which data to predict on
pred <- paste0(lib_max-1, " ", dim(d1)[1])

# Get E for H2
# ....I do wonder if this is enough data to accurately determine E? 
E_h2 <- EmbedDimension(dataFrame = d1, columns = "H2", target = "H2",
                        lib = lib, pred = pred, showPlot = TRUE)
E_h2 <- E_h2 %>% slice_max(rho) # keep the row with max prediction skill
E_h2 <- E_h2[,1] # keep just the value of E

# Get E for H1 
E_h1 <- EmbedDimension(dataFrame = d1, columns = "H1", target = "H1",
                       lib = lib, pred = pred, showPlot = TRUE)
E_h1 <- E_h1 %>% slice_max(rho) 
E_h1 <- E_h1[,1] 




