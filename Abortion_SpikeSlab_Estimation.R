# Abortion, Estimation
# Data is created in Create_Abortion_DataSet.R

rm(list = ls())

# Set working directory
setwd("Insert_own_directory")

Packages <- c('tidyverse', 'plm')
invisible(lapply(Packages, require, character.only = TRUE))


# Load SpikeSlab with Treatment Effects
source('ss_treatment.R')

df <- readRDS(file = 'df_abortion.Rda')

# Controls

controls <- readRDS(file = 'controls.Rda')
tdums <- controls$tdums[3:13]; all_contviol <- controls$violence; all_contprop <- controls$property; all_contmurd <- controls$murder

df_viol <- df[, c('Dyviol','Dviol',all_contviol,tdums)]
df_prop <- df[, c('Dyprop','Dprop',all_contprop,tdums)]
df_murd <- df[, c('Dymurd','Dmurd',all_contmurd,tdums)]

# Variables included in original paper, p.404
var_viol <- c('Dyviol', 'Dviol', 'Lxxprison', 'Lxxpolice', 'Dxxunemp', 'Dxxincome',
              'Dxxpover', 'Dxxafdc15', 'Dxxgunlaw', 'Dxxbeer')

var_viol <- c(var_viol,tdums)

var_prop <- c('Dyprop', 'Dprop', 'Lxxprison', 'Lxxpolice', 'Dxxunemp', 'Dxxincome',
              'Dxxpover', 'Dxxafdc15', 'Dxxgunlaw', 'Dxxbeer')

var_prop <- c(var_prop,tdums)


var_murd <- c('Dymurd', 'Dmurd', 'Lxxprison', 'Lxxpolice', 'Dxxunemp', 'Dxxincome',
              'Dxxpover', 'Dxxafdc15', 'Dxxgunlaw', 'Dxxber')

var_murd <- c(var_murd,tdums)

#######################################################################################################
########################### With Eight Fixed Regressors from Original Paper ########################### 
#######################################################################################################

#### VIOLENCE ####
#### VIOLENCE ####
#### VIOLENCE ####

d <- df_viol$Dviol; X <- df_viol[,-(1:3)]; Y <- df_viol$Dyviol
X <- scale(X, scale = FALSE); d <- scale(d, scale = FALSE) # Center X and d

viol_index <- c()
for(i in var_viol) viol_index <- c(viol_index, which(colnames(X) == i))

fix_regressors <- rep(0,dim(X)[2])
for(i in viol_index) fix_regressors[i] <- 1

ptm <- proc.time()

draws <- ss_treatment(d,X,Y, fix_regr = fix_regressors, iter = 6000, print_iter = TRUE)

proc.time() - ptm

#### PROPERTY ####
#### PROPERTY ####
#### PROPERTY ####

d <- df_prop$Dprop; X <- df_prop[,-(1:3)]; Y <- df_prop$Dyprop
X <- scale(X, scale = FALSE); d <- scale(d, scale = FALSE) # Center X and d

prop_index <- c()
for(i in var_prop) prop_index <- c(prop_index, which(colnames(X) == i))

fix_regressors <- rep(0,dim(X)[2])
for(i in prop_index) fix_regressors[i] <- 1

ptm <- proc.time()

draws_prop <- ss_treatment(d,X,Y, fix_regr = fix_regressors, iter = 6000, print_iter = TRUE)

proc.time() - ptm

#### MURDER ####
#### MURDER ####
#### MURDER ####

d <- df_murd$Dmurd; X <- df_murd[,-(1:3)]; Y <- df_murd$Dymurd
X <- scale(X, scale = FALSE); d <- scale(d, scale = FALSE) # Center X and d

murd_index <- c()
for(i in var_murd) murd_index <- c(murd_index, which(colnames(X) == i))

fix_regressors <- rep(0,dim(X)[2])
for(i in murd_index) fix_regressors[i] <- 1

ptm <- proc.time()

draws_murd <- ss_treatment(d,X,Y, fix_regr = fix_regressors, iter = 6000, print_iter = TRUE)

proc.time() - ptm


#######################################################################################################
########################## No Fixed Regressors: Fully Data-Driven Selection ########################### 
#######################################################################################################

#### VIOLENCE ####
#### VIOLENCE ####
#### VIOLENCE ####

ptm <- proc.time()

d <- df_viol$Dviol; X <- df_viol[,-(1:3)]; Y <- df_viol$Dyviol
X <- scale(X, scale = FALSE); d <- scale(d, scale = FALSE) # Center X and d

draws_nofix <- ss_treatment(d,X,Y, iter = 6000, print_iter = TRUE)

proc.time() - ptm

#### PROPERTY ####
#### PROPERTY ####
#### PROPERTY ####

ptm <- proc.time()

d <- df_prop$Dprop; X <- df_prop[,-(1:3)]; Y <- df_prop$Dyprop
X <- scale(X, scale = FALSE); d <- scale(d, scale = FALSE) # Center X and d

draws_prop_nofix <- ss_treatment(d,X,Y, iter = 6000, print_iter = TRUE)

proc.time() - ptm

#### MURDER ####
#### MURDER ####
#### MURDER ####

ptm <- proc.time()

d <- df_murd$Dmurd; X <- df_murd[,-(1:3)]; Y <- df_murd$Dymurd
X <- scale(X, scale = FALSE); d <- scale(d, scale = FALSE) # Center X and d

draws_murd_nofix <- ss_treatment(d,X,Y, iter = 6000, print_iter = TRUE)

proc.time() - ptm



