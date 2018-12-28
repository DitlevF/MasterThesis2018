# Abortion data in R
rm(list = ls())

setwd("Insert_own_directory")

# Read Publicy Available Data from Christian Hansen homepage
data <- read.delim('levitt_ex.dat')

# Load Packages
Packages <- c('tidyverse', 'fastDummies')
invisible(lapply(Packages, require, character.only = TRUE))

df <- as.tibble(data)

# Drop DC state
df <- subset(df, statenum != 9) # Removes 34 observations
dim(df)

# Drop years not used
df <- subset(df, year >= 85 & year <= 97 )
dim(df)

# Normalized trend variable
df$trend <- (df$year - 85)/12


#########################################################################################################
######################################### Generate Variables  ###########################################
#########################################################################################################

attach(df)

df$xxincome <- xxincome/100
df$xxpover <- xxpover/100
df$xxafdc15 <- xxafdc15/1000
df$xxbeer <-  xxbeer/100

# Generate year dummies
df$year <- as.factor(year) # class(df$year) # It's a factor now
df <- fastDummies::dummy_cols(df, select_columns = 'year')

tdums <- colnames(df)[19:31]
xx <- colnames(df)[10:17]


detach("package:plm", unload=TRUE) # If not detached, the below algorithm will use plm's diff, which doesn't work

# Generate differenced variables
xx_d <- c(NA)
for(i in xx){
  new_col <- paste0("D", i)
  
  df <- df %>%
    group_by(statenum) %>%
    mutate(z = get(i) - lag(get(i))) 
  
  names(df)[names(df) == 'z'] <- new_col # Change column name
  xx_d <- cbind(xx_d,new_col)
  
}
xx_d <- c(xx_d)[-1] # Differenced variables list

# Squared differenced variables

xx_d_sq <- c(NA)
for(i in xx_d){
  new_col <- paste0(i,"2")
  
  df <- df %>%
    group_by(statenum) %>%
    mutate(z = get(i)^2) 
  
  names(df)[names(df) == 'z'] <- new_col
  xx_d_sq <- cbind(xx_d_sq,new_col)
  
}

xx_d_sq <- c(xx_d_sq)[-1]

# Difference interactions

xx_d_int <- c(NA)
for(i in 1:(length(xx_d)-1)){
  for(j in (i+1):length(xx_d)){
    
    new_col <- paste0(xx_d[i],"X",xx_d[j])
    
    df <- df %>%
      group_by(statenum) %>%
      mutate(z = get(xx_d[i]) * get(xx_d[j])) 
    
    names(df)[names(df) == 'z'] <- new_col
    xx_d_int <- cbind(xx_d_int,new_col)
    
  }
}

xx_d_int <- c(xx_d_int)[-1]


# Lags
Lxx <- c(NA)

for(i in xx){
  new_col <- paste0('L',i)
  
  df <- df %>%
    group_by(statenum) %>%
    mutate(z = lag(get(i))) 
  
  names(df)[names(df) == 'z'] <- new_col
  Lxx <- cbind(Lxx,new_col)
  
}

Lxx <- c(Lxx)[-1]

# Squared Lags
Lxx2 <- c(NA)

for(i in Lxx){
  new_col <- paste0(i,"2")
  
  df <- df %>%
    group_by(statenum) %>%
    mutate(z = (get(i))^2) 
  
  names(df)[names(df) == 'z'] <- new_col
  Lxx2 <- cbind(Lxx2,new_col)
  
}

Lxx2 <- c(Lxx2)[-1]

# Means
Mxx <- c(NA)

for(i in xx){
  new_col <- paste0("M",i)
  df <- df %>%
    group_by(statenum) %>%
    mutate(z = mean(get(i)))
  names(df)[names(df) == 'z'] <- new_col
  
  Mxx <- cbind(Mxx,new_col)
  
}

Mxx <- c(Mxx)[-1]

# Squared Means
Mxx2 <- c(NA)

for(i in Mxx){
  new_col <- paste0(i,"2")
  
  df <- df %>%
    group_by(statenum) %>%
    mutate(z = (get(i))^2) 
  
  names(df)[names(df) == 'z'] <- new_col
  Mxx2 <- cbind(Mxx2,new_col)
  
}
Mxx2 <- c(Mxx2)[-1]

# Initial levels

xx0 <- c(NA)

for(i in xx){
  new_col <- paste0(i,"0")
  
  df <- df %>%
    group_by(statenum) %>%
    mutate(z = first(get(i))) 
  
  names(df)[names(df) == 'z'] <- new_col
  xx0 <- cbind(xx0,new_col)
  
}

xx0 <- c(xx0)[-1]

# Squared Initial Level
xx02 <- c(NA)

for(i in xx0){
  new_col <- paste0(i,"2")
  
  df <- df %>%
    group_by(statenum) %>%
    mutate(z = (get(i))^2) 
  
  names(df)[names(df) == 'z'] <- new_col
  xx02 <- cbind(xx02,new_col)
  
}
xx02 <- c(xx02)[-1]

# Initial Differences

Dxx0 <- c(NA)

for(i in xx_d){
  new_col <- paste0(i,"0")
  
  df <- df %>%
    group_by(statenum) %>%
    mutate(z = nth(get(i),2)) 
  
  names(df)[names(df) == 'z'] <- new_col
  Dxx0 <- cbind(Dxx0,new_col)
  
}
Dxx0 <- Dxx0[-1]

# Create Big List
biglist <- c(xx_d,xx_d_sq,xx_d_int,Lxx,Lxx2,Mxx,Mxx2,xx0,xx02,Dxx0)

# Interactions with trends

IntT1 <- c(NA)

for(i in biglist){
  new_col <- paste0(i,"Xt")
  
  df <- df %>%
    group_by(statenum) %>%
    mutate(z = (get(i))*trend) 
  
  names(df)[names(df) == 'z'] <- new_col
  IntT1 <- cbind(IntT1,new_col)
  
}
IntT1 <- c(IntT1)[-1]

# Interactions with trends squared

IntT2 <- c(NA)

for(i in biglist){
  new_col <- paste0(i,"Xt2")
  
  df <- df %>%
    group_by(statenum) %>%
    mutate(z = (get(i))*(trend^2)) 
  
  names(df)[names(df) == 'z'] <- new_col
  IntT2 <- cbind(IntT2,new_col)
  
}
IntT2 <- c(IntT2)[-1]

shared <- c(biglist, IntT1, IntT2)


#### Violence Specific Controls ####
i <- "efaviol"
j <- "Dviol"

df <- df %>%
  group_by(statenum) %>%
  mutate(Dviol = get(i) - lag(get(i))) 

df <- df %>%
  group_by(statenum) %>%
  mutate(viol0 = first(get(i))) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dviol0 = nth(get(j),2)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(viol02 = viol0^2) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dviol02 = Dviol0^2) 

df <- df %>%
  group_by(statenum) %>%
  mutate(viol0Xt = viol0*trend) 

df <- df %>%
  group_by(statenum) %>%
  mutate(viol0Xt2 = viol0*(trend^2)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(viol02Xt = viol02*(trend)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(viol02Xt2 = viol02*(trend^2)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dviol0Xt = Dviol0*(trend)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dviol0Xt2 = Dviol0*(trend^2)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dviol02Xt = Dviol02*(trend)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dviol02Xt2 = Dviol02*(trend^2)) 

contviol <- c('viol0', 'Dviol0', 'viol02', 'Dviol02', 'viol0Xt',
              'viol0Xt2', 'viol02Xt', 'viol02Xt2', 'Dviol0Xt', 'Dviol0Xt2', 'Dviol02Xt', 'Dviol02Xt2'  )

all_contviol <- c(shared, contviol)

#### Property Specific Controls ####
i <- "efaprop"
j <- "Dprop"

df <- df %>%
  group_by(statenum) %>%
  mutate(Dprop = get(i) - lag(get(i))) 

df <- df %>%
  group_by(statenum) %>%
  mutate(prop0 = first(get(i))) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dprop0 = nth(get(j),2)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(prop02 = prop0^2) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dprop02 = Dprop0^2) 

df <- df %>%
  group_by(statenum) %>%
  mutate(prop0Xt = prop0*trend) 

df <- df %>%
  group_by(statenum) %>%
  mutate(prop0Xt2 = prop0*(trend^2)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(prop02Xt = prop02*(trend)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(prop02Xt2 = prop02*(trend^2)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dprop0Xt = Dprop0*(trend)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dprop0Xt2 = Dprop0*(trend^2)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dprop02Xt = Dprop02*(trend)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dprop02Xt2 = Dprop02*(trend^2)) 

contprop <- c('prop0', 'Dprop0', 'prop02', 'Dprop02', 'prop0Xt',
              'prop0Xt2', 'prop02Xt', 'prop02Xt2', 'Dprop0Xt', 'Dprop0Xt2', 'Dprop02Xt', 'Dprop02Xt2')

all_contprop <- c(shared, contprop)

#### Murder Specific Controls ####
i <- "efamurd"
j <- "Dmurd"

df <- df %>%
  group_by(statenum) %>%
  mutate(Dmurd = get(i) - lag(get(i))) 

df <- df %>%
  group_by(statenum) %>%
  mutate(murd0 = first(get(i))) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dmurd0 = nth(get(j),2)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(murd02 = murd0^2) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dmurd02 = Dmurd0^2) 

df <- df %>%
  group_by(statenum) %>%
  mutate(murd0Xt = murd0*trend) 

df <- df %>%
  group_by(statenum) %>%
  mutate(murd0Xt2 = murd0*(trend^2)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(murd02Xt = murd02*(trend)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(murd02Xt2 = murd02*(trend^2)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dmurd0Xt = Dmurd0*(trend)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dmurd0Xt2 = Dmurd0*(trend^2)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dmurd02Xt = Dmurd02*(trend)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dmurd02Xt2 = Dmurd02*(trend^2)) 

contmurd <- c('murd0', 'Dmurd0', 'murd02', 'Dmurd02', 'murd0Xt',
              'murd0Xt2', 'murd02Xt', 'murd02Xt2', 'Dmurd0Xt', 'Dmurd0Xt2', 'Dmurd02Xt', 'Dmurd02Xt2')

all_contmurd <- c(shared, contmurd)

#########################################################################################################
################################################ ESTIMATION #############################################
#########################################################################################################

# Differenced outcomes
df <- df %>%
  group_by(statenum) %>%
  mutate(Dyviol = lpc_viol - lag(lpc_viol)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dyprop = lpc_prop - lag(lpc_prop)) 

df <- df %>%
  group_by(statenum) %>%
  mutate(Dymurd = lpc_murd - lag(lpc_murd)) 

# Remove first year
df <- subset(df, trend != 0)

df_viol <- df[, c('Dyviol','Dviol',tdums,all_contviol)]
df_prop <- df[, c('Dyprop','Dprop',tdums,all_contprop)]
df_murd <- df[, c('Dmurd','Dmurd',tdums,all_contmurd)]

controls <- list('violence' = all_contviol, 'property' = all_contprop, 'murder' = all_contmurd, 'tdums' = tdums)

saveRDS(controls, file = 'controls.Rda')
saveRDS(df, file = 'df_abortion.Rda')


