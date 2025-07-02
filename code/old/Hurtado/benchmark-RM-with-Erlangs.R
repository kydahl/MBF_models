####################################################################################
## Benchmarking the Rosenzweig-MacArthur model numerical solutions
## ----------------------------------------------------------------
##  Case 1: Standard: No maturation delay, exponentially distributed predator lifetime
##  Case 2: Classic LCT (Erlang maturation time and adult life stage duration), hard-coded
##  Case 3: Phase-type/GLCT implementation of Case 2
####################################################################################
#
#   NOTE: This R script was run on a 64-bit Windows 10 machine using R v3.6.3.
#         Running this on other operating systems may require small modifications
#         like replacing the win.graph() calls with x11() or quartz() calls.
#

library(deSolve)
library(rbenchmark)
library(ggplot2)

## Parameterization....
reps<-140
IC = c(N=1000,Y=10); 
params = c(r = 1, K = 1000, a = 5, h = 500, chi = 0.5, mux = 0.5, muy = 1);

####################################################################################
## Function definitions

# Standard Rosenzweig-MacArthur model
RM.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];

  N=z[1] # prey population
  Y=z[2] # predator population
  
  dN = r*N*(1-N/K)-(a/(h+N))*Y*N
  dY = chi*(a/(h+N))*N*Y-Y/muy
  return(list(c(dN,dY)))
}


# Rosenzweig-MacArthur with Erlang maturation and adult-life-stage times, in a GLCT framework
RMpt.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[6];
  #kx = ps[['kx']] # Number of substates in E, and...
  #ky = ps[['ky']] # ... I.
  
  N = z[1] # prey population
  X = z[1+(1:kx)] # immature predator population substates
  Y = z[1+kx+(1:ky)] # mature predator population sub-states
  P = onesyt %*% Y # sum Y_i

  dN = r*N*(1-N/K)- a/(h+N)*P*N
  dX = chi*(a/(h+N))*P[1,1]*N * ax + Axt %*% X
  dY = ay %*% -onesxt %*% Axt %*% X + Ayt %*% Y
  
  return(list(c(dN,as.numeric(dX),as.numeric(dY))))
}

RMpt.init <- function(ps) {
  # Unpack some parameter values...
  mux=ps[["mux"]] # mean maturation time
  muy=ps[["muy"]] # mean time spent in mature life stage
  kx <<- ps[['kx']] # Number of substates in E, and...
  ky <<- ps[['ky']] # ... I.
  
  # These are Erlang distributions framed in a Phase-type distribution context,
  # where vectors ax = (1 0 ... 0) and matrices Ax are as follows...
  
  ax = matrix(0, nrow=kx, ncol=1); ax[1] = 1;
  Ax = kx/mux*(diag(rep(-1,kx),kx)); if(kx>1) for(i in 1:(kx-1)) {Ax[i,i+1] = kx/mux}
  
  ay = matrix(0, nrow=kx, ncol=1); ay[1] = 1;
  Ay = ky/muy*(diag(rep(-1,ky),ky)); if(ky>1) for(i in 1:(ky-1)) {Ay[i,i+1] = ky/muy}

  # Initial conditions
  z0=numeric(1+kx+ky) # initialize the state variable vector with 0s
  z0[1] <- IC[["N"]]  # 1 in the initial exposed class 
  z0[1+kx+1] <- IC[["Y"]]    # susceptibles = (PopSize - 1)/PopSize
  
  # Set some global variables that the RMpt.ode function can access...
  ax  <<- ax
  Axt <<- t(Ax)
  ay  <<- ay
  Ayt <<- t(Ay)
  onesxt <<- matrix(1,ncol=kx,nrow=1)
  onesx <<- matrix(1,nrow=kx,ncol=1)
  onesyt <<- matrix(1,ncol=ky,nrow=1)
  onesy <<- matrix(1,nrow=ky,ncol=1)
  ICs <<- z0
}

# R.-M. with hard-coded Erlang dwell times
RMlct1.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dY1 = kx/mux*z[2] - z[3]*ky/muy

  return(list(c(dN, dX1, dY1)))
}

RMlct2.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dY1 = kx/mux*z[3] - z[4]*ky/muy
  dY2 = ky/muy*z[4] - z[5]*ky/muy
  
  return(list(c(dN, dX1, dX2, dY1, dY2)))
}

RMlct3.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dY1 = kx/mux*z[4] - z[5]*ky/muy
  dY2 = ky/muy*z[5] - z[6]*ky/muy
  dY3 = ky/muy*z[6] - z[7]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dY1, dY2, dY3)))
}

RMlct4.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dY1 = kx/mux*z[5] - z[6]*ky/muy
  dY2 = ky/muy*z[6] - z[7]*ky/muy
  dY3 = ky/muy*z[7] - z[8]*ky/muy
  dY4 = ky/muy*z[8] - z[9]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dY1, dY2, dY3, dY4)))
}

RMlct5.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dY1 = kx/mux*z[6] - z[7]*ky/muy
  dY2 = ky/muy*z[7] - z[8]*ky/muy
  dY3 = ky/muy*z[8] - z[9]*ky/muy
  dY4 = ky/muy*z[9] - z[10]*ky/muy
  dY5 = ky/muy*z[10] - z[11]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dY1, dY2, dY3, dY4, dY5)))
}

RMlct6.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dY1 = kx/mux*z[7] - z[8]*ky/muy
  dY2 = ky/muy*z[8] - z[9]*ky/muy
  dY3 = ky/muy*z[9] - z[10]*ky/muy
  dY4 = ky/muy*z[10] - z[11]*ky/muy
  dY5 = ky/muy*z[11] - z[12]*ky/muy
  dY6 = ky/muy*z[12] - z[13]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dY1, dY2, dY3, dY4, dY5, dY6)))
}

RMlct7.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dY1 = kx/mux*z[8] - z[9]*ky/muy
  dY2 = ky/muy*z[9] - z[10]*ky/muy
  dY3 = ky/muy*z[10] - z[11]*ky/muy
  dY4 = ky/muy*z[11] - z[12]*ky/muy
  dY5 = ky/muy*z[12] - z[13]*ky/muy
  dY6 = ky/muy*z[13] - z[14]*ky/muy
  dY7 = ky/muy*z[14] - z[15]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dY1, dY2, dY3, dY4, dY5, dY6, dY7)))
}

RMlct8.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dY1 = kx/mux*z[9] - z[10]*ky/muy
  dY2 = ky/muy*z[10] - z[11]*ky/muy
  dY3 = ky/muy*z[11] - z[12]*ky/muy
  dY4 = ky/muy*z[12] - z[13]*ky/muy
  dY5 = ky/muy*z[13] - z[14]*ky/muy
  dY6 = ky/muy*z[14] - z[15]*ky/muy
  dY7 = ky/muy*z[15] - z[16]*ky/muy
  dY8 = ky/muy*z[16] - z[17]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8)))
}

RMlct9.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dY1 = kx/mux*z[10] - z[11]*ky/muy
  dY2 = ky/muy*z[11] - z[12]*ky/muy
  dY3 = ky/muy*z[12] - z[13]*ky/muy
  dY4 = ky/muy*z[13] - z[14]*ky/muy
  dY5 = ky/muy*z[14] - z[15]*ky/muy
  dY6 = ky/muy*z[15] - z[16]*ky/muy
  dY7 = ky/muy*z[16] - z[17]*ky/muy
  dY8 = ky/muy*z[17] - z[18]*ky/muy
  dY9 = ky/muy*z[18] - z[19]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9)))
}

RMlct10.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dY1 = kx/mux*z[11] - z[12]*ky/muy
  dY2 = ky/muy*z[12] - z[13]*ky/muy
  dY3 = ky/muy*z[13] - z[14]*ky/muy
  dY4 = ky/muy*z[14] - z[15]*ky/muy
  dY5 = ky/muy*z[15] - z[16]*ky/muy
  dY6 = ky/muy*z[16] - z[17]*ky/muy
  dY7 = ky/muy*z[17] - z[18]*ky/muy
  dY8 = ky/muy*z[18] - z[19]*ky/muy
  dY9 = ky/muy*z[19] - z[20]*ky/muy
  dY10 = ky/muy*z[20] - z[21]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10)))
}

RMlct11.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dY1 = kx/mux*z[12] - z[13]*ky/muy
  dY2 = ky/muy*z[13] - z[14]*ky/muy
  dY3 = ky/muy*z[14] - z[15]*ky/muy
  dY4 = ky/muy*z[15] - z[16]*ky/muy
  dY5 = ky/muy*z[16] - z[17]*ky/muy
  dY6 = ky/muy*z[17] - z[18]*ky/muy
  dY7 = ky/muy*z[18] - z[19]*ky/muy
  dY8 = ky/muy*z[19] - z[20]*ky/muy
  dY9 = ky/muy*z[20] - z[21]*ky/muy
  dY10 = ky/muy*z[21] - z[22]*ky/muy
  dY11 = ky/muy*z[22] - z[23]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11)))
}

RMlct12.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dY1 = kx/mux*z[13] - z[14]*ky/muy
  dY2 = ky/muy*z[14] - z[15]*ky/muy
  dY3 = ky/muy*z[15] - z[16]*ky/muy
  dY4 = ky/muy*z[16] - z[17]*ky/muy
  dY5 = ky/muy*z[17] - z[18]*ky/muy
  dY6 = ky/muy*z[18] - z[19]*ky/muy
  dY7 = ky/muy*z[19] - z[20]*ky/muy
  dY8 = ky/muy*z[20] - z[21]*ky/muy
  dY9 = ky/muy*z[21] - z[22]*ky/muy
  dY10 = ky/muy*z[22] - z[23]*ky/muy
  dY11 = ky/muy*z[23] - z[24]*ky/muy
  dY12 = ky/muy*z[24] - z[25]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12)))
}

RMlct13.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dY1 = kx/mux*z[14] - z[15]*ky/muy
  dY2 = ky/muy*z[15] - z[16]*ky/muy
  dY3 = ky/muy*z[16] - z[17]*ky/muy
  dY4 = ky/muy*z[17] - z[18]*ky/muy
  dY5 = ky/muy*z[18] - z[19]*ky/muy
  dY6 = ky/muy*z[19] - z[20]*ky/muy
  dY7 = ky/muy*z[20] - z[21]*ky/muy
  dY8 = ky/muy*z[21] - z[22]*ky/muy
  dY9 = ky/muy*z[22] - z[23]*ky/muy
  dY10 = ky/muy*z[23] - z[24]*ky/muy
  dY11 = ky/muy*z[24] - z[25]*ky/muy
  dY12 = ky/muy*z[25] - z[26]*ky/muy
  dY13 = ky/muy*z[26] - z[27]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13)))
}

RMlct14.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dX14 = kx/mux*z[14] - z[15]*kx/mux
  dY1 = kx/mux*z[15] - z[16]*ky/muy
  dY2 = ky/muy*z[16] - z[17]*ky/muy
  dY3 = ky/muy*z[17] - z[18]*ky/muy
  dY4 = ky/muy*z[18] - z[19]*ky/muy
  dY5 = ky/muy*z[19] - z[20]*ky/muy
  dY6 = ky/muy*z[20] - z[21]*ky/muy
  dY7 = ky/muy*z[21] - z[22]*ky/muy
  dY8 = ky/muy*z[22] - z[23]*ky/muy
  dY9 = ky/muy*z[23] - z[24]*ky/muy
  dY10 = ky/muy*z[24] - z[25]*ky/muy
  dY11 = ky/muy*z[25] - z[26]*ky/muy
  dY12 = ky/muy*z[26] - z[27]*ky/muy
  dY13 = ky/muy*z[27] - z[28]*ky/muy
  dY14 = ky/muy*z[28] - z[29]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dX14, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13, dY14)))
}

RMlct15.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dX14 = kx/mux*z[14] - z[15]*kx/mux
  dX15 = kx/mux*z[15] - z[16]*kx/mux
  dY1 = kx/mux*z[16] - z[17]*ky/muy
  dY2 = ky/muy*z[17] - z[18]*ky/muy
  dY3 = ky/muy*z[18] - z[19]*ky/muy
  dY4 = ky/muy*z[19] - z[20]*ky/muy
  dY5 = ky/muy*z[20] - z[21]*ky/muy
  dY6 = ky/muy*z[21] - z[22]*ky/muy
  dY7 = ky/muy*z[22] - z[23]*ky/muy
  dY8 = ky/muy*z[23] - z[24]*ky/muy
  dY9 = ky/muy*z[24] - z[25]*ky/muy
  dY10 = ky/muy*z[25] - z[26]*ky/muy
  dY11 = ky/muy*z[26] - z[27]*ky/muy
  dY12 = ky/muy*z[27] - z[28]*ky/muy
  dY13 = ky/muy*z[28] - z[29]*ky/muy
  dY14 = ky/muy*z[29] - z[30]*ky/muy
  dY15 = ky/muy*z[30] - z[31]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dX14, dX15, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13, dY14, dY15)))
}

RMlct16.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dX14 = kx/mux*z[14] - z[15]*kx/mux
  dX15 = kx/mux*z[15] - z[16]*kx/mux
  dX16 = kx/mux*z[16] - z[17]*kx/mux
  dY1 = kx/mux*z[17] - z[18]*ky/muy
  dY2 = ky/muy*z[18] - z[19]*ky/muy
  dY3 = ky/muy*z[19] - z[20]*ky/muy
  dY4 = ky/muy*z[20] - z[21]*ky/muy
  dY5 = ky/muy*z[21] - z[22]*ky/muy
  dY6 = ky/muy*z[22] - z[23]*ky/muy
  dY7 = ky/muy*z[23] - z[24]*ky/muy
  dY8 = ky/muy*z[24] - z[25]*ky/muy
  dY9 = ky/muy*z[25] - z[26]*ky/muy
  dY10 = ky/muy*z[26] - z[27]*ky/muy
  dY11 = ky/muy*z[27] - z[28]*ky/muy
  dY12 = ky/muy*z[28] - z[29]*ky/muy
  dY13 = ky/muy*z[29] - z[30]*ky/muy
  dY14 = ky/muy*z[30] - z[31]*ky/muy
  dY15 = ky/muy*z[31] - z[32]*ky/muy
  dY16 = ky/muy*z[32] - z[33]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dX14, dX15, dX16, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13, dY14, dY15, dY16)))
}

RMlct17.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dX14 = kx/mux*z[14] - z[15]*kx/mux
  dX15 = kx/mux*z[15] - z[16]*kx/mux
  dX16 = kx/mux*z[16] - z[17]*kx/mux
  dX17 = kx/mux*z[17] - z[18]*kx/mux
  dY1 = kx/mux*z[18] - z[19]*ky/muy
  dY2 = ky/muy*z[19] - z[20]*ky/muy
  dY3 = ky/muy*z[20] - z[21]*ky/muy
  dY4 = ky/muy*z[21] - z[22]*ky/muy
  dY5 = ky/muy*z[22] - z[23]*ky/muy
  dY6 = ky/muy*z[23] - z[24]*ky/muy
  dY7 = ky/muy*z[24] - z[25]*ky/muy
  dY8 = ky/muy*z[25] - z[26]*ky/muy
  dY9 = ky/muy*z[26] - z[27]*ky/muy
  dY10 = ky/muy*z[27] - z[28]*ky/muy
  dY11 = ky/muy*z[28] - z[29]*ky/muy
  dY12 = ky/muy*z[29] - z[30]*ky/muy
  dY13 = ky/muy*z[30] - z[31]*ky/muy
  dY14 = ky/muy*z[31] - z[32]*ky/muy
  dY15 = ky/muy*z[32] - z[33]*ky/muy
  dY16 = ky/muy*z[33] - z[34]*ky/muy
  dY17 = ky/muy*z[34] - z[35]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dX14, dX15, dX16, dX17, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13, dY14, dY15, dY16, dY17)))
}

RMlct18.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dX14 = kx/mux*z[14] - z[15]*kx/mux
  dX15 = kx/mux*z[15] - z[16]*kx/mux
  dX16 = kx/mux*z[16] - z[17]*kx/mux
  dX17 = kx/mux*z[17] - z[18]*kx/mux
  dX18 = kx/mux*z[18] - z[19]*kx/mux
  dY1 = kx/mux*z[19] - z[20]*ky/muy
  dY2 = ky/muy*z[20] - z[21]*ky/muy
  dY3 = ky/muy*z[21] - z[22]*ky/muy
  dY4 = ky/muy*z[22] - z[23]*ky/muy
  dY5 = ky/muy*z[23] - z[24]*ky/muy
  dY6 = ky/muy*z[24] - z[25]*ky/muy
  dY7 = ky/muy*z[25] - z[26]*ky/muy
  dY8 = ky/muy*z[26] - z[27]*ky/muy
  dY9 = ky/muy*z[27] - z[28]*ky/muy
  dY10 = ky/muy*z[28] - z[29]*ky/muy
  dY11 = ky/muy*z[29] - z[30]*ky/muy
  dY12 = ky/muy*z[30] - z[31]*ky/muy
  dY13 = ky/muy*z[31] - z[32]*ky/muy
  dY14 = ky/muy*z[32] - z[33]*ky/muy
  dY15 = ky/muy*z[33] - z[34]*ky/muy
  dY16 = ky/muy*z[34] - z[35]*ky/muy
  dY17 = ky/muy*z[35] - z[36]*ky/muy
  dY18 = ky/muy*z[36] - z[37]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dX14, dX15, dX16, dX17, dX18, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13, dY14, dY15, dY16, dY17, dY18)))
}

RMlct19.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dX14 = kx/mux*z[14] - z[15]*kx/mux
  dX15 = kx/mux*z[15] - z[16]*kx/mux
  dX16 = kx/mux*z[16] - z[17]*kx/mux
  dX17 = kx/mux*z[17] - z[18]*kx/mux
  dX18 = kx/mux*z[18] - z[19]*kx/mux
  dX19 = kx/mux*z[19] - z[20]*kx/mux
  dY1 = kx/mux*z[20] - z[21]*ky/muy
  dY2 = ky/muy*z[21] - z[22]*ky/muy
  dY3 = ky/muy*z[22] - z[23]*ky/muy
  dY4 = ky/muy*z[23] - z[24]*ky/muy
  dY5 = ky/muy*z[24] - z[25]*ky/muy
  dY6 = ky/muy*z[25] - z[26]*ky/muy
  dY7 = ky/muy*z[26] - z[27]*ky/muy
  dY8 = ky/muy*z[27] - z[28]*ky/muy
  dY9 = ky/muy*z[28] - z[29]*ky/muy
  dY10 = ky/muy*z[29] - z[30]*ky/muy
  dY11 = ky/muy*z[30] - z[31]*ky/muy
  dY12 = ky/muy*z[31] - z[32]*ky/muy
  dY13 = ky/muy*z[32] - z[33]*ky/muy
  dY14 = ky/muy*z[33] - z[34]*ky/muy
  dY15 = ky/muy*z[34] - z[35]*ky/muy
  dY16 = ky/muy*z[35] - z[36]*ky/muy
  dY17 = ky/muy*z[36] - z[37]*ky/muy
  dY18 = ky/muy*z[37] - z[38]*ky/muy
  dY19 = ky/muy*z[38] - z[39]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dX14, dX15, dX16, dX17, dX18, dX19, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13, dY14, dY15, dY16, dY17, dY18, dY19)))
}

RMlct20.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dX14 = kx/mux*z[14] - z[15]*kx/mux
  dX15 = kx/mux*z[15] - z[16]*kx/mux
  dX16 = kx/mux*z[16] - z[17]*kx/mux
  dX17 = kx/mux*z[17] - z[18]*kx/mux
  dX18 = kx/mux*z[18] - z[19]*kx/mux
  dX19 = kx/mux*z[19] - z[20]*kx/mux
  dX20 = kx/mux*z[20] - z[21]*kx/mux
  dY1 = kx/mux*z[21] - z[22]*ky/muy
  dY2 = ky/muy*z[22] - z[23]*ky/muy
  dY3 = ky/muy*z[23] - z[24]*ky/muy
  dY4 = ky/muy*z[24] - z[25]*ky/muy
  dY5 = ky/muy*z[25] - z[26]*ky/muy
  dY6 = ky/muy*z[26] - z[27]*ky/muy
  dY7 = ky/muy*z[27] - z[28]*ky/muy
  dY8 = ky/muy*z[28] - z[29]*ky/muy
  dY9 = ky/muy*z[29] - z[30]*ky/muy
  dY10 = ky/muy*z[30] - z[31]*ky/muy
  dY11 = ky/muy*z[31] - z[32]*ky/muy
  dY12 = ky/muy*z[32] - z[33]*ky/muy
  dY13 = ky/muy*z[33] - z[34]*ky/muy
  dY14 = ky/muy*z[34] - z[35]*ky/muy
  dY15 = ky/muy*z[35] - z[36]*ky/muy
  dY16 = ky/muy*z[36] - z[37]*ky/muy
  dY17 = ky/muy*z[37] - z[38]*ky/muy
  dY18 = ky/muy*z[38] - z[39]*ky/muy
  dY19 = ky/muy*z[39] - z[40]*ky/muy
  dY20 = ky/muy*z[40] - z[41]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dX14, dX15, dX16, dX17, dX18, dX19, dX20, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13, dY14, dY15, dY16, dY17, dY18, dY19, dY20)))
}

RMlct21.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dX14 = kx/mux*z[14] - z[15]*kx/mux
  dX15 = kx/mux*z[15] - z[16]*kx/mux
  dX16 = kx/mux*z[16] - z[17]*kx/mux
  dX17 = kx/mux*z[17] - z[18]*kx/mux
  dX18 = kx/mux*z[18] - z[19]*kx/mux
  dX19 = kx/mux*z[19] - z[20]*kx/mux
  dX20 = kx/mux*z[20] - z[21]*kx/mux
  dX21 = kx/mux*z[21] - z[22]*kx/mux
  dY1 = kx/mux*z[22] - z[23]*ky/muy
  dY2 = ky/muy*z[23] - z[24]*ky/muy
  dY3 = ky/muy*z[24] - z[25]*ky/muy
  dY4 = ky/muy*z[25] - z[26]*ky/muy
  dY5 = ky/muy*z[26] - z[27]*ky/muy
  dY6 = ky/muy*z[27] - z[28]*ky/muy
  dY7 = ky/muy*z[28] - z[29]*ky/muy
  dY8 = ky/muy*z[29] - z[30]*ky/muy
  dY9 = ky/muy*z[30] - z[31]*ky/muy
  dY10 = ky/muy*z[31] - z[32]*ky/muy
  dY11 = ky/muy*z[32] - z[33]*ky/muy
  dY12 = ky/muy*z[33] - z[34]*ky/muy
  dY13 = ky/muy*z[34] - z[35]*ky/muy
  dY14 = ky/muy*z[35] - z[36]*ky/muy
  dY15 = ky/muy*z[36] - z[37]*ky/muy
  dY16 = ky/muy*z[37] - z[38]*ky/muy
  dY17 = ky/muy*z[38] - z[39]*ky/muy
  dY18 = ky/muy*z[39] - z[40]*ky/muy
  dY19 = ky/muy*z[40] - z[41]*ky/muy
  dY20 = ky/muy*z[41] - z[42]*ky/muy
  dY21 = ky/muy*z[42] - z[43]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dX14, dX15, dX16, dX17, dX18, dX19, dX20, dX21, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13, dY14, dY15, dY16, dY17, dY18, dY19, dY20, dY21)))
}

RMlct22.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dX14 = kx/mux*z[14] - z[15]*kx/mux
  dX15 = kx/mux*z[15] - z[16]*kx/mux
  dX16 = kx/mux*z[16] - z[17]*kx/mux
  dX17 = kx/mux*z[17] - z[18]*kx/mux
  dX18 = kx/mux*z[18] - z[19]*kx/mux
  dX19 = kx/mux*z[19] - z[20]*kx/mux
  dX20 = kx/mux*z[20] - z[21]*kx/mux
  dX21 = kx/mux*z[21] - z[22]*kx/mux
  dX22 = kx/mux*z[22] - z[23]*kx/mux
  dY1 = kx/mux*z[23] - z[24]*ky/muy
  dY2 = ky/muy*z[24] - z[25]*ky/muy
  dY3 = ky/muy*z[25] - z[26]*ky/muy
  dY4 = ky/muy*z[26] - z[27]*ky/muy
  dY5 = ky/muy*z[27] - z[28]*ky/muy
  dY6 = ky/muy*z[28] - z[29]*ky/muy
  dY7 = ky/muy*z[29] - z[30]*ky/muy
  dY8 = ky/muy*z[30] - z[31]*ky/muy
  dY9 = ky/muy*z[31] - z[32]*ky/muy
  dY10 = ky/muy*z[32] - z[33]*ky/muy
  dY11 = ky/muy*z[33] - z[34]*ky/muy
  dY12 = ky/muy*z[34] - z[35]*ky/muy
  dY13 = ky/muy*z[35] - z[36]*ky/muy
  dY14 = ky/muy*z[36] - z[37]*ky/muy
  dY15 = ky/muy*z[37] - z[38]*ky/muy
  dY16 = ky/muy*z[38] - z[39]*ky/muy
  dY17 = ky/muy*z[39] - z[40]*ky/muy
  dY18 = ky/muy*z[40] - z[41]*ky/muy
  dY19 = ky/muy*z[41] - z[42]*ky/muy
  dY20 = ky/muy*z[42] - z[43]*ky/muy
  dY21 = ky/muy*z[43] - z[44]*ky/muy
  dY22 = ky/muy*z[44] - z[45]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dX14, dX15, dX16, dX17, dX18, dX19, dX20, dX21, dX22, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13, dY14, dY15, dY16, dY17, dY18, dY19, dY20, dY21, dY22)))
}

RMlct23.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dX14 = kx/mux*z[14] - z[15]*kx/mux
  dX15 = kx/mux*z[15] - z[16]*kx/mux
  dX16 = kx/mux*z[16] - z[17]*kx/mux
  dX17 = kx/mux*z[17] - z[18]*kx/mux
  dX18 = kx/mux*z[18] - z[19]*kx/mux
  dX19 = kx/mux*z[19] - z[20]*kx/mux
  dX20 = kx/mux*z[20] - z[21]*kx/mux
  dX21 = kx/mux*z[21] - z[22]*kx/mux
  dX22 = kx/mux*z[22] - z[23]*kx/mux
  dX23 = kx/mux*z[23] - z[24]*kx/mux
  dY1 = kx/mux*z[24] - z[25]*ky/muy
  dY2 = ky/muy*z[25] - z[26]*ky/muy
  dY3 = ky/muy*z[26] - z[27]*ky/muy
  dY4 = ky/muy*z[27] - z[28]*ky/muy
  dY5 = ky/muy*z[28] - z[29]*ky/muy
  dY6 = ky/muy*z[29] - z[30]*ky/muy
  dY7 = ky/muy*z[30] - z[31]*ky/muy
  dY8 = ky/muy*z[31] - z[32]*ky/muy
  dY9 = ky/muy*z[32] - z[33]*ky/muy
  dY10 = ky/muy*z[33] - z[34]*ky/muy
  dY11 = ky/muy*z[34] - z[35]*ky/muy
  dY12 = ky/muy*z[35] - z[36]*ky/muy
  dY13 = ky/muy*z[36] - z[37]*ky/muy
  dY14 = ky/muy*z[37] - z[38]*ky/muy
  dY15 = ky/muy*z[38] - z[39]*ky/muy
  dY16 = ky/muy*z[39] - z[40]*ky/muy
  dY17 = ky/muy*z[40] - z[41]*ky/muy
  dY18 = ky/muy*z[41] - z[42]*ky/muy
  dY19 = ky/muy*z[42] - z[43]*ky/muy
  dY20 = ky/muy*z[43] - z[44]*ky/muy
  dY21 = ky/muy*z[44] - z[45]*ky/muy
  dY22 = ky/muy*z[45] - z[46]*ky/muy
  dY23 = ky/muy*z[46] - z[47]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dX14, dX15, dX16, dX17, dX18, dX19, dX20, dX21, dX22, dX23, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13, dY14, dY15, dY16, dY17, dY18, dY19, dY20, dY21, dY22, dY23)))
}

RMlct24.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dX14 = kx/mux*z[14] - z[15]*kx/mux
  dX15 = kx/mux*z[15] - z[16]*kx/mux
  dX16 = kx/mux*z[16] - z[17]*kx/mux
  dX17 = kx/mux*z[17] - z[18]*kx/mux
  dX18 = kx/mux*z[18] - z[19]*kx/mux
  dX19 = kx/mux*z[19] - z[20]*kx/mux
  dX20 = kx/mux*z[20] - z[21]*kx/mux
  dX21 = kx/mux*z[21] - z[22]*kx/mux
  dX22 = kx/mux*z[22] - z[23]*kx/mux
  dX23 = kx/mux*z[23] - z[24]*kx/mux
  dX24 = kx/mux*z[24] - z[25]*kx/mux
  dY1 = kx/mux*z[25] - z[26]*ky/muy
  dY2 = ky/muy*z[26] - z[27]*ky/muy
  dY3 = ky/muy*z[27] - z[28]*ky/muy
  dY4 = ky/muy*z[28] - z[29]*ky/muy
  dY5 = ky/muy*z[29] - z[30]*ky/muy
  dY6 = ky/muy*z[30] - z[31]*ky/muy
  dY7 = ky/muy*z[31] - z[32]*ky/muy
  dY8 = ky/muy*z[32] - z[33]*ky/muy
  dY9 = ky/muy*z[33] - z[34]*ky/muy
  dY10 = ky/muy*z[34] - z[35]*ky/muy
  dY11 = ky/muy*z[35] - z[36]*ky/muy
  dY12 = ky/muy*z[36] - z[37]*ky/muy
  dY13 = ky/muy*z[37] - z[38]*ky/muy
  dY14 = ky/muy*z[38] - z[39]*ky/muy
  dY15 = ky/muy*z[39] - z[40]*ky/muy
  dY16 = ky/muy*z[40] - z[41]*ky/muy
  dY17 = ky/muy*z[41] - z[42]*ky/muy
  dY18 = ky/muy*z[42] - z[43]*ky/muy
  dY19 = ky/muy*z[43] - z[44]*ky/muy
  dY20 = ky/muy*z[44] - z[45]*ky/muy
  dY21 = ky/muy*z[45] - z[46]*ky/muy
  dY22 = ky/muy*z[46] - z[47]*ky/muy
  dY23 = ky/muy*z[47] - z[48]*ky/muy
  dY24 = ky/muy*z[48] - z[49]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dX14, dX15, dX16, dX17, dX18, dX19, dX20, dX21, dX22, dX23, dX24, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13, dY14, dY15, dY16, dY17, dY18, dY19, dY20, dY21, dY22, dY23, dY24)))
}

RMlct25.ode <- function(tm,z,ps) {
  r = ps[1]; K = ps[2]; a = ps[3]; h = ps[4]; chi = ps[5]; mux = ps[6]; muy = ps[7];
  #kx = ps[['kx']] 
  #ky = ps[['ky']] 
  
  N = z[1]
  P = onesyt %*% z[1+kx+(1:ky)] # total mature predators; sum Y_i
  
  dN  = r*N*(1-N/K) - a/(h+N)*P*N
  dX1 = chi*a/(h+N)*P*N - z[2]*kx/mux
  dX2 = kx/mux*z[2] - z[3]*kx/mux
  dX3 = kx/mux*z[3] - z[4]*kx/mux
  dX4 = kx/mux*z[4] - z[5]*kx/mux
  dX5 = kx/mux*z[5] - z[6]*kx/mux
  dX6 = kx/mux*z[6] - z[7]*kx/mux
  dX7 = kx/mux*z[7] - z[8]*kx/mux
  dX8 = kx/mux*z[8] - z[9]*kx/mux
  dX9 = kx/mux*z[9] - z[10]*kx/mux
  dX10 = kx/mux*z[10] - z[11]*kx/mux
  dX11 = kx/mux*z[11] - z[12]*kx/mux
  dX12 = kx/mux*z[12] - z[13]*kx/mux
  dX13 = kx/mux*z[13] - z[14]*kx/mux
  dX14 = kx/mux*z[14] - z[15]*kx/mux
  dX15 = kx/mux*z[15] - z[16]*kx/mux
  dX16 = kx/mux*z[16] - z[17]*kx/mux
  dX17 = kx/mux*z[17] - z[18]*kx/mux
  dX18 = kx/mux*z[18] - z[19]*kx/mux
  dX19 = kx/mux*z[19] - z[20]*kx/mux
  dX20 = kx/mux*z[20] - z[21]*kx/mux
  dX21 = kx/mux*z[21] - z[22]*kx/mux
  dX22 = kx/mux*z[22] - z[23]*kx/mux
  dX23 = kx/mux*z[23] - z[24]*kx/mux
  dX24 = kx/mux*z[24] - z[25]*kx/mux
  dX25 = kx/mux*z[25] - z[26]*kx/mux
  dY1 = kx/mux*z[26] - z[27]*ky/muy
  dY2 = ky/muy*z[27] - z[28]*ky/muy
  dY3 = ky/muy*z[28] - z[29]*ky/muy
  dY4 = ky/muy*z[29] - z[30]*ky/muy
  dY5 = ky/muy*z[30] - z[31]*ky/muy
  dY6 = ky/muy*z[31] - z[32]*ky/muy
  dY7 = ky/muy*z[32] - z[33]*ky/muy
  dY8 = ky/muy*z[33] - z[34]*ky/muy
  dY9 = ky/muy*z[34] - z[35]*ky/muy
  dY10 = ky/muy*z[35] - z[36]*ky/muy
  dY11 = ky/muy*z[36] - z[37]*ky/muy
  dY12 = ky/muy*z[37] - z[38]*ky/muy
  dY13 = ky/muy*z[38] - z[39]*ky/muy
  dY14 = ky/muy*z[39] - z[40]*ky/muy
  dY15 = ky/muy*z[40] - z[41]*ky/muy
  dY16 = ky/muy*z[41] - z[42]*ky/muy
  dY17 = ky/muy*z[42] - z[43]*ky/muy
  dY18 = ky/muy*z[43] - z[44]*ky/muy
  dY19 = ky/muy*z[44] - z[45]*ky/muy
  dY20 = ky/muy*z[45] - z[46]*ky/muy
  dY21 = ky/muy*z[46] - z[47]*ky/muy
  dY22 = ky/muy*z[47] - z[48]*ky/muy
  dY23 = ky/muy*z[48] - z[49]*ky/muy
  dY24 = ky/muy*z[49] - z[50]*ky/muy
  dY25 = ky/muy*z[50] - z[51]*ky/muy
  
  return(list(c(dN, dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11, dX12, dX13, dX14, dX15, dX16, dX17, dX18, dX19, dX20, dX21, dX22, dX23, dX24, dX25, dY1, dY2, dY3, dY4, dY5, dY6, dY7, dY8, dY9, dY10, dY11, dY12, dY13, dY14, dY15, dY16, dY17, dY18, dY19, dY20, dY21, dY22, dY23, dY24, dY25)))
}

################################################################################
## Benchmark...
################################################################################

parms <- params

RM.ode <- compiler::cmpfun(RM.ode)
RMpt.ode <- compiler::cmpfun(RMpt.ode)
RMlct1.ode <- compiler::cmpfun(RMlct1.ode)
RMlct2.ode <- compiler::cmpfun(RMlct2.ode)
RMlct3.ode <- compiler::cmpfun(RMlct3.ode)
RMlct4.ode <- compiler::cmpfun(RMlct4.ode)
RMlct5.ode <- compiler::cmpfun(RMlct5.ode)
RMlct6.ode <- compiler::cmpfun(RMlct6.ode)
RMlct7.ode <- compiler::cmpfun(RMlct7.ode)
RMlct8.ode <- compiler::cmpfun(RMlct8.ode)
RMlct9.ode <- compiler::cmpfun(RMlct9.ode)
RMlct10.ode <- compiler::cmpfun(RMlct10.ode)
RMlct11.ode <- compiler::cmpfun(RMlct11.ode)
RMlct12.ode <- compiler::cmpfun(RMlct12.ode)
RMlct13.ode <- compiler::cmpfun(RMlct13.ode)
RMlct14.ode <- compiler::cmpfun(RMlct14.ode)
RMlct15.ode <- compiler::cmpfun(RMlct15.ode)
RMlct16.ode <- compiler::cmpfun(RMlct16.ode)
RMlct17.ode <- compiler::cmpfun(RMlct17.ode)
RMlct18.ode <- compiler::cmpfun(RMlct18.ode)
RMlct19.ode <- compiler::cmpfun(RMlct19.ode)
RMlct20.ode <- compiler::cmpfun(RMlct20.ode)
RMlct21.ode <- compiler::cmpfun(RMlct21.ode)
RMlct22.ode <- compiler::cmpfun(RMlct22.ode)
RMlct23.ode <- compiler::cmpfun(RMlct23.ode)
RMlct24.ode <- compiler::cmpfun(RMlct24.ode)
RMlct25.ode <- compiler::cmpfun(RMlct25.ode)

mthd = "ode45"
atol= 1e-6
Tmax = 500
tms=seq(0,Tmax,length=500)


parms1 <- parms; parms1['kx']=1; parms1['ky']=1; RMpt.init(parms1); 
b1=benchmark(RM.1 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
             RMlct.1={ode(y=ICs, times=tms, func=RMlct1.ode, parms = parms1, method = mthd, atol=atol)},
             RMpt.1 ={ode(y=ICs, times=tms, func=RMpt.ode,   parms = parms1, method = mthd, atol=atol)}, 
             replications = reps)

parms2 <- parms; parms2['kx']=2; parms2['ky']=2; parms2; RMpt.init(parms2);
b2=benchmark(RM.2 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.2={ode(y=ICs, times=tms, func=RMlct2.ode, parms = parms2, method = mthd, atol=atol)},
          RMpt.2 ={ode(y=ICs, times=tms, func=RMpt.ode,   parms = parms2, method = mthd, atol=atol)}, 
          replications = reps)
          
parms3 <- parms; parms3['kx']=3; parms3['ky']=3; parms3; RMpt.init(parms3);
b3=benchmark(RM.3 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.3={ode(y=ICs, times=tms, func=RMlct3.ode, parms = parms3, method = mthd, atol=atol)},
          RMpt.3 ={ode(y=ICs, times=tms, func=RMpt.ode,   parms = parms3, method = mthd, atol=atol)}, 
          replications = reps)
          
parms4 <- parms; parms4['kx']=4; parms4['ky']=4; parms4; RMpt.init(parms4);
b4=benchmark(RM.4 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.4={ode(y=ICs, times=tms, func=RMlct4.ode, parms = parms4, method = mthd, atol=atol)},
          RMpt.4 ={ode(y=ICs, times=tms, func=RMpt.ode,   parms = parms4, method = mthd, atol=atol)}, 
          replications = reps)
          
parms5 <- parms; parms5['kx']=5; parms5['ky']=5; parms5; RMpt.init(parms5);
b5=benchmark(RM.5 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.5={ode(y=ICs, times=tms, func=RMlct5.ode, parms = parms5, method = mthd, atol=atol)},
          RMpt.5 ={ode(y=ICs, times=tms, func=RMpt.ode,   parms = parms5, method = mthd, atol=atol)}, 
          replications = reps)
          
parms6 <- parms; parms6['kx']=6; parms6['ky']=6; parms6; RMpt.init(parms6);
b6=benchmark(RM.6 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.6={ode(y=ICs, times=tms, func=RMlct6.ode, parms = parms6, method = mthd, atol=atol)},
          RMpt.6 ={ode(y=ICs, times=tms, func=RMpt.ode,   parms = parms6, method = mthd, atol=atol)}, 
          replications = reps)
          
parms7 <- parms; parms7['kx']=7; parms7['ky']=7; parms7; RMpt.init(parms7);
b7=benchmark(RM.7 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.7={ode(y=ICs, times=tms, func=RMlct7.ode, parms = parms7, method = mthd, atol=atol)},
          RMpt.7 ={ode(y=ICs, times=tms, func=RMpt.ode,   parms = parms7, method = mthd, atol=atol)}, 
          replications = reps)
          
parms8 <- parms; parms8['kx']=8; parms8['ky']=8; parms8; RMpt.init(parms8);
b8=benchmark(RM.8 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.8={ode(y=ICs, times=tms, func=RMlct8.ode, parms = parms8, method = mthd, atol=atol)},
          RMpt.8 ={ode(y=ICs, times=tms, func=RMpt.ode,   parms = parms8, method = mthd, atol=atol)}, 
          replications = reps) 
          
parms9 <- parms; parms9['kx']=9; parms9['ky']=9; parms9; RMpt.init(parms9);
b9=benchmark(RM.9 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.9={ode(y=ICs, times=tms, func=RMlct9.ode, parms = parms9, method = mthd, atol=atol)},
          RMpt.9 ={ode(y=ICs, times=tms, func=RMpt.ode,   parms = parms9, method = mthd, atol=atol)}, 
          replications = reps)    
          
parms10 <- parms; parms10['kx']=10; parms10['ky']=10; parms10; RMpt.init(parms10);
b10=benchmark(RM.10 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.10={ode(y=ICs, times=tms, func=RMlct10.ode, parms = parms10, method = mthd, atol=atol)},
          RMpt.10 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms10, method = mthd, atol=atol)}, 
          replications = reps) 
          
parms11 <- parms; parms11['kx']=11; parms11['ky']=11; parms11; RMpt.init(parms11);
b11=benchmark(RM.11 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.11={ode(y=ICs, times=tms, func=RMlct11.ode, parms = parms11, method = mthd, atol=atol)},
          RMpt.11 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms11, method = mthd, atol=atol)}, 
          replications = reps) 
          
parms12 <- parms; parms12['kx']=12; parms12['ky']=12; parms12; RMpt.init(parms12);
b12=benchmark(RM.12 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.12={ode(y=ICs, times=tms, func=RMlct12.ode, parms = parms12, method = mthd, atol=atol)},
          RMpt.12 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms12, method = mthd, atol=atol)}, 
          replications = reps) 
          
parms13 <- parms; parms13['kx']=13; parms13['ky']=13; parms13; RMpt.init(parms13);
b13=benchmark(RM.13 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.13={ode(y=ICs, times=tms, func=RMlct13.ode, parms = parms13, method = mthd, atol=atol)},
          RMpt.13 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms13, method = mthd, atol=atol)}, 
          replications = reps) 
          
parms14 <- parms; parms14['kx']=14; parms14['ky']=14; parms14; RMpt.init(parms14);
b14=benchmark(RM.14 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.14={ode(y=ICs, times=tms, func=RMlct14.ode, parms = parms14, method = mthd, atol=atol)},
          RMpt.14 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms14, method = mthd, atol=atol)}, 
          replications = reps) 
          
parms15 <- parms; parms15['kx']=15; parms15['ky']=15; parms15; RMpt.init(parms15);
b15=benchmark(RM.15 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.15={ode(y=ICs, times=tms, func=RMlct15.ode, parms = parms15, method = mthd, atol=atol)},
          RMpt.15 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms15, method = mthd, atol=atol)}, 
          replications = reps) 
          
parms16 <- parms; parms16['kx']=16; parms16['ky']=16; parms16; RMpt.init(parms16);
b16=benchmark(RM.16 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.16={ode(y=ICs, times=tms, func=RMlct16.ode, parms = parms16, method = mthd, atol=atol)},
          RMpt.16 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms16, method = mthd, atol=atol)}, 
          replications = reps) 
          
parms17 <- parms; parms17['kx']=17; parms17['ky']=17; parms17; RMpt.init(parms17);
b17=benchmark(RM.17 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.17={ode(y=ICs, times=tms, func=RMlct17.ode, parms = parms17, method = mthd, atol=atol)},
          RMpt.17 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms17, method = mthd, atol=atol)}, 
          replications = reps) 
          
parms18 <- parms; parms18['kx']=18; parms18['ky']=18; parms18; RMpt.init(parms18);
b18=benchmark(RM.18 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.18={ode(y=ICs, times=tms, func=RMlct18.ode, parms = parms18, method = mthd, atol=atol)},
          RMpt.18 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms18, method = mthd, atol=atol)}, 
          replications = reps) 
          
parms19 <- parms; parms19['kx']=19; parms19['ky']=19; parms19; RMpt.init(parms19);
b19=benchmark(RM.19 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.19={ode(y=ICs, times=tms, func=RMlct19.ode, parms = parms19, method = mthd, atol=atol)},
          RMpt.19 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms19, method = mthd, atol=atol)}, 
          replications = reps) 
          
parms20 <- parms; parms20['kx']=20; parms20['ky']=20; parms20; RMpt.init(parms20);
b20=benchmark(RM.20 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
          RMlct.20={ode(y=ICs, times=tms, func=RMlct20.ode, parms = parms20, method = mthd, atol=atol)},
          RMpt.20 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms20, method = mthd, atol=atol)}, 
          replications = reps) 
          
# parms21 <- parms; parms21['kx']=21; parms21['ky']=21; parms21; RMpt.init(parms21);
# b21=benchmark(RM.21 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
#           RMlct.21={ode(y=ICs, times=tms, func=RMlct21.ode, parms = parms21, method = mthd, atol=atol)},
#           RMpt.21 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms21, method = mthd, atol=atol)}, 
#           replications = reps)
#           
# parms22 <- parms; parms22['kx']=22; parms22['ky']=22; parms22; RMpt.init(parms22);
# b22=benchmark(RM.22 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
#           RMlct.22={ode(y=ICs, times=tms, func=RMlct22.ode, parms = parms22, method = mthd, atol=atol)},
#           RMpt.22 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms22, method = mthd, atol=atol)}, 
#           replications = reps)
#           
# parms23 <- parms; parms23['kx']=23; parms23['ky']=23; parms23; RMpt.init(parms23);
# b23=benchmark(RM.23 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
#           RMlct.23={ode(y=ICs, times=tms, func=RMlct23.ode, parms = parms23, method = mthd, atol=atol)},
#           RMpt.23 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms23, method = mthd, atol=atol)}, 
#           replications = reps)
#           
# parms24 <- parms; parms24['kx']=24; parms24['ky']=24; parms24; RMpt.init(parms24);
# b24=benchmark(RM.24 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
#           RMlct.24={ode(y=ICs, times=tms, func=RMlct24.ode, parms = parms24, method = mthd, atol=atol)},
#           RMpt.24 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms24, method = mthd, atol=atol)}, 
#           replications = reps)
#           
# parms25 <- parms; parms25['kx']=25; parms25['ky']=25; parms25; RMpt.init(parms25);
# b25=benchmark(RM.25 ={ode(y=IC, times = tms, func = RM.ode, parms = parms, method = mthd, atol=atol)}, 
#           RMlct.25={ode(y=ICs, times=tms, func=RMlct25.ode, parms = parms25, method = mthd, atol=atol)},
#           RMpt.25 ={ode(y=ICs, times=tms, func=RMpt.ode, parms = parms25, method = mthd, atol=atol)}, 
#           replications = reps)
                                                                                                                                                                                                                              
b1; b2; b3; b4; b5; b6; b7; b8; b9; b10; b11; b12; b13; b14; b15; b16; b17; b18; b19; b20;# b21; b22; b23; b24; b25

# First, reorganize the data and change up some labeling...
# 1. Rearrange rows to the right order, save a new copy...

out <- rbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20)#,b21,b22,b23,b24,b25)

out$N <- stringr::str_split_fixed(out$test,'\\.',2)[,2]
out$Model <- stringr::str_split_fixed(out$test,'\\.',2)[,1]
# out$Model <- gsub("RM$","RM w/o Maturation (Exponential)",out$Model)
# out$Model <- gsub("lct"," w/ Maturation (Erlang / LCT)",out$Model)
# out$Model <- gsub("pt"," w/ Maturation (Erlang as Phase-Type / GLCT)",out$Model)
out$Model <- gsub("RM$","RM w/ Standard: No maturation delay & exponential time as adult predators (M independent)",out$Model)
out$Model <- gsub("lct"," w/ LCT:    Erlang maturation delay & time as adult predators",out$Model)
out$Model <- gsub("pt"," w/ GLCT: Erlang maturation delay & time as adult predators (phase-type formulation)",out$Model)
out$Model <- gsub("RM w/ ",'',out$Model)
out$Model <- factor(out$Model, levels = unique(out$Model)[] )
out$N <- as.numeric(out$N)
head(out)

# 2. Plot using "N" to group the trios...

#win.graph(12,6); 
pdf("RMbenchmark.pdf",12,6);
ggplot(out, aes(x=factor(N), y=user.self, fill=Model, group=test)) + 
  geom_bar(stat="identity", position="dodge", width=0.6) + theme_bw(base_size=18) + 
  xlab("M (Model dimension is 1 + 2M)") + ylab("Elapsed Time (seconds)") + 
  scale_fill_grey(start=.75, end=0.2, name="Rosenzweig-MacArthur Model")+
  theme(legend.position=c(.375,.85), 
        legend.background = element_rect(fill="white", size = 0.5, color = "darkgray")) 
dev.off()

## Compare the trajectories for the two approaches to implementing the LCT
# outLCT  = ode(y=ICs, times=tms, func=RMlct2.ode, parms = parms2, method = mthd, atol=atol)
# outGLCT = ode(y=ICs, times=tms, func=RMpt.ode,   parms = parms2, method = mthd, atol=atol) 
# matplot(tms, cbind(N=outLCT[,2],X=colSums(t(outLCT[,3:4])),I=colSums(t(outLCT[,5:6]))),
#         type="l",lty=1, lwd=2, col=c("green","blue","red"),
#         xlab="Time", ylab="Population")
# matplot(tms, cbind(N=outGLCT[,2],X=colSums(t(outGLCT[,3:4])),I=colSums(t(outGLCT[,5:6]))),
#         type="l",lty=2, lwd=2, col=c("gray","gray","gray"), add=TRUE, 
#         xlab="Time", ylab="Population")


## Compare the trajectories for the two approaches to implementing the LCT
pdf("RM-compare.pdf",12,6);
  outLCT  = ode(y=ICs, times=tms, func=RMlct20.ode, parms = parms20, method = mthd, atol=atol)
  outGLCT = ode(y=ICs, times=tms, func=RMpt.ode,   parms = parms20, method = mthd, atol=atol) 
  matplot(tms, cbind(N=outLCT[,2],X=colSums(t(outLCT[,3:22])),I=colSums(t(outLCT[,23:42]))),
          type="l",lty=1, lwd=1, col="black", # c("green","blue","red"),
          xlab="Time", ylab="Population")
  matplot(tms, cbind(N=outGLCT[,2],X=colSums(t(outGLCT[,3:22])),I=colSums(t(outGLCT[,23:42]))),
          type="l",lty=3, lwd=1, col=c("gray","gray","gray"), add=TRUE, 
          xlab="Time", ylab="Population")
dev.off()


