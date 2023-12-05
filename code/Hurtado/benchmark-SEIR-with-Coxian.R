####################################################################################
## Benchmarking SEIR model numerical solutions
## ----------------------------------------------
##  Case 1: Exponential dwell time distributions (dimension = 3 if R omitted)
##  Case 2: Classic LCT (Coxian dwell times), hard-coded number of substates
##  Case 3: SEIR w/ phase-type distributions, parameterized w/ Coxian distributions
####################################################################################
#
#   NOTE: This R script was run on a 64-bit Windows 10 machine using R v3.6.3.
#         Running this on other operating systems may require small modifications
#         like replacing the win.graph() calls with x11() or quartz() calls.
#

library(deSolve)

## Parameterization....
reps<-500
parms = c(b=1, muE=4, muI=7, cvE=1, cvI=1, kE=1, kI=1, rhoE=0.99, rhoI=0.99) 


####################################################################################
## Function definitions

# SEIR with Exponential dwell times
SEIR.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I

  dS = -b*z[1]*z[3]
  dE =  b*z[1]*z[3] - z[2]/muE
  dI =  z[2]/muE - z[3]/muI
  dR = z[3]/muI
  
  return(list(c(dS, dE, dI, dR)))
}

# SEIR with hard-coded Coxian dwell times
SEIRlctCoxian1.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE =1
  kI =1
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot  - z[2]*kE/muE
  dI1 =  z[2]*kE/muE  - z[3]*kI/muI
  dR  =  z[3]*kI/muI
  
  return(list(c(dS, dE1, dI1, dR)))
}
SEIRlctCoxian2.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE=2
  kI=2
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
    
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE - z[3]*rE
  dI1 =  z[3]*rE + (1-rhoE)*z[2]*rE - z[4]*rI
  dI2 =  rhoI*z[4]*rI - z[5]*rI
  dR  =  z[5]*rI + (1-rhoI)*z[4]*rI
  
  return(list(c(dS, dE1, dE2, dI1, dI2, dR)))
}
SEIRlctCoxian3.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 3
  kI = 3
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dI1 =  z[4]*rE + (1-rhoE)*sum(z[2:3])*rE - z[5]*rI
  dI2 =  rhoI*z[5]*rI  - z[6]*rI
  dI3 =  rhoI*z[6]*rI  - z[7]*rI
  dR  =  z[7]*rI + (1-rhoI)*sum(z[5:6])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dI1, dI2, dI3, dR)))
}
SEIRlctCoxian4.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 4
  kI = 4
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dI1 =  z[5]*rE + (1-rhoE)*sum(z[2:4])*rE - z[6]*rI
  dI2 =  rhoI*z[6]*rI  - z[7]*rI
  dI3 =  rhoI*z[7]*rI  - z[8]*rI
  dI4 =  rhoI*z[8]*rI  - z[9]*rI
  dR  =  z[9]*rI + (1-rhoI)*sum(z[6:8])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dI1, dI2, dI3, dI4, dR)))
}
SEIRlctCoxian5.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 5
  kI = 5
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dI1 =  z[6]*rE + (1-rhoE)*sum(z[2:5])*rE - z[7]*rI
  dI2 =  rhoI*z[7]*rI  - z[8]*rI
  dI3 =  rhoI*z[8]*rI  - z[9]*rI
  dI4 =  rhoI*z[9]*rI  - z[10]*rI
  dI5 =  rhoI*z[10]*rI  - z[11]*rI
  dR  =  z[11]*rI + (1-rhoI)*sum(z[7:10])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dI1, dI2, dI3, dI4, dI5, dR)))
}
SEIRlctCoxian6.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 6
  kI = 6
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dI1 =  z[7]*rE + (1-rhoE)*sum(z[2:6])*rE - z[8]*rI
  dI2 =  rhoI*z[8]*rI  - z[9]*rI
  dI3 =  rhoI*z[9]*rI  - z[10]*rI
  dI4 =  rhoI*z[10]*rI  - z[11]*rI
  dI5 =  rhoI*z[11]*rI  - z[12]*rI
  dI6 =  rhoI*z[12]*rI  - z[13]*rI
  dR  =  z[13]*rI + (1-rhoI)*sum(z[8:12])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dI1, dI2, dI3, dI4, dI5, dI6, dR)))
}
SEIRlctCoxian7.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 7
  kI = 7
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dI1 =  z[8]*rE + (1-rhoE)*sum(z[2:7])*rE - z[9]*rI
  dI2 =  rhoI*z[9]*rI  - z[10]*rI
  dI3 =  rhoI*z[10]*rI  - z[11]*rI
  dI4 =  rhoI*z[11]*rI  - z[12]*rI
  dI5 =  rhoI*z[12]*rI  - z[13]*rI
  dI6 =  rhoI*z[13]*rI  - z[14]*rI
  dI7 =  rhoI*z[14]*rI  - z[15]*rI
  dR  =  z[15]*rI + (1-rhoI)*sum(z[9:14])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dR)))
}
SEIRlctCoxian8.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 8
  kI = 8
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dI1 =  z[9]*rE + (1-rhoE)*sum(z[2:8])*rE - z[10]*rI
  dI2 =  rhoI*z[10]*rI  - z[11]*rI
  dI3 =  rhoI*z[11]*rI  - z[12]*rI
  dI4 =  rhoI*z[12]*rI  - z[13]*rI
  dI5 =  rhoI*z[13]*rI  - z[14]*rI
  dI6 =  rhoI*z[14]*rI  - z[15]*rI
  dI7 =  rhoI*z[15]*rI  - z[16]*rI
  dI8 =  rhoI*z[16]*rI  - z[17]*rI
  dR  =  z[17]*rI + (1-rhoI)*sum(z[10:16])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dR)))
}
SEIRlctCoxian9.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 9
  kI = 9
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dI1 =  z[10]*rE + (1-rhoE)*sum(z[2:9])*rE - z[11]*rI
  dI2 =  rhoI*z[11]*rI  - z[12]*rI
  dI3 =  rhoI*z[12]*rI  - z[13]*rI
  dI4 =  rhoI*z[13]*rI  - z[14]*rI
  dI5 =  rhoI*z[14]*rI  - z[15]*rI
  dI6 =  rhoI*z[15]*rI  - z[16]*rI
  dI7 =  rhoI*z[16]*rI  - z[17]*rI
  dI8 =  rhoI*z[17]*rI  - z[18]*rI
  dI9 =  rhoI*z[18]*rI  - z[19]*rI
  dR  =  z[19]*rI + (1-rhoI)*sum(z[11:18])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dR)))
}

SEIRlctCoxian10.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 10
  kI = 10
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 = rhoE*z[10]*rE  - z[11]*rE
  dI1 =  z[11]*rE + (1-rhoE)*sum(z[2:10])*rE - z[12]*rI
  dI2 =  rhoI*z[12]*rI  - z[13]*rI
  dI3 =  rhoI*z[13]*rI  - z[14]*rI
  dI4 =  rhoI*z[14]*rI  - z[15]*rI
  dI5 =  rhoI*z[15]*rI  - z[16]*rI
  dI6 =  rhoI*z[16]*rI  - z[17]*rI
  dI7 =  rhoI*z[17]*rI  - z[18]*rI
  dI8 =  rhoI*z[18]*rI  - z[19]*rI
  dI9 =  rhoI*z[19]*rI  - z[20]*rI
  dI10 =  rhoI*z[20]*rI  - z[21]*rI
  dR  =  z[21]*rI + (1-rhoI)*sum(z[12:20])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dR)))
}

SEIRlctCoxian11.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 11
  kI = 11
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dI1 =  z[12]*rE + (1-rhoE)*sum(z[2:11])*rE  - z[13]*rI
  dI2 =  rhoI*z[13]*rI  - z[14]*rI
  dI3 =  rhoI*z[14]*rI  - z[15]*rI
  dI4 =  rhoI*z[15]*rI  - z[16]*rI
  dI5 =  rhoI*z[16]*rI  - z[17]*rI
  dI6 =  rhoI*z[17]*rI  - z[18]*rI
  dI7 =  rhoI*z[18]*rI  - z[19]*rI
  dI8 =  rhoI*z[19]*rI  - z[20]*rI
  dI9 =  rhoI*z[20]*rI  - z[21]*rI
  dI10 =  rhoI*z[21]*rI  - z[22]*rI
  dI11 =  rhoI*z[22]*rI  - z[23]*rI
  dR  =  z[23]*rI + (1-rhoI)*sum(z[13:22])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dR)))
}

SEIRlctCoxian12.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 12
  kI = 12
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dI1 =  z[13]*rE + (1-rhoE)*sum(z[2:12])*rE - z[14]*rI
  dI2 =  rhoI*z[14]*rI  - z[15]*rI
  dI3 =  rhoI*z[15]*rI  - z[16]*rI
  dI4 =  rhoI*z[16]*rI  - z[17]*rI
  dI5 =  rhoI*z[17]*rI  - z[18]*rI
  dI6 =  rhoI*z[18]*rI  - z[19]*rI
  dI7 =  rhoI*z[19]*rI  - z[20]*rI
  dI8 =  rhoI*z[20]*rI  - z[21]*rI
  dI9 =  rhoI*z[21]*rI  - z[22]*rI
  dI10 =  rhoI*z[22]*rI  - z[23]*rI
  dI11 =  rhoI*z[23]*rI  - z[24]*rI
  dI12 =  rhoI*z[24]*rI  - z[25]*rI
  dR  =  z[25]*rI + (1-rhoI)*sum(z[14:24])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dR)))
}

SEIRlctCoxian13.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 13
  kI = 13
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dI1 =  z[14]*rE + (1-rhoE)*sum(z[2:13])*rE - z[15]*rI
  dI2 =  rhoI*z[15]*rI  - z[16]*rI
  dI3 =  rhoI*z[16]*rI  - z[17]*rI
  dI4 =  rhoI*z[17]*rI  - z[18]*rI
  dI5 =  rhoI*z[18]*rI  - z[19]*rI
  dI6 =  rhoI*z[19]*rI  - z[20]*rI
  dI7 =  rhoI*z[20]*rI  - z[21]*rI
  dI8 =  rhoI*z[21]*rI  - z[22]*rI
  dI9 =  rhoI*z[22]*rI  - z[23]*rI
  dI10 =  rhoI*z[23]*rI  - z[24]*rI
  dI11 =  rhoI*z[24]*rI  - z[25]*rI
  dI12 =  rhoI*z[25]*rI  - z[26]*rI
  dI13 =  rhoI*z[26]*rI  - z[27]*rI
  dR  =  z[27]*rI + (1-rhoI)*sum(z[15:26])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dR)))
}

SEIRlctCoxian14.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 14
  kI = 14
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dE14 =  rhoE*z[14]*rE  - z[15]*rE
  dI1 =  z[15]*rE + (1-rhoE)*sum(z[2:14])*rE - z[16]*rI
  dI2 =  rhoI*z[16]*rI  - z[17]*rI
  dI3 =  rhoI*z[17]*rI  - z[18]*rI
  dI4 =  rhoI*z[18]*rI  - z[19]*rI
  dI5 =  rhoI*z[19]*rI  - z[20]*rI
  dI6 =  rhoI*z[20]*rI  - z[21]*rI
  dI7 =  rhoI*z[21]*rI  - z[22]*rI
  dI8 =  rhoI*z[22]*rI  - z[23]*rI
  dI9 =  rhoI*z[23]*rI  - z[24]*rI
  dI10 =  rhoI*z[24]*rI  - z[25]*rI
  dI11 =  rhoI*z[25]*rI  - z[26]*rI
  dI12 =  rhoI*z[26]*rI  - z[27]*rI
  dI13 =  rhoI*z[27]*rI  - z[28]*rI
  dI14 =  rhoI*z[28]*rI  - z[29]*rI
  dR  =  z[29]*rI + (1-rhoI)*sum(z[16:28])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dR)))
}

SEIRlctCoxian15.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 15
  kI = 15
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dE14 =  rhoE*z[14]*rE  - z[15]*rE
  dE15 =  rhoE*z[15]*rE  - z[16]*rE
  dI1 =  z[16]*rE + (1-rhoE)*sum(z[2:15])*rE - z[17]*rI
  dI2 =  rhoI*z[17]*rI  - z[18]*rI
  dI3 =  rhoI*z[18]*rI  - z[19]*rI
  dI4 =  rhoI*z[19]*rI  - z[20]*rI
  dI5 =  rhoI*z[20]*rI  - z[21]*rI
  dI6 =  rhoI*z[21]*rI  - z[22]*rI
  dI7 =  rhoI*z[22]*rI  - z[23]*rI
  dI8 =  rhoI*z[23]*rI  - z[24]*rI
  dI9 =  rhoI*z[24]*rI  - z[25]*rI
  dI10 =  rhoI*z[25]*rI  - z[26]*rI
  dI11 =  rhoI*z[26]*rI  - z[27]*rI
  dI12 =  rhoI*z[27]*rI  - z[28]*rI
  dI13 =  rhoI*z[28]*rI  - z[29]*rI
  dI14 =  rhoI*z[29]*rI  - z[30]*rI
  dI15 =  rhoI*z[30]*rI  - z[31]*rI
  dR  =  z[31]*rI + (1-rhoI)*sum(z[17:30])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dR)))
}

SEIRlctCoxian16.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 16
  kI = 16
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dE14 =  rhoE*z[14]*rE  - z[15]*rE
  dE15 =  rhoE*z[15]*rE  - z[16]*rE
  dE16 =  rhoE*z[16]*rE  - z[17]*rE
  dI1 =  z[17]*rE + (1-rhoE)*sum(z[2:16])*rE - z[18]*rI
  dI2 =  rhoI*z[18]*rI  - z[19]*rI
  dI3 =  rhoI*z[19]*rI  - z[20]*rI
  dI4 =  rhoI*z[20]*rI  - z[21]*rI
  dI5 =  rhoI*z[21]*rI  - z[22]*rI
  dI6 =  rhoI*z[22]*rI  - z[23]*rI
  dI7 =  rhoI*z[23]*rI  - z[24]*rI
  dI8 =  rhoI*z[24]*rI  - z[25]*rI
  dI9 =  rhoI*z[25]*rI  - z[26]*rI
  dI10 =  rhoI*z[26]*rI  - z[27]*rI
  dI11 =  rhoI*z[27]*rI  - z[28]*rI
  dI12 =  rhoI*z[28]*rI  - z[29]*rI
  dI13 =  rhoI*z[29]*rI  - z[30]*rI
  dI14 =  rhoI*z[30]*rI  - z[31]*rI
  dI15 =  rhoI*z[31]*rI  - z[32]*rI
  dI16 =  rhoI*z[32]*rI  - z[33]*rI
  dR  =  z[33]*rI + (1-rhoI)*sum(z[18:32])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dR)))
}

SEIRlctCoxian17.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 17
  kI = 17
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dE14 =  rhoE*z[14]*rE  - z[15]*rE
  dE15 =  rhoE*z[15]*rE  - z[16]*rE
  dE16 =  rhoE*z[16]*rE  - z[17]*rE
  dE17 =  rhoE*z[17]*rE  - z[18]*rE
  dI1 =  z[18]*rE + (1-rhoE)*sum(z[2:17])*rE - z[19]*rI
  dI2 =  rhoI*z[19]*rI  - z[20]*rI
  dI3 =  rhoI*z[20]*rI  - z[21]*rI
  dI4 =  rhoI*z[21]*rI  - z[22]*rI
  dI5 =  rhoI*z[22]*rI  - z[23]*rI
  dI6 =  rhoI*z[23]*rI  - z[24]*rI
  dI7 =  rhoI*z[24]*rI  - z[25]*rI
  dI8 =  rhoI*z[25]*rI  - z[26]*rI
  dI9 =  rhoI*z[26]*rI  - z[27]*rI
  dI10 =  rhoI*z[27]*rI  - z[28]*rI
  dI11 =  rhoI*z[28]*rI  - z[29]*rI
  dI12 =  rhoI*z[29]*rI  - z[30]*rI
  dI13 =  rhoI*z[30]*rI  - z[31]*rI
  dI14 =  rhoI*z[31]*rI  - z[32]*rI
  dI15 =  rhoI*z[32]*rI  - z[33]*rI
  dI16 =  rhoI*z[33]*rI  - z[34]*rI
  dI17 =  rhoI*z[34]*rI  - z[35]*rI
  dR  =  z[35]*rI + (1-rhoI)*sum(z[19:34])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dR)))
}

SEIRlctCoxian18.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 18
  kI = 18
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dE14 =  rhoE*z[14]*rE  - z[15]*rE
  dE15 =  rhoE*z[15]*rE  - z[16]*rE
  dE16 =  rhoE*z[16]*rE  - z[17]*rE
  dE17 =  rhoE*z[17]*rE  - z[18]*rE
  dE18 =  rhoE*z[18]*rE  - z[19]*rE
  dI1 =  z[19]*rE + (1-rhoE)*sum(z[2:18])*rE - z[20]*rI
  dI2 =  rhoI*z[20]*rI  - z[21]*rI
  dI3 =  rhoI*z[21]*rI  - z[22]*rI
  dI4 =  rhoI*z[22]*rI  - z[23]*rI
  dI5 =  rhoI*z[23]*rI  - z[24]*rI
  dI6 =  rhoI*z[24]*rI  - z[25]*rI
  dI7 =  rhoI*z[25]*rI  - z[26]*rI
  dI8 =  rhoI*z[26]*rI  - z[27]*rI
  dI9 =  rhoI*z[27]*rI  - z[28]*rI
  dI10 =  rhoI*z[28]*rI  - z[29]*rI
  dI11 =  rhoI*z[29]*rI  - z[30]*rI
  dI12 =  rhoI*z[30]*rI  - z[31]*rI
  dI13 =  rhoI*z[31]*rI  - z[32]*rI
  dI14 =  rhoI*z[32]*rI  - z[33]*rI
  dI15 =  rhoI*z[33]*rI  - z[34]*rI
  dI16 =  rhoI*z[34]*rI  - z[35]*rI
  dI17 =  rhoI*z[35]*rI  - z[36]*rI
  dI18 =  rhoI*z[36]*rI  - z[37]*rI
  dR  =  z[37]*rI + (1-rhoI)*sum(z[20:36])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dR)))
}

SEIRlctCoxian19.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 19
  kI = 19
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dE14 =  rhoE*z[14]*rE  - z[15]*rE
  dE15 =  rhoE*z[15]*rE  - z[16]*rE
  dE16 =  rhoE*z[16]*rE  - z[17]*rE
  dE17 =  rhoE*z[17]*rE  - z[18]*rE
  dE18 =  rhoE*z[18]*rE  - z[19]*rE
  dE19 =  rhoE*z[19]*rE  - z[20]*rE
  dI1 =  z[20]*rE + (1-rhoE)*sum(z[2:19])*rE - z[21]*rI
  dI2 =  rhoI*z[21]*rI  - z[22]*rI
  dI3 =  rhoI*z[22]*rI  - z[23]*rI
  dI4 =  rhoI*z[23]*rI  - z[24]*rI
  dI5 =  rhoI*z[24]*rI  - z[25]*rI
  dI6 =  rhoI*z[25]*rI  - z[26]*rI
  dI7 =  rhoI*z[26]*rI  - z[27]*rI
  dI8 =  rhoI*z[27]*rI  - z[28]*rI
  dI9 =  rhoI*z[28]*rI  - z[29]*rI
  dI10 =  rhoI*z[29]*rI  - z[30]*rI
  dI11 =  rhoI*z[30]*rI  - z[31]*rI
  dI12 =  rhoI*z[31]*rI  - z[32]*rI
  dI13 =  rhoI*z[32]*rI  - z[33]*rI
  dI14 =  rhoI*z[33]*rI  - z[34]*rI
  dI15 =  rhoI*z[34]*rI  - z[35]*rI
  dI16 =  rhoI*z[35]*rI  - z[36]*rI
  dI17 =  rhoI*z[36]*rI  - z[37]*rI
  dI18 =  rhoI*z[37]*rI  - z[38]*rI
  dI19 =  rhoI*z[38]*rI  - z[39]*rI
  dR  =  z[39]*rI + (1-rhoI)*sum(z[21:38])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dR)))
}

SEIRlctCoxian20.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 20
  kI = 20
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dE14 =  rhoE*z[14]*rE  - z[15]*rE
  dE15 =  rhoE*z[15]*rE  - z[16]*rE
  dE16 =  rhoE*z[16]*rE  - z[17]*rE
  dE17 =  rhoE*z[17]*rE  - z[18]*rE
  dE18 =  rhoE*z[18]*rE  - z[19]*rE
  dE19 =  rhoE*z[19]*rE  - z[20]*rE
  dE20 =  rhoE*z[20]*rE  - z[21]*rE
  dI1 =  z[21]*rE + (1-rhoE)*sum(z[2:20])*rE - z[22]*rI
  dI2 =  rhoI*z[22]*rI  - z[23]*rI
  dI3 =  rhoI*z[23]*rI  - z[24]*rI
  dI4 =  rhoI*z[24]*rI  - z[25]*rI
  dI5 =  rhoI*z[25]*rI  - z[26]*rI
  dI6 =  rhoI*z[26]*rI  - z[27]*rI
  dI7 =  rhoI*z[27]*rI  - z[28]*rI
  dI8 =  rhoI*z[28]*rI  - z[29]*rI
  dI9 =  rhoI*z[29]*rI  - z[30]*rI
  dI10 =  rhoI*z[30]*rI  - z[31]*rI
  dI11 =  rhoI*z[31]*rI  - z[32]*rI
  dI12 =  rhoI*z[32]*rI  - z[33]*rI
  dI13 =  rhoI*z[33]*rI  - z[34]*rI
  dI14 =  rhoI*z[34]*rI  - z[35]*rI
  dI15 =  rhoI*z[35]*rI  - z[36]*rI
  dI16 =  rhoI*z[36]*rI  - z[37]*rI
  dI17 =  rhoI*z[37]*rI  - z[38]*rI
  dI18 =  rhoI*z[38]*rI  - z[39]*rI
  dI19 =  rhoI*z[39]*rI  - z[40]*rI
  dI20 =  rhoI*z[40]*rI  - z[41]*rI
  dR  =  z[41]*rI + (1-rhoI)*sum(z[22:40])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dE20, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI20, dR)))
}

SEIRlctCoxian21.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 21
  kI = 21
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dE14 =  rhoE*z[14]*rE  - z[15]*rE
  dE15 =  rhoE*z[15]*rE  - z[16]*rE
  dE16 =  rhoE*z[16]*rE  - z[17]*rE
  dE17 =  rhoE*z[17]*rE  - z[18]*rE
  dE18 =  rhoE*z[18]*rE  - z[19]*rE
  dE19 =  rhoE*z[19]*rE  - z[20]*rE
  dE20 =  rhoE*z[20]*rE  - z[21]*rE
  dE21 =  rhoE*z[21]*rE  - z[22]*rE
  dI1 =  z[22]*rE + (1-rhoE)*sum(z[2:21])*rE - z[23]*rI
  dI2 =  rhoI*z[23]*rI  - z[24]*rI
  dI3 =  rhoI*z[24]*rI  - z[25]*rI
  dI4 =  rhoI*z[25]*rI  - z[26]*rI
  dI5 =  rhoI*z[26]*rI  - z[27]*rI
  dI6 =  rhoI*z[27]*rI  - z[28]*rI
  dI7 =  rhoI*z[28]*rI  - z[29]*rI
  dI8 =  rhoI*z[29]*rI  - z[30]*rI
  dI9 =  rhoI*z[30]*rI  - z[31]*rI
  dI10 =  rhoI*z[31]*rI  - z[32]*rI
  dI11 =  rhoI*z[32]*rI  - z[33]*rI
  dI12 =  rhoI*z[33]*rI  - z[34]*rI
  dI13 =  rhoI*z[34]*rI  - z[35]*rI
  dI14 =  rhoI*z[35]*rI  - z[36]*rI
  dI15 =  rhoI*z[36]*rI  - z[37]*rI
  dI16 =  rhoI*z[37]*rI  - z[38]*rI
  dI17 =  rhoI*z[38]*rI  - z[39]*rI
  dI18 =  rhoI*z[39]*rI  - z[40]*rI
  dI19 =  rhoI*z[40]*rI  - z[41]*rI
  dI20 =  rhoI*z[41]*rI  - z[42]*rI
  dI21 =  rhoI*z[42]*rI  - z[43]*rI
  dR  =  z[43]*rI + (1-rhoI)*sum(z[23:42])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dE20, dE21, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI20, dI21, dR)))
}

SEIRlctCoxian22.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 22
  kI = 22
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 = rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dE14 =  rhoE*z[14]*rE  - z[15]*rE
  dE15 =  rhoE*z[15]*rE  - z[16]*rE
  dE16 =  rhoE*z[16]*rE  - z[17]*rE
  dE17 =  rhoE*z[17]*rE  - z[18]*rE
  dE18 =  rhoE*z[18]*rE  - z[19]*rE
  dE19 =  rhoE*z[19]*rE  - z[20]*rE
  dE20 =  rhoE*z[20]*rE  - z[21]*rE
  dE21 =  rhoE*z[21]*rE  - z[22]*rE
  dE22 =  rhoE*z[22]*rE  - z[23]*rE
  dI1 =  z[23]*rE + (1-rhoE)*sum(z[2:22])*rE - z[24]*rI
  dI2 =  rhoI*z[24]*rI  - z[25]*rI
  dI3 =  rhoI*z[25]*rI  - z[26]*rI
  dI4 =  rhoI*z[26]*rI  - z[27]*rI
  dI5 =  rhoI*z[27]*rI  - z[28]*rI
  dI6 =  rhoI*z[28]*rI  - z[29]*rI
  dI7 =  rhoI*z[29]*rI  - z[30]*rI
  dI8 =  rhoI*z[30]*rI  - z[31]*rI
  dI9 =  rhoI*z[31]*rI  - z[32]*rI
  dI10 =  rhoI*z[32]*rI  - z[33]*rI
  dI11 =  rhoI*z[33]*rI  - z[34]*rI
  dI12 =  rhoI*z[34]*rI  - z[35]*rI
  dI13 =  rhoI*z[35]*rI  - z[36]*rI
  dI14 =  rhoI*z[36]*rI  - z[37]*rI
  dI15 =  rhoI*z[37]*rI  - z[38]*rI
  dI16 =  rhoI*z[38]*rI  - z[39]*rI
  dI17 =  rhoI*z[39]*rI  - z[40]*rI
  dI18 =  rhoI*z[40]*rI  - z[41]*rI
  dI19 =  rhoI*z[41]*rI  - z[42]*rI
  dI20 =  rhoI*z[42]*rI  - z[43]*rI
  dI21 =  rhoI*z[43]*rI  - z[44]*rI
  dI22 =  rhoI*z[44]*rI  - z[45]*rI
  dR  =  z[45]*rI + (1-rhoI)*sum(z[24:44])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dE20, dE21, dE22, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI20, dI21, dI22, dR)))
}

SEIRlctCoxian23.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 23
  kI = 23
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dE14 =  rhoE*z[14]*rE  - z[15]*rE
  dE15 =  rhoE*z[15]*rE  - z[16]*rE
  dE16 =  rhoE*z[16]*rE  - z[17]*rE
  dE17 =  rhoE*z[17]*rE  - z[18]*rE
  dE18 =  rhoE*z[18]*rE  - z[19]*rE
  dE19 =  rhoE*z[19]*rE  - z[20]*rE
  dE20 =  rhoE*z[20]*rE  - z[21]*rE
  dE21 =  rhoE*z[21]*rE  - z[22]*rE
  dE22 =  rhoE*z[22]*rE  - z[23]*rE
  dE23 =  rhoE*z[23]*rE  - z[24]*rE
  dI1 =  z[24]*rE + (1-rhoE)*sum(z[2:23])*rE - z[25]*rI
  dI2 =  rhoI*z[25]*rI  - z[26]*rI
  dI3 =  rhoI*z[26]*rI  - z[27]*rI
  dI4 =  rhoI*z[27]*rI  - z[28]*rI
  dI5 =  rhoI*z[28]*rI  - z[29]*rI
  dI6 =  rhoI*z[29]*rI  - z[30]*rI
  dI7 =  rhoI*z[30]*rI  - z[31]*rI
  dI8 =  rhoI*z[31]*rI  - z[32]*rI
  dI9 =  rhoI*z[32]*rI  - z[33]*rI
  dI10 =  rhoI*z[33]*rI  - z[34]*rI
  dI11 =  rhoI*z[34]*rI  - z[35]*rI
  dI12 =  rhoI*z[35]*rI  - z[36]*rI
  dI13 =  rhoI*z[36]*rI  - z[37]*rI
  dI14 =  rhoI*z[37]*rI  - z[38]*rI
  dI15 =  rhoI*z[38]*rI  - z[39]*rI
  dI16 =  rhoI*z[39]*rI  - z[40]*rI
  dI17 =  rhoI*z[40]*rI  - z[41]*rI
  dI18 =  rhoI*z[41]*rI  - z[42]*rI
  dI19 =  rhoI*z[42]*rI  - z[43]*rI
  dI20 =  rhoI*z[43]*rI  - z[44]*rI
  dI21 =  rhoI*z[44]*rI  - z[45]*rI
  dI22 =  rhoI*z[45]*rI  - z[46]*rI
  dI23 =  rhoI*z[46]*rI  - z[47]*rI
  dR  =  z[47]*rI + (1-rhoI)*sum(z[25:46])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dE20, dE21, dE22, dE23, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI20, dI21, dI22, dI23, dR)))
}

SEIRlctCoxian24.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 24
  kI = 24
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dE14 =  rhoE*z[14]*rE  - z[15]*rE
  dE15 =  rhoE*z[15]*rE  - z[16]*rE
  dE16 =  rhoE*z[16]*rE  - z[17]*rE
  dE17 =  rhoE*z[17]*rE  - z[18]*rE
  dE18 =  rhoE*z[18]*rE  - z[19]*rE
  dE19 =  rhoE*z[19]*rE  - z[20]*rE
  dE20 =  rhoE*z[20]*rE  - z[21]*rE
  dE21 =  rhoE*z[21]*rE  - z[22]*rE
  dE22 =  rhoE*z[22]*rE  - z[23]*rE
  dE23 =  rhoE*z[23]*rE  - z[24]*rE
  dE24 =  rhoE*z[24]*rE  - z[25]*rE
  dI1 =  z[25]*rE + (1-rhoE)*sum(z[2:24])*rE - z[26]*rI
  dI2 =  rhoI*z[26]*rI  - z[27]*rI
  dI3 =  rhoI*z[27]*rI  - z[28]*rI
  dI4 =  rhoI*z[28]*rI  - z[29]*rI
  dI5 =  rhoI*z[29]*rI  - z[30]*rI
  dI6 =  rhoI*z[30]*rI  - z[31]*rI
  dI7 =  rhoI*z[31]*rI  - z[32]*rI
  dI8 =  rhoI*z[32]*rI  - z[33]*rI
  dI9 =  rhoI*z[33]*rI  - z[34]*rI
  dI10 =  rhoI*z[34]*rI  - z[35]*rI
  dI11 =  rhoI*z[35]*rI  - z[36]*rI
  dI12 =  rhoI*z[36]*rI  - z[37]*rI
  dI13 =  rhoI*z[37]*rI  - z[38]*rI
  dI14 =  rhoI*z[38]*rI  - z[39]*rI
  dI15 =  rhoI*z[39]*rI  - z[40]*rI
  dI16 =  rhoI*z[40]*rI  - z[41]*rI
  dI17 =  rhoI*z[41]*rI  - z[42]*rI
  dI18 =  rhoI*z[42]*rI  - z[43]*rI
  dI19 =  rhoI*z[43]*rI  - z[44]*rI
  dI20 =  rhoI*z[44]*rI  - z[45]*rI
  dI21 =  rhoI*z[45]*rI  - z[46]*rI
  dI22 =  rhoI*z[46]*rI  - z[47]*rI
  dI23 =  rhoI*z[47]*rI  - z[48]*rI
  dI24 =  rhoI*z[48]*rI  - z[49]*rI
  dR  =  z[49]*rI + (1-rhoI)*sum(z[26:48])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dE20, dE21, dE22, dE23, dE24, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI20, dI21, dI22, dI23, dI24, dR)))
}

SEIRlctCoxian25.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 25
  kI = 25
  rhoE=ps[["rhoE"]]# parameter for Coxian distribution. rhoE=1 yields Erlang case
  rhoI=ps[["rhoI"]]# parameter for Coxian distribution. rhoI=1 yields Erlang case
  # Coxian mean muE = (1 + p + ... + p^kE-1))/rE where rE is the loss rate from each state.
  # Rearranging this, we can fix the mean at muE by definint rE appropriately:
  rE=sum(rhoE^(0:(kE-1)))/muE # Coxian distribution w/ this rate yields mean=muE
  rI=sum(rhoI^(0:(kI-1)))/muI 
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*rE
  dE2 =  rhoE*z[2]*rE  - z[3]*rE
  dE3 =  rhoE*z[3]*rE  - z[4]*rE
  dE4 =  rhoE*z[4]*rE  - z[5]*rE
  dE5 =  rhoE*z[5]*rE  - z[6]*rE
  dE6 =  rhoE*z[6]*rE  - z[7]*rE
  dE7 =  rhoE*z[7]*rE  - z[8]*rE
  dE8 =  rhoE*z[8]*rE  - z[9]*rE
  dE9 =  rhoE*z[9]*rE  - z[10]*rE
  dE10 =  rhoE*z[10]*rE  - z[11]*rE
  dE11 =  rhoE*z[11]*rE  - z[12]*rE
  dE12 =  rhoE*z[12]*rE  - z[13]*rE
  dE13 =  rhoE*z[13]*rE  - z[14]*rE
  dE14 =  rhoE*z[14]*rE  - z[15]*rE
  dE15 =  rhoE*z[15]*rE  - z[16]*rE
  dE16 =  rhoE*z[16]*rE  - z[17]*rE
  dE17 =  rhoE*z[17]*rE  - z[18]*rE
  dE18 =  rhoE*z[18]*rE  - z[19]*rE
  dE19 =  rhoE*z[19]*rE  - z[20]*rE
  dE20 =  rhoE*z[20]*rE  - z[21]*rE
  dE21 =  rhoE*z[21]*rE  - z[22]*rE
  dE22 =  rhoE*z[22]*rE  - z[23]*rE
  dE23 =  rhoE*z[23]*rE  - z[24]*rE
  dE24 =  rhoE*z[24]*rE  - z[25]*rE
  dE25 =  rhoE*z[25]*rE  - z[26]*rE
  dI1 =  z[26]*rE + (1-rhoE)*sum(z[2:25])*rE - z[27]*rI
  dI2 =  rhoI*z[27]*rI  - z[28]*rI
  dI3 =  rhoI*z[28]*rI  - z[29]*rI
  dI4 =  rhoI*z[29]*rI  - z[30]*rI
  dI5 =  rhoI*z[30]*rI  - z[31]*rI
  dI6 =  rhoI*z[31]*rI  - z[32]*rI
  dI7 =  rhoI*z[32]*rI  - z[33]*rI
  dI8 =  rhoI*z[33]*rI  - z[34]*rI
  dI9 =  rhoI*z[34]*rI  - z[35]*rI
  dI10 =  rhoI*z[35]*rI  - z[36]*rI
  dI11 =  rhoI*z[36]*rI  - z[37]*rI
  dI12 =  rhoI*z[37]*rI  - z[38]*rI
  dI13 =  rhoI*z[38]*rI  - z[39]*rI
  dI14 =  rhoI*z[39]*rI  - z[40]*rI
  dI15 =  rhoI*z[40]*rI  - z[41]*rI
  dI16 =  rhoI*z[41]*rI  - z[42]*rI
  dI17 =  rhoI*z[42]*rI  - z[43]*rI
  dI18 =  rhoI*z[43]*rI  - z[44]*rI
  dI19 =  rhoI*z[44]*rI  - z[45]*rI
  dI20 =  rhoI*z[45]*rI  - z[46]*rI
  dI21 =  rhoI*z[46]*rI  - z[47]*rI
  dI22 =  rhoI*z[47]*rI  - z[48]*rI
  dI23 =  rhoI*z[48]*rI  - z[49]*rI
  dI24 =  rhoI*z[49]*rI  - z[50]*rI
  dI25 =  rhoI*z[50]*rI  - z[51]*rI
  dR  =  z[51]*rI + (1-rhoI)*sum(z[27:50])*rI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dE20, dE21, dE22, dE23, dE24, dE25, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI20, dI21, dI22, dI23, dI24, dI25, dR)))
}

# SEIR with Coxian dwell times in a GLCT framework
SEIRpt.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE =ps[['kE']]
  kI =ps[['kI']]
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  dS = -b*z[1]*Itot
  dE =  b*z[1]*Itot                        *  aE + AEt %*% z[1+(1:kE)]  # t(AE) %*% Evec
  dI =  aI %*% (-OnesE%*%AEt %*% z[1+(1:kE)]) + AIt %*% z[1+kE+(1:kI)] #  t(AI) %*% Ivec
  dR = as.numeric(-OnesI%*%AIt %*% z[1+kE+(1:kI)]) 
  
  return(list(c(dS, as.numeric(dE), as.numeric(dI), dR)))
}

SEIRpt.init <- function(ps) {

  #print("Defining ICs, AEt, AIt, OnesE, and OnesI ... ")
  
  # Unpack some parameter values...
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = ps[['kE']] # Number of substates in E, and...
  kI = ps[['kI']] # ... I.
  rhoE=ps[['rhoE']] # fraction advancing to the next "phase" (for Coxian distribution)
  rhoI=ps[['rhoI']]
  
  # These are Coxian distributions framed in a Phase-type distribution context,
  # where vector a = (1 0 ... 0) and matrix A is as follows...
  
  aE = matrix(0,nrow = kE, ncol=1); aE[1] = 1;
  AE = sum(rhoE^(0:(kE-1)))/muE*(diag(rep(-1,kE),kE)); 
  if(kE>1) for(i in 1:(kE-1)) {AE[i,i+1] = rhoE*sum(rhoE^(0:(kE-1)))/muE}
  
  aI = matrix(0,nrow = kI, ncol=1); aI[1] = 1;
  AI = sum(rhoI^(0:(kI-1)))/muI*(diag(rep(-1,kI),kI)); 
  if(kI>1) for(i in 1:(kI-1)) {AI[i,i+1] = rhoI*sum(rhoI^(0:(kI-1)))/muI}
  
  
  # Initial conditions
  PopSize=10000
  z0=numeric(kE+kI+2) # initialize the 1+kE+kI+1 state variables
  z0[2] <- 1/PopSize  # 1 in the initial exposed class 
  z0[1] <- 1-z0[2]    # susceptibles = (PopSize - 1)/PopSize
  
  # Set some global variables...
  aE  <<- aE
  AEt <<- t(AE)
  aI  <<- aI
  AIt <<- t(AI)
  OnesE <<- matrix(1,ncol=kE,nrow=1)
  OnesI <<- matrix(1,ncol=kI,nrow=1)
  ICs <<- z0
}

############################################################################
## Example plot

library(ggplot2)

Tmax = 100
tms=seq(0,Tmax,length=200)
SEIRpt.init(parms) # set's some initial conditions in ICs, and some matrices needed for SEIRpt.ode()

out = ode(y=ICs, times = tms, func = SEIRpt.ode, parms = parms, method = "ode23");
c(Pos=sum(out[,-1]>=0), Neg=sum(out[,-1]<0) )

matplot(tms, out[,-1],type="l",lty=1)

out = ode(y=ICs, times = tms, func = SEIRpt.ode, parms = parms, method = "ode45");
c(Pos=sum(out[,-1]>=0), Neg=sum(out[,-1]<0) )

matplot(tms, out[,-1],type="l",lty=2,add=TRUE)

out = ode(y=ICs, times = tms, func = SEIRpt.ode, parms = parms, method = "lsodes");
c(Pos=sum(out[,-1]>=0), Neg=sum(out[,-1]<0) )

matplot(tms, out[,-1],type="l",lty=2,add=TRUE)

out = ode(y=ICs, times = tms, func = SEIRpt.ode, parms = parms, method = "vode");
c(Pos=sum(out[,-1]>=0), Neg=sum(out[,-1]<0) )

matplot(tms, out[,-1],type="l",lty=2,add=TRUE)

# Some basic benchmarking...
library(rbenchmark)
benchmark(
  ode23 = ode(y=ICs, times = tms, func = SEIRpt.ode, parms = parms, method = "ode23", atol=1e-6),
  ode45 = ode(y=ICs, times = tms, func = SEIRpt.ode, parms = parms, method = "ode45", atol=1e-6),
  lsoda = ode(y=ICs, times = tms, func = SEIRpt.ode, parms = parms, method = "lsoda", atol=1e-6),
  lsode = ode(y=ICs, times = tms, func = SEIRpt.ode, parms = parms, method = "lsode", atol=1e-6),
  lsodes = ode(y=ICs, times = tms, func = SEIRpt.ode, parms = parms, method = "lsodes", atol=1e-6),
  vode = ode(y=ICs, times = tms, func = SEIRpt.ode, parms = parms, method = "vode", atol=1e-6),
  replications = 100)

############################################################################
## Benchmark...

SEIR.ode <- compiler::cmpfun(SEIR.ode)
SEIRpt.ode <- compiler::cmpfun(SEIRpt.ode)
SEIRlctCoxian1.ode <- compiler::cmpfun(SEIRlctCoxian1.ode)
SEIRlctCoxian2.ode <- compiler::cmpfun(SEIRlctCoxian2.ode)
SEIRlctCoxian3.ode <- compiler::cmpfun(SEIRlctCoxian3.ode)
SEIRlctCoxian4.ode <- compiler::cmpfun(SEIRlctCoxian4.ode)
SEIRlctCoxian5.ode <- compiler::cmpfun(SEIRlctCoxian5.ode)
SEIRlctCoxian6.ode <- compiler::cmpfun(SEIRlctCoxian6.ode)
SEIRlctCoxian7.ode <- compiler::cmpfun(SEIRlctCoxian7.ode)
SEIRlctCoxian8.ode <- compiler::cmpfun(SEIRlctCoxian8.ode)
SEIRlctCoxian9.ode <- compiler::cmpfun(SEIRlctCoxian9.ode)
SEIRlctCoxian10.ode <- compiler::cmpfun(SEIRlctCoxian10.ode)
SEIRlctCoxian11.ode <- compiler::cmpfun(SEIRlctCoxian11.ode)
SEIRlctCoxian12.ode <- compiler::cmpfun(SEIRlctCoxian12.ode)
SEIRlctCoxian13.ode <- compiler::cmpfun(SEIRlctCoxian13.ode)
SEIRlctCoxian14.ode <- compiler::cmpfun(SEIRlctCoxian14.ode)
SEIRlctCoxian15.ode <- compiler::cmpfun(SEIRlctCoxian15.ode)
SEIRlctCoxian16.ode <- compiler::cmpfun(SEIRlctCoxian16.ode)
SEIRlctCoxian17.ode <- compiler::cmpfun(SEIRlctCoxian17.ode)
SEIRlctCoxian18.ode <- compiler::cmpfun(SEIRlctCoxian18.ode)
SEIRlctCoxian19.ode <- compiler::cmpfun(SEIRlctCoxian19.ode)
SEIRlctCoxian20.ode <- compiler::cmpfun(SEIRlctCoxian20.ode)
SEIRlctCoxian21.ode <- compiler::cmpfun(SEIRlctCoxian21.ode)
SEIRlctCoxian22.ode <- compiler::cmpfun(SEIRlctCoxian22.ode)
SEIRlctCoxian23.ode <- compiler::cmpfun(SEIRlctCoxian23.ode)
SEIRlctCoxian24.ode <- compiler::cmpfun(SEIRlctCoxian24.ode)
SEIRlctCoxian25.ode <- compiler::cmpfun(SEIRlctCoxian25.ode)

mthd = "ode45"
atol= 1e-6
IC=c(S=0.9999, E=0.0001, I=0, R=0) # for SEIR.ode()

parms1 <- parms; parms1['kE']=1; parms1['kI']=1; parms1; SEIRpt.init(parms1)
b1=benchmark(SEIR.1 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.1={ode(y=ICs, times=tms, func=SEIRlctCoxian1.ode, parms = parms1, method = mthd, atol=atol)},
             SEIRpt.1 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms1, method = mthd, atol=atol)}, 
             replications = reps)

parms2 <- parms; parms2['kE']=2; parms2['kI']=2; parms2; SEIRpt.init(parms2)
b2=benchmark(SEIR.2 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
          SEIRlctCoxian.2={ode(y=ICs, times=tms, func=SEIRlctCoxian2.ode, parms = parms2, method = mthd, atol=atol)},
          SEIRpt.2 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms2, method = mthd, atol=atol)}, 
          replications = reps)

parms3 <- parms; parms3['kE']=3; parms3['kI']=3; parms3; SEIRpt.init(parms3)
b3=benchmark(SEIR.3 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
          SEIRlctCoxian.3={ode(y=ICs, times=tms, func=SEIRlctCoxian3.ode, parms = parms3, method = mthd, atol=atol)},
          SEIRpt.3 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms3, method = mthd, atol=atol)}, 
          replications = reps)

parms4 <- parms; parms4['kE']=4; parms4['kI']=4; parms4; SEIRpt.init(parms4)
b4=benchmark(SEIR.4 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
          SEIRlctCoxian.4={ode(y=ICs, times=tms, func=SEIRlctCoxian4.ode, parms = parms4, method = mthd, atol=atol)},
          SEIRpt.4 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms4, method = mthd, atol=atol)}, 
          replications = reps)

parms5 <- parms; parms5['kE']=5; parms5['kI']=5; parms5; SEIRpt.init(parms5)
b5=benchmark(SEIR.5 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.5={ode(y=ICs, times=tms, func=SEIRlctCoxian5.ode, parms = parms5, method = mthd, atol=atol)},
             SEIRpt.5 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms5, method = mthd, atol=atol)}, 
             replications = reps)

parms6 <- parms; parms6['kE']=6; parms6['kI']=6; parms6; SEIRpt.init(parms6)
b6=benchmark(SEIR.6 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.6={ode(y=ICs, times=tms, func=SEIRlctCoxian6.ode, parms = parms6, method = mthd, atol=atol)},
             SEIRpt.6 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms6, method = mthd, atol=atol)}, 
             replications = reps)

parms7 <- parms; parms7['kE']=7; parms7['kI']=7; parms7; SEIRpt.init(parms7)
b7=benchmark(SEIR.7 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.7={ode(y=ICs, times=tms, func=SEIRlctCoxian7.ode, parms = parms7, method = mthd, atol=atol)},
             SEIRpt.7 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms7, method = mthd, atol=atol)}, 
             replications = reps)

parms8 <- parms; parms8['kE']=8; parms8['kI']=8; parms8; SEIRpt.init(parms8)
b8=benchmark(SEIR.8 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.8={ode(y=ICs, times=tms, func=SEIRlctCoxian8.ode, parms = parms8, method = mthd, atol=atol)},
             SEIRpt.8 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms8, method = mthd, atol=atol)}, 
             replications = reps)
             
parms9 <- parms; parms9['kE']=9; parms9['kI']=9; parms9; SEIRpt.init(parms9)
b9=benchmark(SEIR.9 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.9={ode(y=ICs, times=tms, func=SEIRlctCoxian9.ode, parms = parms9, method = mthd, atol=atol)},
             SEIRpt.9 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms9, method = mthd, atol=atol)}, 
             replications = reps)
             
parms10 <- parms; parms10['kE']=10; parms10['kI']=10; parms10; SEIRpt.init(parms10)
b10=benchmark(SEIR.10 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.10={ode(y=ICs, times=tms, func=SEIRlctCoxian10.ode, parms = parms10, method = mthd, atol=atol)},
             SEIRpt.10 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms10, method = mthd, atol=atol)}, 
             replications = reps)
   
parms11 <- parms; parms11['kE']=11; parms11['kI']=11; parms11; SEIRpt.init(parms11)
b11=benchmark(SEIR.11 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.11={ode(y=ICs, times=tms, func=SEIRlctCoxian11.ode, parms = parms11, method = mthd, atol=atol)},
             SEIRpt.11 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms11, method = mthd, atol=atol)}, 
             replications = reps)

parms12 <- parms; parms12['kE']=12; parms12['kI']=12; parms12; SEIRpt.init(parms12)
b12=benchmark(SEIR.12 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.12={ode(y=ICs, times=tms, func=SEIRlctCoxian12.ode, parms = parms12, method = mthd, atol=atol)},
             SEIRpt.12 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms12, method = mthd, atol=atol)}, 
             replications = reps)

parms13 <- parms; parms13['kE']=13; parms13['kI']=13; parms13; SEIRpt.init(parms13)
b13=benchmark(SEIR.13 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.13={ode(y=ICs, times=tms, func=SEIRlctCoxian13.ode, parms = parms13, method = mthd, atol=atol)},
             SEIRpt.13 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms13, method = mthd, atol=atol)}, 
             replications = reps)
      
parms14 <- parms; parms14['kE']=14; parms14['kI']=14; parms14; SEIRpt.init(parms14)
b14=benchmark(SEIR.14 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.14={ode(y=ICs, times=tms, func=SEIRlctCoxian14.ode, parms = parms14, method = mthd, atol=atol)},
             SEIRpt.14 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms14, method = mthd, atol=atol)}, 
             replications = reps)

parms15 <- parms; parms15['kE']=15; parms15['kI']=15; parms15; SEIRpt.init(parms15)
b15=benchmark(SEIR.15 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.15={ode(y=ICs, times=tms, func=SEIRlctCoxian15.ode, parms = parms15, method = mthd, atol=atol)},
             SEIRpt.15 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms15, method = mthd, atol=atol)}, 
             replications = reps)

parms16 <- parms; parms16['kE']=16; parms16['kI']=16; parms16; SEIRpt.init(parms16)
b16=benchmark(SEIR.16 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.16={ode(y=ICs, times=tms, func=SEIRlctCoxian16.ode, parms = parms16, method = mthd, atol=atol)},
             SEIRpt.16 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms16, method = mthd, atol=atol)}, 
             replications = reps)
             
parms17 <- parms; parms17['kE']=17; parms17['kI']=17; parms17; SEIRpt.init(parms17)
b17=benchmark(SEIR.17 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.17={ode(y=ICs, times=tms, func=SEIRlctCoxian17.ode, parms = parms17, method = mthd, atol=atol)},
             SEIRpt.17 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms17, method = mthd, atol=atol)}, 
             replications = reps)

parms18 <- parms; parms18['kE']=18; parms18['kI']=18; parms18; SEIRpt.init(parms18)
b18=benchmark(SEIR.18 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.18={ode(y=ICs, times=tms, func=SEIRlctCoxian18.ode, parms = parms18, method = mthd, atol=atol)},
             SEIRpt.18 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms18, method = mthd, atol=atol)}, 
             replications = reps)

parms19 <- parms; parms19['kE']=19; parms19['kI']=19; parms19; SEIRpt.init(parms19)
b19=benchmark(SEIR.19 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.19={ode(y=ICs, times=tms, func=SEIRlctCoxian19.ode, parms = parms19, method = mthd, atol=atol)},
             SEIRpt.19 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms19, method = mthd, atol=atol)}, 
             replications = reps)

parms20 <- parms; parms20['kE']=20; parms20['kI']=20; parms20; SEIRpt.init(parms20)
b20=benchmark(SEIR.20 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.20={ode(y=ICs, times=tms, func=SEIRlctCoxian20.ode, parms = parms20, method = mthd, atol=atol)},
             SEIRpt.20 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms20, method = mthd, atol=atol)}, 
             replications = reps)

parms21 <- parms; parms21['kE']=21; parms21['kI']=21; parms21; SEIRpt.init(parms21)
b21=benchmark(SEIR.21 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.21={ode(y=ICs, times=tms, func=SEIRlctCoxian21.ode, parms = parms21, method = mthd, atol=atol)},
             SEIRpt.21 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms21, method = mthd, atol=atol)}, 
             replications = reps)

parms22 <- parms; parms22['kE']=22; parms22['kI']=22; parms22; SEIRpt.init(parms22)
b22=benchmark(SEIR.22 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.22={ode(y=ICs, times=tms, func=SEIRlctCoxian22.ode, parms = parms22, method = mthd, atol=atol)},
             SEIRpt.22 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms22, method = mthd, atol=atol)}, 
             replications = reps)

parms23 <- parms; parms23['kE']=23; parms23['kI']=23; parms23; SEIRpt.init(parms23)
b23=benchmark(SEIR.23 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.23={ode(y=ICs, times=tms, func=SEIRlctCoxian23.ode, parms = parms23, method = mthd, atol=atol)},
             SEIRpt.23 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms23, method = mthd, atol=atol)}, 
             replications = reps)

parms24 <- parms; parms24['kE']=24; parms24['kI']=24; parms24; SEIRpt.init(parms24)
b24=benchmark(SEIR.24 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.24={ode(y=ICs, times=tms, func=SEIRlctCoxian24.ode, parms = parms24, method = mthd, atol=atol)},
             SEIRpt.24 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms24, method = mthd, atol=atol)}, 
             replications = reps)
             
parms25 <- parms; parms25['kE']=25; parms25['kI']=25; parms25; SEIRpt.init(parms25)
b25=benchmark(SEIR.25 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlctCoxian.25={ode(y=ICs, times=tms, func=SEIRlctCoxian25.ode, parms = parms25, method = mthd, atol=atol)},
             SEIRpt.25 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms25, method = mthd, atol=atol)}, 
             replications = reps)
                                                                                                       
b1; b2; b3; b4; b5; b6; b7; b8; b9; b10; b11; b12; b13; b14; b15; b16; b17; b18; b19; b20; b21; b22; b23; b24; b25


# First, reorganize the data and change up some labeling...
# 1. Rearrange rows to the right order, save a new copy...
out <- rbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25)

out$N <- stringr::str_split_fixed(out$test,'\\.',2)[,2]
out$Model <- stringr::str_split_fixed(out$test,'\\.',2)[,1]
# out$Model <- gsub("SEIR$","SEIR (Exponential)",out$Model)
# out$Model <- gsub("lct"," (Coxian / LCT)",out$Model)
# out$Model <- gsub("pt"," (Coxian as Phase-Type / GLCT)",out$Model)
out$Model <- gsub("SEIR$","Standard: Exponential latent & infectious periods (N independent)",out$Model)
out$Model <- gsub("lctCoxian","LCT:    Coxian latent & infectious periods",out$Model)
out$Model <- gsub("pt","GLCT: Coxian latent & infectious periods (phase-type formulation)",out$Model)
out$Model <- gsub("SEIR",'',out$Model)
out$Model <- factor(out$Model, levels = unique(out$Model) )
out$N <- as.numeric(out$N)
head(out)

# 2. Plot using "N" to group the trios...

#x11(12,6); 
pdf("SEIRCoxianbenchmark.pdf",12,6);
ggplot(out, aes(x=factor(N), y=user.self, fill=Model, group=test)) + 
  geom_bar(stat="identity", position="dodge", width=0.6) + theme_bw(base_size=18) + 
  xlab("N (Substates of E and I; Model dimension = 2 + 2N)") + ylab("Elapsed Time (seconds)") + 
  scale_fill_grey(start=.75, end=0.2, name="SEIR Model (Coxian)")+
  theme(legend.position=c(.296,.85), 
        legend.background = element_rect(fill="white", size = 0.5, color = "darkgray")) 
dev.off()

# Compare the trajectories for the two approaches to implementing the LCT
outLCT  = ode(y=ICs, times=tms, func=SEIRlctCoxian25.ode, parms = parms25, method = mthd, atol=atol)
outGLCT = ode(y=ICs, times=tms, func=SEIRpt.ode,    parms = parms25, method = mthd, atol=atol) 

pdf("SEIRCoxian-compare.pdf",12,6);
matplot(tms, cbind(S=outLCT[,2],E=colSums(t(outLCT[,3:27])),I=colSums(t(outLCT[,28:52])), R=outLCT[,53]),
        type="l",lty=1, lwd=2, col=c("green","lightblue","red","gray"),
        xlab="Time", ylab="% Population")
matplot(tms, cbind(S=outGLCT[,2],E=colSums(t(outGLCT[,3:27])),I=colSums(t(outGLCT[,28:52])), R=outGLCT[,53]),
        type="l",lty=3, lwd=2, col="black",add=TRUE,
        xlab="Time", ylab="% Population")
dev.off()


