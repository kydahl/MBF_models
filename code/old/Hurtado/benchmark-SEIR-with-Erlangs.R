####################################################################################
## Benchmarking SEIR model numerical solutions
## ----------------------------------------------
##  Case 1: Exponential dwell time distributions (dimension = 3 if R omitted)
##  Case 2: Classic LCT (Erlang dwell times), hard-coded number of substates
##  Case 3: SEIR w/ phase-type distributions, parameterized w/ Erlang distributions
####################################################################################
#
#   NOTE: This R script was run on a 64-bit Windows 10 machine using R v3.6.3.
#         Running this on other operating systems may require small modifications
#         like replacing the win.graph() calls with x11() or quartz() calls.
#

library(deSolve)

## Parameterization....
reps<-500
parms = c(b=1, muE=4, muI=7, cvE=1, cvI=1, kE=1, kI=1) 


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

# SEIR with hard-coded Erlang dwell times
SEIRlct1.ode <- function(tm,z,ps) {
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
SEIRlct2.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE =2
  kI =2
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot  - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dI1 =  z[3]*kE/muE  - z[4]*kI/muI
  dI2 =  z[4]*kI/muI  - z[5]*kI/muI
  dR  =  z[5]*kI/muI
  
  return(list(c(dS, dE1, dE2, dI1, dI2, dR)))
}
SEIRlct3.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 3
  kI = 3
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dI1 =  z[4]*kE/muE  - z[5]*kI/muI
  dI2 =  z[5]*kI/muI  - z[6]*kI/muI
  dI3 =  z[6]*kI/muI  - z[7]*kI/muI
  dR  =  z[7]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dI1, dI2, dI3, dR)))
}
SEIRlct4.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 4
  kI = 4
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dI1 =  z[5]*kE/muE  - z[6]*kI/muI
  dI2 =  z[6]*kI/muI  - z[7]*kI/muI
  dI3 =  z[7]*kI/muI  - z[8]*kI/muI
  dI4 =  z[8]*kI/muI  - z[9]*kI/muI
  dR  =  z[9]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dI1, dI2, dI3, dI4, dR)))
}
SEIRlct5.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 5
  kI = 5
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dI1 =  z[6]*kE/muE  - z[7]*kI/muI
  dI2 =  z[7]*kI/muI  - z[8]*kI/muI
  dI3 =  z[8]*kI/muI  - z[9]*kI/muI
  dI4 =  z[9]*kI/muI  - z[10]*kI/muI
  dI5 =  z[10]*kI/muI  - z[11]*kI/muI
  dR  =  z[11]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dI1, dI2, dI3, dI4, dI5, dR)))
}
SEIRlct6.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 6
  kI = 6
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dI1 =  z[7]*kE/muE  - z[8]*kI/muI
  dI2 =  z[8]*kI/muI  - z[9]*kI/muI
  dI3 =  z[9]*kI/muI  - z[10]*kI/muI
  dI4 =  z[10]*kI/muI  - z[11]*kI/muI
  dI5 =  z[11]*kI/muI  - z[12]*kI/muI
  dI6 =  z[12]*kI/muI  - z[13]*kI/muI
  dR  =  z[13]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dI1, dI2, dI3, dI4, dI5, dI6, dR)))
}
SEIRlct7.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 7
  kI = 7
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dI1 =  z[8]*kE/muE  - z[9]*kI/muI
  dI2 =  z[9]*kI/muI  - z[10]*kI/muI
  dI3 =  z[10]*kI/muI  - z[11]*kI/muI
  dI4 =  z[11]*kI/muI  - z[12]*kI/muI
  dI5 =  z[12]*kI/muI  - z[13]*kI/muI
  dI6 =  z[13]*kI/muI  - z[14]*kI/muI
  dI7 =  z[14]*kI/muI  - z[15]*kI/muI
  dR  =  z[15]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dR)))
}
SEIRlct8.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 8
  kI = 8
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dI1 =  z[9]*kE/muE  - z[10]*kI/muI
  dI2 =  z[10]*kI/muI  - z[11]*kI/muI
  dI3 =  z[11]*kI/muI  - z[12]*kI/muI
  dI4 =  z[12]*kI/muI  - z[13]*kI/muI
  dI5 =  z[13]*kI/muI  - z[14]*kI/muI
  dI6 =  z[14]*kI/muI  - z[15]*kI/muI
  dI7 =  z[15]*kI/muI  - z[16]*kI/muI
  dI8 =  z[16]*kI/muI  - z[17]*kI/muI
  dR  =  z[17]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dR)))
}
SEIRlct9.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 9
  kI = 9
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dI1 =  z[10]*kE/muE  - z[11]*kI/muI
  dI2 =  z[11]*kI/muI  - z[12]*kI/muI
  dI3 =  z[12]*kI/muI  - z[13]*kI/muI
  dI4 =  z[13]*kI/muI  - z[14]*kI/muI
  dI5 =  z[14]*kI/muI  - z[15]*kI/muI
  dI6 =  z[15]*kI/muI  - z[16]*kI/muI
  dI7 =  z[16]*kI/muI  - z[17]*kI/muI
  dI8 =  z[17]*kI/muI  - z[18]*kI/muI
  dI9 =  z[18]*kI/muI  - z[19]*kI/muI
  dR  =  z[19]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dR)))
}

SEIRlct10.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 10
  kI = 10
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dI1 =  z[11]*kE/muE  - z[12]*kI/muI
  dI2 =  z[12]*kI/muI  - z[13]*kI/muI
  dI3 =  z[13]*kI/muI  - z[14]*kI/muI
  dI4 =  z[14]*kI/muI  - z[15]*kI/muI
  dI5 =  z[15]*kI/muI  - z[16]*kI/muI
  dI6 =  z[16]*kI/muI  - z[17]*kI/muI
  dI7 =  z[17]*kI/muI  - z[18]*kI/muI
  dI8 =  z[18]*kI/muI  - z[19]*kI/muI
  dI9 =  z[19]*kI/muI  - z[20]*kI/muI
  dI10 =  z[20]*kI/muI  - z[21]*kI/muI
  dR  =  z[21]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dR)))
}

SEIRlct11.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 11
  kI = 11
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dI1 =  z[12]*kE/muE  - z[13]*kI/muI
  dI2 =  z[13]*kI/muI  - z[14]*kI/muI
  dI3 =  z[14]*kI/muI  - z[15]*kI/muI
  dI4 =  z[15]*kI/muI  - z[16]*kI/muI
  dI5 =  z[16]*kI/muI  - z[17]*kI/muI
  dI6 =  z[17]*kI/muI  - z[18]*kI/muI
  dI7 =  z[18]*kI/muI  - z[19]*kI/muI
  dI8 =  z[19]*kI/muI  - z[20]*kI/muI
  dI9 =  z[20]*kI/muI  - z[21]*kI/muI
  dI10 =  z[21]*kI/muI  - z[22]*kI/muI
  dI11 =  z[22]*kI/muI  - z[23]*kI/muI
  dR  =  z[23]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dR)))
}

SEIRlct12.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 12
  kI = 12
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dI1 =  z[13]*kE/muE  - z[14]*kI/muI
  dI2 =  z[14]*kI/muI  - z[15]*kI/muI
  dI3 =  z[15]*kI/muI  - z[16]*kI/muI
  dI4 =  z[16]*kI/muI  - z[17]*kI/muI
  dI5 =  z[17]*kI/muI  - z[18]*kI/muI
  dI6 =  z[18]*kI/muI  - z[19]*kI/muI
  dI7 =  z[19]*kI/muI  - z[20]*kI/muI
  dI8 =  z[20]*kI/muI  - z[21]*kI/muI
  dI9 =  z[21]*kI/muI  - z[22]*kI/muI
  dI10 =  z[22]*kI/muI  - z[23]*kI/muI
  dI11 =  z[23]*kI/muI  - z[24]*kI/muI
  dI12 =  z[24]*kI/muI  - z[25]*kI/muI
  dR  =  z[25]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dR)))
}

SEIRlct13.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 13
  kI = 13
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dI1 =  z[14]*kE/muE  - z[15]*kI/muI
  dI2 =  z[15]*kI/muI  - z[16]*kI/muI
  dI3 =  z[16]*kI/muI  - z[17]*kI/muI
  dI4 =  z[17]*kI/muI  - z[18]*kI/muI
  dI5 =  z[18]*kI/muI  - z[19]*kI/muI
  dI6 =  z[19]*kI/muI  - z[20]*kI/muI
  dI7 =  z[20]*kI/muI  - z[21]*kI/muI
  dI8 =  z[21]*kI/muI  - z[22]*kI/muI
  dI9 =  z[22]*kI/muI  - z[23]*kI/muI
  dI10 =  z[23]*kI/muI  - z[24]*kI/muI
  dI11 =  z[24]*kI/muI  - z[25]*kI/muI
  dI12 =  z[25]*kI/muI  - z[26]*kI/muI
  dI13 =  z[26]*kI/muI  - z[27]*kI/muI
  dR  =  z[27]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dR)))
}

SEIRlct14.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 14
  kI = 14
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dE14 =  z[14]*kE/muE  - z[15]*kE/muE
  dI1 =  z[15]*kE/muE  - z[16]*kI/muI
  dI2 =  z[16]*kI/muI  - z[17]*kI/muI
  dI3 =  z[17]*kI/muI  - z[18]*kI/muI
  dI4 =  z[18]*kI/muI  - z[19]*kI/muI
  dI5 =  z[19]*kI/muI  - z[20]*kI/muI
  dI6 =  z[20]*kI/muI  - z[21]*kI/muI
  dI7 =  z[21]*kI/muI  - z[22]*kI/muI
  dI8 =  z[22]*kI/muI  - z[23]*kI/muI
  dI9 =  z[23]*kI/muI  - z[24]*kI/muI
  dI10 =  z[24]*kI/muI  - z[25]*kI/muI
  dI11 =  z[25]*kI/muI  - z[26]*kI/muI
  dI12 =  z[26]*kI/muI  - z[27]*kI/muI
  dI13 =  z[27]*kI/muI  - z[28]*kI/muI
  dI14 =  z[28]*kI/muI  - z[29]*kI/muI
  dR  =  z[29]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dR)))
}

SEIRlct15.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 15
  kI = 15
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dE14 =  z[14]*kE/muE  - z[15]*kE/muE
  dE15 =  z[15]*kE/muE  - z[16]*kE/muE
  dI1 =  z[16]*kE/muE  - z[17]*kI/muI
  dI2 =  z[17]*kI/muI  - z[18]*kI/muI
  dI3 =  z[18]*kI/muI  - z[19]*kI/muI
  dI4 =  z[19]*kI/muI  - z[20]*kI/muI
  dI5 =  z[20]*kI/muI  - z[21]*kI/muI
  dI6 =  z[21]*kI/muI  - z[22]*kI/muI
  dI7 =  z[22]*kI/muI  - z[23]*kI/muI
  dI8 =  z[23]*kI/muI  - z[24]*kI/muI
  dI9 =  z[24]*kI/muI  - z[25]*kI/muI
  dI10 =  z[25]*kI/muI  - z[26]*kI/muI
  dI11 =  z[26]*kI/muI  - z[27]*kI/muI
  dI12 =  z[27]*kI/muI  - z[28]*kI/muI
  dI13 =  z[28]*kI/muI  - z[29]*kI/muI
  dI14 =  z[29]*kI/muI  - z[30]*kI/muI
  dI15 =  z[30]*kI/muI  - z[31]*kI/muI
  dR  =  z[31]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dR)))
}

SEIRlct16.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 16
  kI = 16
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dE14 =  z[14]*kE/muE  - z[15]*kE/muE
  dE15 =  z[15]*kE/muE  - z[16]*kE/muE
  dE16 =  z[16]*kE/muE  - z[17]*kE/muE
  dI1 =  z[17]*kE/muE  - z[18]*kI/muI
  dI2 =  z[18]*kI/muI  - z[19]*kI/muI
  dI3 =  z[19]*kI/muI  - z[20]*kI/muI
  dI4 =  z[20]*kI/muI  - z[21]*kI/muI
  dI5 =  z[21]*kI/muI  - z[22]*kI/muI
  dI6 =  z[22]*kI/muI  - z[23]*kI/muI
  dI7 =  z[23]*kI/muI  - z[24]*kI/muI
  dI8 =  z[24]*kI/muI  - z[25]*kI/muI
  dI9 =  z[25]*kI/muI  - z[26]*kI/muI
  dI10 =  z[26]*kI/muI  - z[27]*kI/muI
  dI11 =  z[27]*kI/muI  - z[28]*kI/muI
  dI12 =  z[28]*kI/muI  - z[29]*kI/muI
  dI13 =  z[29]*kI/muI  - z[30]*kI/muI
  dI14 =  z[30]*kI/muI  - z[31]*kI/muI
  dI15 =  z[31]*kI/muI  - z[32]*kI/muI
  dI16 =  z[32]*kI/muI  - z[33]*kI/muI
  dR  =  z[33]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dR)))
}

SEIRlct17.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 17
  kI = 17
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dE14 =  z[14]*kE/muE  - z[15]*kE/muE
  dE15 =  z[15]*kE/muE  - z[16]*kE/muE
  dE16 =  z[16]*kE/muE  - z[17]*kE/muE
  dE17 =  z[17]*kE/muE  - z[18]*kE/muE
  dI1 =  z[18]*kE/muE  - z[19]*kI/muI
  dI2 =  z[19]*kI/muI  - z[20]*kI/muI
  dI3 =  z[20]*kI/muI  - z[21]*kI/muI
  dI4 =  z[21]*kI/muI  - z[22]*kI/muI
  dI5 =  z[22]*kI/muI  - z[23]*kI/muI
  dI6 =  z[23]*kI/muI  - z[24]*kI/muI
  dI7 =  z[24]*kI/muI  - z[25]*kI/muI
  dI8 =  z[25]*kI/muI  - z[26]*kI/muI
  dI9 =  z[26]*kI/muI  - z[27]*kI/muI
  dI10 =  z[27]*kI/muI  - z[28]*kI/muI
  dI11 =  z[28]*kI/muI  - z[29]*kI/muI
  dI12 =  z[29]*kI/muI  - z[30]*kI/muI
  dI13 =  z[30]*kI/muI  - z[31]*kI/muI
  dI14 =  z[31]*kI/muI  - z[32]*kI/muI
  dI15 =  z[32]*kI/muI  - z[33]*kI/muI
  dI16 =  z[33]*kI/muI  - z[34]*kI/muI
  dI17 =  z[34]*kI/muI  - z[35]*kI/muI
  dR  =  z[35]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dR)))
}

SEIRlct18.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 18
  kI = 18
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dE14 =  z[14]*kE/muE  - z[15]*kE/muE
  dE15 =  z[15]*kE/muE  - z[16]*kE/muE
  dE16 =  z[16]*kE/muE  - z[17]*kE/muE
  dE17 =  z[17]*kE/muE  - z[18]*kE/muE
  dE18 =  z[18]*kE/muE  - z[19]*kE/muE
  dI1 =  z[19]*kE/muE  - z[20]*kI/muI
  dI2 =  z[20]*kI/muI  - z[21]*kI/muI
  dI3 =  z[21]*kI/muI  - z[22]*kI/muI
  dI4 =  z[22]*kI/muI  - z[23]*kI/muI
  dI5 =  z[23]*kI/muI  - z[24]*kI/muI
  dI6 =  z[24]*kI/muI  - z[25]*kI/muI
  dI7 =  z[25]*kI/muI  - z[26]*kI/muI
  dI8 =  z[26]*kI/muI  - z[27]*kI/muI
  dI9 =  z[27]*kI/muI  - z[28]*kI/muI
  dI10 =  z[28]*kI/muI  - z[29]*kI/muI
  dI11 =  z[29]*kI/muI  - z[30]*kI/muI
  dI12 =  z[30]*kI/muI  - z[31]*kI/muI
  dI13 =  z[31]*kI/muI  - z[32]*kI/muI
  dI14 =  z[32]*kI/muI  - z[33]*kI/muI
  dI15 =  z[33]*kI/muI  - z[34]*kI/muI
  dI16 =  z[34]*kI/muI  - z[35]*kI/muI
  dI17 =  z[35]*kI/muI  - z[36]*kI/muI
  dI18 =  z[36]*kI/muI  - z[37]*kI/muI
  dR  =  z[37]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dR)))
}

SEIRlct19.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 19
  kI = 19
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dE14 =  z[14]*kE/muE  - z[15]*kE/muE
  dE15 =  z[15]*kE/muE  - z[16]*kE/muE
  dE16 =  z[16]*kE/muE  - z[17]*kE/muE
  dE17 =  z[17]*kE/muE  - z[18]*kE/muE
  dE18 =  z[18]*kE/muE  - z[19]*kE/muE
  dE19 =  z[19]*kE/muE  - z[20]*kE/muE
  dI1 =  z[20]*kE/muE  - z[21]*kI/muI
  dI2 =  z[21]*kI/muI  - z[22]*kI/muI
  dI3 =  z[22]*kI/muI  - z[23]*kI/muI
  dI4 =  z[23]*kI/muI  - z[24]*kI/muI
  dI5 =  z[24]*kI/muI  - z[25]*kI/muI
  dI6 =  z[25]*kI/muI  - z[26]*kI/muI
  dI7 =  z[26]*kI/muI  - z[27]*kI/muI
  dI8 =  z[27]*kI/muI  - z[28]*kI/muI
  dI9 =  z[28]*kI/muI  - z[29]*kI/muI
  dI10 =  z[29]*kI/muI  - z[30]*kI/muI
  dI11 =  z[30]*kI/muI  - z[31]*kI/muI
  dI12 =  z[31]*kI/muI  - z[32]*kI/muI
  dI13 =  z[32]*kI/muI  - z[33]*kI/muI
  dI14 =  z[33]*kI/muI  - z[34]*kI/muI
  dI15 =  z[34]*kI/muI  - z[35]*kI/muI
  dI16 =  z[35]*kI/muI  - z[36]*kI/muI
  dI17 =  z[36]*kI/muI  - z[37]*kI/muI
  dI18 =  z[37]*kI/muI  - z[38]*kI/muI
  dI19 =  z[38]*kI/muI  - z[39]*kI/muI
  dR  =  z[39]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dR)))
}

SEIRlct20.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 20
  kI = 20
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dE14 =  z[14]*kE/muE  - z[15]*kE/muE
  dE15 =  z[15]*kE/muE  - z[16]*kE/muE
  dE16 =  z[16]*kE/muE  - z[17]*kE/muE
  dE17 =  z[17]*kE/muE  - z[18]*kE/muE
  dE18 =  z[18]*kE/muE  - z[19]*kE/muE
  dE19 =  z[19]*kE/muE  - z[20]*kE/muE
  dE20 =  z[20]*kE/muE  - z[21]*kE/muE
  dI1 =  z[21]*kE/muE  - z[22]*kI/muI
  dI2 =  z[22]*kI/muI  - z[23]*kI/muI
  dI3 =  z[23]*kI/muI  - z[24]*kI/muI
  dI4 =  z[24]*kI/muI  - z[25]*kI/muI
  dI5 =  z[25]*kI/muI  - z[26]*kI/muI
  dI6 =  z[26]*kI/muI  - z[27]*kI/muI
  dI7 =  z[27]*kI/muI  - z[28]*kI/muI
  dI8 =  z[28]*kI/muI  - z[29]*kI/muI
  dI9 =  z[29]*kI/muI  - z[30]*kI/muI
  dI10 =  z[30]*kI/muI  - z[31]*kI/muI
  dI11 =  z[31]*kI/muI  - z[32]*kI/muI
  dI12 =  z[32]*kI/muI  - z[33]*kI/muI
  dI13 =  z[33]*kI/muI  - z[34]*kI/muI
  dI14 =  z[34]*kI/muI  - z[35]*kI/muI
  dI15 =  z[35]*kI/muI  - z[36]*kI/muI
  dI16 =  z[36]*kI/muI  - z[37]*kI/muI
  dI17 =  z[37]*kI/muI  - z[38]*kI/muI
  dI18 =  z[38]*kI/muI  - z[39]*kI/muI
  dI19 =  z[39]*kI/muI  - z[40]*kI/muI
  dI20 =  z[40]*kI/muI  - z[41]*kI/muI
  dR  =  z[41]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dE20, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI20, dR)))
}

SEIRlct21.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 21
  kI = 21
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dE14 =  z[14]*kE/muE  - z[15]*kE/muE
  dE15 =  z[15]*kE/muE  - z[16]*kE/muE
  dE16 =  z[16]*kE/muE  - z[17]*kE/muE
  dE17 =  z[17]*kE/muE  - z[18]*kE/muE
  dE18 =  z[18]*kE/muE  - z[19]*kE/muE
  dE19 =  z[19]*kE/muE  - z[20]*kE/muE
  dE20 =  z[20]*kE/muE  - z[21]*kE/muE
  dE21 =  z[21]*kE/muE  - z[22]*kE/muE
  dI1 =  z[22]*kE/muE  - z[23]*kI/muI
  dI2 =  z[23]*kI/muI  - z[24]*kI/muI
  dI3 =  z[24]*kI/muI  - z[25]*kI/muI
  dI4 =  z[25]*kI/muI  - z[26]*kI/muI
  dI5 =  z[26]*kI/muI  - z[27]*kI/muI
  dI6 =  z[27]*kI/muI  - z[28]*kI/muI
  dI7 =  z[28]*kI/muI  - z[29]*kI/muI
  dI8 =  z[29]*kI/muI  - z[30]*kI/muI
  dI9 =  z[30]*kI/muI  - z[31]*kI/muI
  dI10 =  z[31]*kI/muI  - z[32]*kI/muI
  dI11 =  z[32]*kI/muI  - z[33]*kI/muI
  dI12 =  z[33]*kI/muI  - z[34]*kI/muI
  dI13 =  z[34]*kI/muI  - z[35]*kI/muI
  dI14 =  z[35]*kI/muI  - z[36]*kI/muI
  dI15 =  z[36]*kI/muI  - z[37]*kI/muI
  dI16 =  z[37]*kI/muI  - z[38]*kI/muI
  dI17 =  z[38]*kI/muI  - z[39]*kI/muI
  dI18 =  z[39]*kI/muI  - z[40]*kI/muI
  dI19 =  z[40]*kI/muI  - z[41]*kI/muI
  dI20 =  z[41]*kI/muI  - z[42]*kI/muI
  dI21 =  z[42]*kI/muI  - z[43]*kI/muI
  dR  =  z[43]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dE20, dE21, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI20, dI21, dR)))
}

SEIRlct22.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 22
  kI = 22
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dE14 =  z[14]*kE/muE  - z[15]*kE/muE
  dE15 =  z[15]*kE/muE  - z[16]*kE/muE
  dE16 =  z[16]*kE/muE  - z[17]*kE/muE
  dE17 =  z[17]*kE/muE  - z[18]*kE/muE
  dE18 =  z[18]*kE/muE  - z[19]*kE/muE
  dE19 =  z[19]*kE/muE  - z[20]*kE/muE
  dE20 =  z[20]*kE/muE  - z[21]*kE/muE
  dE21 =  z[21]*kE/muE  - z[22]*kE/muE
  dE22 =  z[22]*kE/muE  - z[23]*kE/muE
  dI1 =  z[23]*kE/muE  - z[24]*kI/muI
  dI2 =  z[24]*kI/muI  - z[25]*kI/muI
  dI3 =  z[25]*kI/muI  - z[26]*kI/muI
  dI4 =  z[26]*kI/muI  - z[27]*kI/muI
  dI5 =  z[27]*kI/muI  - z[28]*kI/muI
  dI6 =  z[28]*kI/muI  - z[29]*kI/muI
  dI7 =  z[29]*kI/muI  - z[30]*kI/muI
  dI8 =  z[30]*kI/muI  - z[31]*kI/muI
  dI9 =  z[31]*kI/muI  - z[32]*kI/muI
  dI10 =  z[32]*kI/muI  - z[33]*kI/muI
  dI11 =  z[33]*kI/muI  - z[34]*kI/muI
  dI12 =  z[34]*kI/muI  - z[35]*kI/muI
  dI13 =  z[35]*kI/muI  - z[36]*kI/muI
  dI14 =  z[36]*kI/muI  - z[37]*kI/muI
  dI15 =  z[37]*kI/muI  - z[38]*kI/muI
  dI16 =  z[38]*kI/muI  - z[39]*kI/muI
  dI17 =  z[39]*kI/muI  - z[40]*kI/muI
  dI18 =  z[40]*kI/muI  - z[41]*kI/muI
  dI19 =  z[41]*kI/muI  - z[42]*kI/muI
  dI20 =  z[42]*kI/muI  - z[43]*kI/muI
  dI21 =  z[43]*kI/muI  - z[44]*kI/muI
  dI22 =  z[44]*kI/muI  - z[45]*kI/muI
  dR  =  z[45]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dE20, dE21, dE22, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI20, dI21, dI22, dR)))
}

SEIRlct23.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 23
  kI = 23
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dE14 =  z[14]*kE/muE  - z[15]*kE/muE
  dE15 =  z[15]*kE/muE  - z[16]*kE/muE
  dE16 =  z[16]*kE/muE  - z[17]*kE/muE
  dE17 =  z[17]*kE/muE  - z[18]*kE/muE
  dE18 =  z[18]*kE/muE  - z[19]*kE/muE
  dE19 =  z[19]*kE/muE  - z[20]*kE/muE
  dE20 =  z[20]*kE/muE  - z[21]*kE/muE
  dE21 =  z[21]*kE/muE  - z[22]*kE/muE
  dE22 =  z[22]*kE/muE  - z[23]*kE/muE
  dE23 =  z[23]*kE/muE  - z[24]*kE/muE
  dI1 =  z[24]*kE/muE  - z[25]*kI/muI
  dI2 =  z[25]*kI/muI  - z[26]*kI/muI
  dI3 =  z[26]*kI/muI  - z[27]*kI/muI
  dI4 =  z[27]*kI/muI  - z[28]*kI/muI
  dI5 =  z[28]*kI/muI  - z[29]*kI/muI
  dI6 =  z[29]*kI/muI  - z[30]*kI/muI
  dI7 =  z[30]*kI/muI  - z[31]*kI/muI
  dI8 =  z[31]*kI/muI  - z[32]*kI/muI
  dI9 =  z[32]*kI/muI  - z[33]*kI/muI
  dI10 =  z[33]*kI/muI  - z[34]*kI/muI
  dI11 =  z[34]*kI/muI  - z[35]*kI/muI
  dI12 =  z[35]*kI/muI  - z[36]*kI/muI
  dI13 =  z[36]*kI/muI  - z[37]*kI/muI
  dI14 =  z[37]*kI/muI  - z[38]*kI/muI
  dI15 =  z[38]*kI/muI  - z[39]*kI/muI
  dI16 =  z[39]*kI/muI  - z[40]*kI/muI
  dI17 =  z[40]*kI/muI  - z[41]*kI/muI
  dI18 =  z[41]*kI/muI  - z[42]*kI/muI
  dI19 =  z[42]*kI/muI  - z[43]*kI/muI
  dI20 =  z[43]*kI/muI  - z[44]*kI/muI
  dI21 =  z[44]*kI/muI  - z[45]*kI/muI
  dI22 =  z[45]*kI/muI  - z[46]*kI/muI
  dI23 =  z[46]*kI/muI  - z[47]*kI/muI
  dR  =  z[47]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dE20, dE21, dE22, dE23, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI20, dI21, dI22, dI23, dR)))
}

SEIRlct24.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 24
  kI = 24
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dE14 =  z[14]*kE/muE  - z[15]*kE/muE
  dE15 =  z[15]*kE/muE  - z[16]*kE/muE
  dE16 =  z[16]*kE/muE  - z[17]*kE/muE
  dE17 =  z[17]*kE/muE  - z[18]*kE/muE
  dE18 =  z[18]*kE/muE  - z[19]*kE/muE
  dE19 =  z[19]*kE/muE  - z[20]*kE/muE
  dE20 =  z[20]*kE/muE  - z[21]*kE/muE
  dE21 =  z[21]*kE/muE  - z[22]*kE/muE
  dE22 =  z[22]*kE/muE  - z[23]*kE/muE
  dE23 =  z[23]*kE/muE  - z[24]*kE/muE
  dE24 =  z[24]*kE/muE  - z[25]*kE/muE
  dI1 =  z[25]*kE/muE  - z[26]*kI/muI
  dI2 =  z[26]*kI/muI  - z[27]*kI/muI
  dI3 =  z[27]*kI/muI  - z[28]*kI/muI
  dI4 =  z[28]*kI/muI  - z[29]*kI/muI
  dI5 =  z[29]*kI/muI  - z[30]*kI/muI
  dI6 =  z[30]*kI/muI  - z[31]*kI/muI
  dI7 =  z[31]*kI/muI  - z[32]*kI/muI
  dI8 =  z[32]*kI/muI  - z[33]*kI/muI
  dI9 =  z[33]*kI/muI  - z[34]*kI/muI
  dI10 =  z[34]*kI/muI  - z[35]*kI/muI
  dI11 =  z[35]*kI/muI  - z[36]*kI/muI
  dI12 =  z[36]*kI/muI  - z[37]*kI/muI
  dI13 =  z[37]*kI/muI  - z[38]*kI/muI
  dI14 =  z[38]*kI/muI  - z[39]*kI/muI
  dI15 =  z[39]*kI/muI  - z[40]*kI/muI
  dI16 =  z[40]*kI/muI  - z[41]*kI/muI
  dI17 =  z[41]*kI/muI  - z[42]*kI/muI
  dI18 =  z[42]*kI/muI  - z[43]*kI/muI
  dI19 =  z[43]*kI/muI  - z[44]*kI/muI
  dI20 =  z[44]*kI/muI  - z[45]*kI/muI
  dI21 =  z[45]*kI/muI  - z[46]*kI/muI
  dI22 =  z[46]*kI/muI  - z[47]*kI/muI
  dI23 =  z[47]*kI/muI  - z[48]*kI/muI
  dI24 =  z[48]*kI/muI  - z[49]*kI/muI
  dR  =  z[49]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dE20, dE21, dE22, dE23, dE24, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI20, dI21, dI22, dI23, dI24, dR)))
}

SEIRlct25.ode <- function(tm,z,ps) {
  b=ps[["b"]]     # beta
  muE=ps[["muE"]] # mean time spent in E
  muI=ps[["muI"]] # mean time in I
  kE = 25
  kI = 25
  
  Itot = sum(z[1+kE+(1:kI)]) # sum(Ivec) 
  
  dS  = -b*z[1]*Itot
  dE1 =  b*z[1]*Itot - z[2]*kE/muE
  dE2 =  z[2]*kE/muE  - z[3]*kE/muE
  dE3 =  z[3]*kE/muE  - z[4]*kE/muE
  dE4 =  z[4]*kE/muE  - z[5]*kE/muE
  dE5 =  z[5]*kE/muE  - z[6]*kE/muE
  dE6 =  z[6]*kE/muE  - z[7]*kE/muE
  dE7 =  z[7]*kE/muE  - z[8]*kE/muE
  dE8 =  z[8]*kE/muE  - z[9]*kE/muE
  dE9 =  z[9]*kE/muE  - z[10]*kE/muE
  dE10 =  z[10]*kE/muE  - z[11]*kE/muE
  dE11 =  z[11]*kE/muE  - z[12]*kE/muE
  dE12 =  z[12]*kE/muE  - z[13]*kE/muE
  dE13 =  z[13]*kE/muE  - z[14]*kE/muE
  dE14 =  z[14]*kE/muE  - z[15]*kE/muE
  dE15 =  z[15]*kE/muE  - z[16]*kE/muE
  dE16 =  z[16]*kE/muE  - z[17]*kE/muE
  dE17 =  z[17]*kE/muE  - z[18]*kE/muE
  dE18 =  z[18]*kE/muE  - z[19]*kE/muE
  dE19 =  z[19]*kE/muE  - z[20]*kE/muE
  dE20 =  z[20]*kE/muE  - z[21]*kE/muE
  dE21 =  z[21]*kE/muE  - z[22]*kE/muE
  dE22 =  z[22]*kE/muE  - z[23]*kE/muE
  dE23 =  z[23]*kE/muE  - z[24]*kE/muE
  dE24 =  z[24]*kE/muE  - z[25]*kE/muE
  dE25 =  z[25]*kE/muE  - z[26]*kE/muE
  dI1 =  z[26]*kE/muE  - z[27]*kI/muI
  dI2 =  z[27]*kI/muI  - z[28]*kI/muI
  dI3 =  z[28]*kI/muI  - z[29]*kI/muI
  dI4 =  z[29]*kI/muI  - z[30]*kI/muI
  dI5 =  z[30]*kI/muI  - z[31]*kI/muI
  dI6 =  z[31]*kI/muI  - z[32]*kI/muI
  dI7 =  z[32]*kI/muI  - z[33]*kI/muI
  dI8 =  z[33]*kI/muI  - z[34]*kI/muI
  dI9 =  z[34]*kI/muI  - z[35]*kI/muI
  dI10 =  z[35]*kI/muI  - z[36]*kI/muI
  dI11 =  z[36]*kI/muI  - z[37]*kI/muI
  dI12 =  z[37]*kI/muI  - z[38]*kI/muI
  dI13 =  z[38]*kI/muI  - z[39]*kI/muI
  dI14 =  z[39]*kI/muI  - z[40]*kI/muI
  dI15 =  z[40]*kI/muI  - z[41]*kI/muI
  dI16 =  z[41]*kI/muI  - z[42]*kI/muI
  dI17 =  z[42]*kI/muI  - z[43]*kI/muI
  dI18 =  z[43]*kI/muI  - z[44]*kI/muI
  dI19 =  z[44]*kI/muI  - z[45]*kI/muI
  dI20 =  z[45]*kI/muI  - z[46]*kI/muI
  dI21 =  z[46]*kI/muI  - z[47]*kI/muI
  dI22 =  z[47]*kI/muI  - z[48]*kI/muI
  dI23 =  z[48]*kI/muI  - z[49]*kI/muI
  dI24 =  z[49]*kI/muI  - z[50]*kI/muI
  dI25 =  z[50]*kI/muI  - z[51]*kI/muI
  dR  =  z[51]*kI/muI
  
  return(list(c(dS, dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8, dE9, dE10, dE11, dE12, dE13, dE14, dE15, dE16, dE17, dE18, dE19, dE20, dE21, dE22, dE23, dE24, dE25, dI1, dI2, dI3, dI4, dI5, dI6, dI7, dI8, dI9, dI10, dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18, dI19, dI20, dI21, dI22, dI23, dI24, dI25, dR)))
}

# SEIR with Erlang dwell times in a GLCT framework
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
  
  # These are Erlang distributions framed in a Phase-type distribution context,
  # where vector a = (1 0 ... 0) and matrix A are as follows...
  
  aE = matrix(0,nrow = kE, ncol=1); aE[1] = 1;
  AE = kE/muE*(diag(rep(-1,kE),kE)); if(kE>1) for(i in 1:(kE-1)) {AE[i,i+1] = kE/muE}
  
  aI = matrix(0,nrow = kI, ncol=1); aI[1] = 1;
  AI = kI/muI*(diag(rep(-1,kI),kI)); if(kI>1) for(i in 1:(kI-1)) {AI[i,i+1] = kI/muI}
  
  
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
SEIRlct1.ode <- compiler::cmpfun(SEIRlct1.ode)
SEIRlct2.ode <- compiler::cmpfun(SEIRlct2.ode)
SEIRlct3.ode <- compiler::cmpfun(SEIRlct3.ode)
SEIRlct4.ode <- compiler::cmpfun(SEIRlct4.ode)
SEIRlct5.ode <- compiler::cmpfun(SEIRlct5.ode)
SEIRlct6.ode <- compiler::cmpfun(SEIRlct6.ode)
SEIRlct7.ode <- compiler::cmpfun(SEIRlct7.ode)
SEIRlct8.ode <- compiler::cmpfun(SEIRlct8.ode)
SEIRlct9.ode <- compiler::cmpfun(SEIRlct9.ode)
SEIRlct10.ode <- compiler::cmpfun(SEIRlct10.ode)
SEIRlct11.ode <- compiler::cmpfun(SEIRlct11.ode)
SEIRlct12.ode <- compiler::cmpfun(SEIRlct12.ode)
SEIRlct13.ode <- compiler::cmpfun(SEIRlct13.ode)
SEIRlct14.ode <- compiler::cmpfun(SEIRlct14.ode)
SEIRlct15.ode <- compiler::cmpfun(SEIRlct15.ode)
SEIRlct16.ode <- compiler::cmpfun(SEIRlct16.ode)
SEIRlct17.ode <- compiler::cmpfun(SEIRlct17.ode)
SEIRlct18.ode <- compiler::cmpfun(SEIRlct18.ode)
SEIRlct19.ode <- compiler::cmpfun(SEIRlct19.ode)
SEIRlct20.ode <- compiler::cmpfun(SEIRlct20.ode)
SEIRlct21.ode <- compiler::cmpfun(SEIRlct21.ode)
SEIRlct22.ode <- compiler::cmpfun(SEIRlct22.ode)
SEIRlct23.ode <- compiler::cmpfun(SEIRlct23.ode)
SEIRlct24.ode <- compiler::cmpfun(SEIRlct24.ode)
SEIRlct25.ode <- compiler::cmpfun(SEIRlct25.ode)

mthd = "ode45"
atol= 1e-6
IC=c(S=0.9999, E=0.0001, I=0, R=0) # for SEIR.ode()

parms1 <- parms; parms1['kE']=1; parms1['kI']=1; parms1; SEIRpt.init(parms1)
b1=benchmark(SEIR.1 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.1={ode(y=ICs, times=tms, func=SEIRlct1.ode, parms = parms1, method = mthd, atol=atol)},
             SEIRpt.1 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms1, method = mthd, atol=atol)}, 
             replications = reps)

parms2 <- parms; parms2['kE']=2; parms2['kI']=2; parms2; SEIRpt.init(parms2)
b2=benchmark(SEIR.2 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
          SEIRlct.2={ode(y=ICs, times=tms, func=SEIRlct2.ode, parms = parms2, method = mthd, atol=atol)},
          SEIRpt.2 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms2, method = mthd, atol=atol)}, 
          replications = reps)

parms3 <- parms; parms3['kE']=3; parms3['kI']=3; parms3; SEIRpt.init(parms3)
b3=benchmark(SEIR.3 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
          SEIRlct.3={ode(y=ICs, times=tms, func=SEIRlct3.ode, parms = parms3, method = mthd, atol=atol)},
          SEIRpt.3 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms3, method = mthd, atol=atol)}, 
          replications = reps)

parms4 <- parms; parms4['kE']=4; parms4['kI']=4; parms4; SEIRpt.init(parms4)
b4=benchmark(SEIR.4 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
          SEIRlct.4={ode(y=ICs, times=tms, func=SEIRlct4.ode, parms = parms4, method = mthd, atol=atol)},
          SEIRpt.4 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms4, method = mthd, atol=atol)}, 
          replications = reps)

parms5 <- parms; parms5['kE']=5; parms5['kI']=5; parms5; SEIRpt.init(parms5)
b5=benchmark(SEIR.5 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.5={ode(y=ICs, times=tms, func=SEIRlct5.ode, parms = parms5, method = mthd, atol=atol)},
             SEIRpt.5 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms5, method = mthd, atol=atol)}, 
             replications = reps)

parms6 <- parms; parms6['kE']=6; parms6['kI']=6; parms6; SEIRpt.init(parms6)
b6=benchmark(SEIR.6 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.6={ode(y=ICs, times=tms, func=SEIRlct6.ode, parms = parms6, method = mthd, atol=atol)},
             SEIRpt.6 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms6, method = mthd, atol=atol)}, 
             replications = reps)

parms7 <- parms; parms7['kE']=7; parms7['kI']=7; parms7; SEIRpt.init(parms7)
b7=benchmark(SEIR.7 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.7={ode(y=ICs, times=tms, func=SEIRlct7.ode, parms = parms7, method = mthd, atol=atol)},
             SEIRpt.7 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms7, method = mthd, atol=atol)}, 
             replications = reps)

parms8 <- parms; parms8['kE']=8; parms8['kI']=8; parms8; SEIRpt.init(parms8)
b8=benchmark(SEIR.8 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.8={ode(y=ICs, times=tms, func=SEIRlct8.ode, parms = parms8, method = mthd, atol=atol)},
             SEIRpt.8 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms8, method = mthd, atol=atol)}, 
             replications = reps)
             
parms9 <- parms; parms9['kE']=9; parms9['kI']=9; parms9; SEIRpt.init(parms9)
b9=benchmark(SEIR.9 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.9={ode(y=ICs, times=tms, func=SEIRlct9.ode, parms = parms9, method = mthd, atol=atol)},
             SEIRpt.9 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms9, method = mthd, atol=atol)}, 
             replications = reps)
             
parms10 <- parms; parms10['kE']=10; parms10['kI']=10; parms10; SEIRpt.init(parms10)
b10=benchmark(SEIR.10 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.10={ode(y=ICs, times=tms, func=SEIRlct10.ode, parms = parms10, method = mthd, atol=atol)},
             SEIRpt.10 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms10, method = mthd, atol=atol)}, 
             replications = reps)
   
parms11 <- parms; parms11['kE']=11; parms11['kI']=11; parms11; SEIRpt.init(parms11)
b11=benchmark(SEIR.11 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.11={ode(y=ICs, times=tms, func=SEIRlct11.ode, parms = parms11, method = mthd, atol=atol)},
             SEIRpt.11 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms11, method = mthd, atol=atol)}, 
             replications = reps)

parms12 <- parms; parms12['kE']=12; parms12['kI']=12; parms12; SEIRpt.init(parms12)
b12=benchmark(SEIR.12 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.12={ode(y=ICs, times=tms, func=SEIRlct12.ode, parms = parms12, method = mthd, atol=atol)},
             SEIRpt.12 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms12, method = mthd, atol=atol)}, 
             replications = reps)

parms13 <- parms; parms13['kE']=13; parms13['kI']=13; parms13; SEIRpt.init(parms13)
b13=benchmark(SEIR.13 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.13={ode(y=ICs, times=tms, func=SEIRlct13.ode, parms = parms13, method = mthd, atol=atol)},
             SEIRpt.13 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms13, method = mthd, atol=atol)}, 
             replications = reps)
      
parms14 <- parms; parms14['kE']=14; parms14['kI']=14; parms14; SEIRpt.init(parms14)
b14=benchmark(SEIR.14 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.14={ode(y=ICs, times=tms, func=SEIRlct14.ode, parms = parms14, method = mthd, atol=atol)},
             SEIRpt.14 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms14, method = mthd, atol=atol)}, 
             replications = reps)

parms15 <- parms; parms15['kE']=15; parms15['kI']=15; parms15; SEIRpt.init(parms15)
b15=benchmark(SEIR.15 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.15={ode(y=ICs, times=tms, func=SEIRlct15.ode, parms = parms15, method = mthd, atol=atol)},
             SEIRpt.15 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms15, method = mthd, atol=atol)}, 
             replications = reps)

parms16 <- parms; parms16['kE']=16; parms16['kI']=16; parms16; SEIRpt.init(parms16)
b16=benchmark(SEIR.16 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.16={ode(y=ICs, times=tms, func=SEIRlct16.ode, parms = parms16, method = mthd, atol=atol)},
             SEIRpt.16 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms16, method = mthd, atol=atol)}, 
             replications = reps)
             
parms17 <- parms; parms17['kE']=17; parms17['kI']=17; parms17; SEIRpt.init(parms17)
b17=benchmark(SEIR.17 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.17={ode(y=ICs, times=tms, func=SEIRlct17.ode, parms = parms17, method = mthd, atol=atol)},
             SEIRpt.17 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms17, method = mthd, atol=atol)}, 
             replications = reps)

parms18 <- parms; parms18['kE']=18; parms18['kI']=18; parms18; SEIRpt.init(parms18)
b18=benchmark(SEIR.18 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.18={ode(y=ICs, times=tms, func=SEIRlct18.ode, parms = parms18, method = mthd, atol=atol)},
             SEIRpt.18 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms18, method = mthd, atol=atol)}, 
             replications = reps)

parms19 <- parms; parms19['kE']=19; parms19['kI']=19; parms19; SEIRpt.init(parms19)
b19=benchmark(SEIR.19 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.19={ode(y=ICs, times=tms, func=SEIRlct19.ode, parms = parms19, method = mthd, atol=atol)},
             SEIRpt.19 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms19, method = mthd, atol=atol)}, 
             replications = reps)

parms20 <- parms; parms20['kE']=20; parms20['kI']=20; parms20; SEIRpt.init(parms20)
b20=benchmark(SEIR.20 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.20={ode(y=ICs, times=tms, func=SEIRlct20.ode, parms = parms20, method = mthd, atol=atol)},
             SEIRpt.20 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms20, method = mthd, atol=atol)}, 
             replications = reps)

parms21 <- parms; parms21['kE']=21; parms21['kI']=21; parms21; SEIRpt.init(parms21)
b21=benchmark(SEIR.21 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.21={ode(y=ICs, times=tms, func=SEIRlct21.ode, parms = parms21, method = mthd, atol=atol)},
             SEIRpt.21 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms21, method = mthd, atol=atol)}, 
             replications = reps)

parms22 <- parms; parms22['kE']=22; parms22['kI']=22; parms22; SEIRpt.init(parms22)
b22=benchmark(SEIR.22 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.22={ode(y=ICs, times=tms, func=SEIRlct22.ode, parms = parms22, method = mthd, atol=atol)},
             SEIRpt.22 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms22, method = mthd, atol=atol)}, 
             replications = reps)

parms23 <- parms; parms23['kE']=23; parms23['kI']=23; parms23; SEIRpt.init(parms23)
b23=benchmark(SEIR.23 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.23={ode(y=ICs, times=tms, func=SEIRlct23.ode, parms = parms23, method = mthd, atol=atol)},
             SEIRpt.23 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms23, method = mthd, atol=atol)}, 
             replications = reps)

parms24 <- parms; parms24['kE']=24; parms24['kI']=24; parms24; SEIRpt.init(parms24)
b24=benchmark(SEIR.24 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.24={ode(y=ICs, times=tms, func=SEIRlct24.ode, parms = parms24, method = mthd, atol=atol)},
             SEIRpt.24 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms24, method = mthd, atol=atol)}, 
             replications = reps)
             
parms25 <- parms; parms25['kE']=25; parms25['kI']=25; parms25; SEIRpt.init(parms25)
b25=benchmark(SEIR.25 ={ode(y=IC, times = tms, func = SEIR.ode, parms = parms, method = mthd, atol=atol)}, 
             SEIRlct.25={ode(y=ICs, times=tms, func=SEIRlct25.ode, parms = parms25, method = mthd, atol=atol)},
             SEIRpt.25 ={ode(y=ICs, times=tms, func=SEIRpt.ode,   parms = parms25, method = mthd, atol=atol)}, 
             replications = reps)
                                                                                                       
b1; b2; b3; b4; b5; b6; b7; b8; b9; b10; b11; b12; b13; b14; b15; b16; b17; b18; b19; b20; b21; b22; b23; b24; b25


# First, reorganize the data and change up some labeling...
# 1. Rearrange rows to the right order, save a new copy...
out <- rbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25)

out$N <- stringr::str_split_fixed(out$test,'\\.',2)[,2]
out$Model <- stringr::str_split_fixed(out$test,'\\.',2)[,1]
# out$Model <- gsub("SEIR$","SEIR (Exponential)",out$Model)
# out$Model <- gsub("lct"," (Erlang / LCT)",out$Model)
# out$Model <- gsub("pt"," (Erlang as Phase-Type / GLCT)",out$Model)
out$Model <- gsub("SEIR$","Standard: Exponential latent & infectious periods (N independent)",out$Model)
out$Model <- gsub("lct","LCT:    Erlang latent & infectious periods",out$Model)
out$Model <- gsub("pt","GLCT: Erlang latent & infectious periods (phase-type formulation)",out$Model)
out$Model <- gsub("SEIR",'',out$Model)
out$Model <- factor(out$Model, levels = unique(out$Model) )
out$N <- as.numeric(out$N)
head(out)

# 2. Plot using "N" to group the trios...

#win.graph(12,6); 
pdf("SEIRbenchmark.pdf",12,6);
ggplot(out, aes(x=factor(N), y=user.self, fill=Model, group=test)) + 
  geom_bar(stat="identity", position="dodge", width=0.6) + theme_bw(base_size=18) + 
  xlab("N (Substates of E and I; Model dimension = 2 + 2N)") + ylab("Elapsed Time (seconds)") + 
  scale_fill_grey(start=.75, end=0.2, name="SEIR Model (Erlang)")+
  theme(legend.position=c(.296,.85), 
        legend.background = element_rect(fill="white", size = 0.5, color = "darkgray")) 
dev.off()

# Compare the trajectories for the two approaches to implementing the LCT
outLCT  = ode(y=ICs, times=tms, func=SEIRlct25.ode, parms = parms25, method = mthd, atol=atol)
outGLCT = ode(y=ICs, times=tms, func=SEIRpt.ode,    parms = parms25, method = mthd, atol=atol) 

pdf("SEIR-compare.pdf",12,6);
matplot(tms, cbind(S=outLCT[,2],E=colSums(t(outLCT[,3:27])),I=colSums(t(outLCT[,28:52])), R=outLCT[,53]),
        type="l",lty=1, lwd=2, col="black", #c("green","blue","red","gray"),
        xlab="Time", ylab="% Population")
matplot(tms, cbind(S=outGLCT[,2],E=colSums(t(outGLCT[,3:27])),I=colSums(t(outGLCT[,28:52])), R=outGLCT[,53]),
        type="l",lty=3, lwd=2, col="gray",add=TRUE,
        xlab="Time", ylab="% Population")
dev.off()


