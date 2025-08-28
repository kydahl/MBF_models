using Symbolics
using LinearAlgebra

@variables lQ lL lP lG sigma pQ pL pP pG

A = [
    -lQ+(1-pQ)*lQ       pQ*lQ               0f0   0f0   
    (1-sigma)*(1-pL)*lL -lL+sigma*(1-pL)*lL pL*lL 0f0 
    (1-sigma)*(1-pP)*lP sigma*(1-pP)*lP     -lP   pP*lP
    (1-sigma)*(1-pG)*lG sigma*(1-pG)*lG     0f0   -lG
]

@variables mu eta GV gR gV GR

W = (mu+eta)*I - transpose(A)
Z = mu*I - transpose(A)

Y = [-eta 0    0    -GV*((1-gR)*gV + gR)*GR*pG*lG;
     0    -eta 0    0;
     0    0    -eta 0;
     0    0    0    -eta]

temp_inv = inv(Z)*(-Y)*inv(W)

temp_RVH = temp_inv[3,4]

temp_RVH2 = simplify(temp_RVH;expand = true, threaded = true)