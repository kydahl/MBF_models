using Symbolics
using LinearAlgebra

function matrix_simplify(mat)
    # go through entry by entry and simplify
    out_mat = mat
    for i in 1:size(mat, 1), j in 1:size(mat, 2)
        out_mat[i,j] = simplify(mat[i,j]; expand = true)
    end
    return out_mat
end


@variables lQ lL lP lG sigma pQ pL pP pG

A = [
    -lQ+(1-pQ)*lQ       pQ*lQ               0   0   
    (1-sigma)*(1-pL)*lL -lL+sigma*(1-pL)*lL pL*lL 0 
    (1-sigma)*(1-pP)*lP sigma*(1-pP)*lP     -lP   pP*lP
    (1-sigma)*(1-pG)*lG sigma*(1-pG)*lG     0   -lG
]

spec_mat = matrix_simplify(-[1;0;0;0]*transpose(A*[1; 1; 1; 1]))


@variables mu eta gV gR

GV = gV/(mu+gV)
GV_E = gV/(mu+gV+eta)
GR = gR/(mu+gR)
GR_E = gR/(mu+gR+eta)

W = matrix_simplify((mu+eta)*I - transpose(A) - GV_E * GR_E * spec_mat)
GammaE = matrix_simplify(W)
Z = matrix_simplify(mu*I - transpose(A) - GV * GR * spec_mat)
GammaI = matrix_simplify(Z)
tau_pre = eta*I + GV*((1-(eta/(mu+gR+eta)))*(eta/(mu+gV+eta)) + (eta/(mu+gR+eta)))*GR*spec_mat
tauE = matrix_simplify(tau_pre * GammaE)

temp_inv = matrix_simplify(GammaI*tauE)

temp_RVH = simplify(temp_inv[3,4]; expand = true)
using SymPy
temp_RVH2 = SymPy.simplify(temp_RVH)