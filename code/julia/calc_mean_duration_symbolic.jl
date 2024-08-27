using Symbolics

@variables f mu lQ lL lP lG pQ pL pP pG

# Note to self: this level of detail is not necessary. 
# We can just use the mean blood-feeding stage duration neglecting mortality.
# i.e. mu = 0
# For the parameter values I've been using (mu = 1/(20 days)), including mu amounts to ~60 minute reduction in the duration when it is set to 1 day


mat = [mu+lQ -lQ 0 0
       -f*(1-pL)*lL mu+lL-(1-f)*(1-pL)*lL -pL*lL 0
       -f*(1-pP)*lP -(1-f)*(1-pP)*lP mu+lP -pP*lP
       -f*(1-pG)*lG -(1-f)*(1-pG)*lG 0 mu+lG]

neg_mat = [-(mu+lQ) lQ 0 0
f*(1-pL)*lL -(mu+lL)+(1-f)*(1-pL)*lL pL*lL 0
f*(1-pP)*lP (1-f)*(1-pP)*lP -(mu+lP) pP*lP
f*(1-pG)*lG (1-f)*(1-pG)*lG 0 -(mu+lG)]

temp_inv = simplify(inv(neg_mat))

alpha = transpose([1 0 0 0])
one_vec = transpose([1 1 1 1])

pre_vec = transpose(alpha) * temp_inv

simple_pre_vec = pre_vec
simple_pre_vec[1] = simplify(pre_vec[1])
simple_pre_vec[2] = simplify(pre_vec[2])
simple_pre_vec[3] = simplify(pre_vec[3])
simple_pre_vec[4] = simplify(pre_vec[4])


mean_dur = simple_pre_vec * one_vec


mean_dur2 = transpose(one_vec) * transpose(temp_inv) * alpha