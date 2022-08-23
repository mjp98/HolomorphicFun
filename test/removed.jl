using HolomorphicFun

s = Sqrt{false}(1,0.5)
e  = Exponential(1.)
S = SFun(s)
a = e*S

regular(a,1+10im)
removedsingularity(singularities(a)[1],1+10im)



a = SFun(Sqrt{false}(1.0+0im,-0.5))
b = SFun(ScalarPole(-1.0,2))
c = a*b
s = singularities(c)
