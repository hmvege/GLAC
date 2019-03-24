import sympy as sp


# Setting general matrix A
u,w,xi0 = sp.symbols("u,w,xi0",real=True)

h0 = (u**2 - w**2)*sp.exp(2*sp.I*u) + sp.exp(-sp.I*u)*(8*u**2*sp.cos(w) + 2*sp.I*u*(3*u**2 + w**2)*xi0)
h1 = 2*u*sp.exp(2*sp.I*u) - sp.exp(-sp.I*u)*(2*u*sp.cos(w) - sp.I*(3*u**2 - w**2)*xi0)
h2 = sp.exp(2*sp.I*u) - sp.exp(-sp.I*u)*(sp.cos(w) + 3*sp.I*u*xi0)

print "\nh0: "
sp.pprint(h0)
print "\nReal part of h0: "
sp.pprint(sp.re(h0))
print "\nImag part of h0: "
sp.pprint(sp.im(h0))

print "\nh1: "
sp.pprint(h1)
print "\nReal part of h1: "
sp.pprint(sp.re(h1))
print "\nImag part of h1: "
sp.pprint(sp.im(h1))

print "\nh2: "
sp.pprint(h2)
print "\nReal part of h2: "
sp.pprint(sp.re(h2))
print "\nImag part of h2: "
sp.pprint(sp.im(h2))
