import numpy as np, matplotlib.pyplot as plt

def getLatticeSpacing(beta):
	# if beta < 5.7: raise ValueError("Beta should be larger than 5.7!")
	r = 0.5
	bval = (beta - 6)
	a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*0.5
	return a

print getLatticeSpacing(6.0)

X = np.linspace(5.7,6.5,1000)
Y = getLatticeSpacing(X)

plt.plot(X,Y,label=r"$a = \frac{1}{2}\exp(-1.6805 - 1.7139\beta + 0.8155\beta^2 - 0.6667\beta^3)$")
plt.xlabel(r"$\beta$")
plt.ylabel(r"$a[fm]$")
plt.title("Beta value vs lattice spacing")
plt.legend()
plt.show()