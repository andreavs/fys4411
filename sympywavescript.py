from sympy import *
import numpy as np
import numpy.linalg as npl

def gradient(phi):
	grad = [diff(phi,x), diff(phi,y), diff(phi,z)]
	return grad

def cleanup(phi):
	return phi.factor().subs(r,R)

def laplace(phi):
	return diff(diff(phi,x),x) + diff(diff(phi,y),y) + diff(diff(phi,z),z)

if __name__ == '__main__':
	x,y,z,r,a = symbols('x y z r a')
	phi1s, phi2s = symbols('phi1s phi2s')
	r = sqrt(x*x + y*y + z*z)
	R = Symbol('r')
	phi1s = exp(-a*r)
	s1grad = gradient(phi1s)

	phi2s = (1-a*r/2)*exp(-a*r/2)
	s2grad = gradient(phi2s)
	s2lapl = laplace(phi2s)
	#print (cleanup(s2grad))
	phi2px = Symbol('phi2px')
	phi2px = a*x*exp(-a*r/2)
	px2grad = gradient(phi2px)

	print cleanup(s1grad[0])
	print (s2lapl.factor().subs(r,R))
	# print cleanup(px2grad[0])
	# print cleanup(px2grad[1])
	# print cleanup(px2grad[2])
	# print cleanup(laplace(phi2px))
	Rx, Ry, Rz = symbols('Rx Ry Rz')
	Rn = Symbol('Rn')
	Rn = sqrt(Rx*Rx + Ry*Ry + Rz*Rz)
	Rnn = Symbol('Rn')
	r1p1, r1p2 = symbols('r1p1 r1p2')
	r1p1 = sqrt((x+Rx/2)**2 + (y+Ry/2)**2 + (z+Rz/2)**2);
	r1p2 = sqrt((x-Rx/2)**2 + (y-Ry/2)**2 + (z-Rz/2)**2);
	r1p1n = Symbol('r1p1')
	r1p2n = Symbol('r1p2')
	psi = Symbol('psi');
	psi = exp(-a*r1p1)# + exp(-a*r1p2)
	print laplace(psi).factor().subs(r,R).subs(Rn, Rnn).subs(r1p1, r1p1n).subs(r1p2, r1p2n)

	# x1,x2,y1,y2,z1,z2,r1,r2 = symbols('')