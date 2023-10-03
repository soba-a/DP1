import sympy as sp
import numpy as np
import math as m

warm_z, Fbx, Fpx, theta, warm_fz, factM = sp.symbols('warm_z, Fbx, Fpx, theta, warm_fz, factM')
wabre_fz = 2

r_loc = sp.Matrix([0, 0, 0])
fact_loc = sp.Matrix([205, 0, warm_z])
Fb_loc = sp.Matrix([410, 65, warm_z+25])
Fp_loc = sp.Matrix([410, 100, warm_z+25])
warm_loc = sp.Matrix([205, 0, warm_z])
wabre_loc = sp.Matrix([410, 0, warm_z+25])

# r = sp.Matrix([rx, ry, rz])
fact = sp.Matrix([-factM * sp.cos(theta), 0, factM * sp.sin(theta)])
Fb = sp.Matrix([Fbx, 0, 0])
Fp = sp.Matrix([Fpx, 0, 0])
warm = sp.Matrix([0, 0, -warm_fz])
wbre = sp.Matrix([0, 0, -wabre_fz])
r =  -(fact + Fb + Fp + warm + wbre)
ract = -(r + fact + Fb + Fp + warm + wbre)

Mact = fact_loc.cross(fact)
# print(Mact)
MFb = Fb_loc.cross(Fb)
# print(MFb)
MFp = Fp_loc.cross(Fp)
# print(MFp)
Mwarm = warm_loc.cross(warm)
# print(Mwarm)
Mwabre = wabre_loc.cross(wbre)
equations = [ract, Mact, MFb, MFp, Mwarm, Mwabre]
print(sp.solve([ract, Mact, MFb, MFp, Mwarm, Mwabre], (Fbx, Fpx, factM, warm_fz)))