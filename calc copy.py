import sympy as sp
import numpy as np
import math as m

rx, ry, rz, rax, ray, raz, fact, theta, Fbx, Fbz, warm, Fp =  sp.symbols('rx, ry, rz, rax, ray, raz, fact, theta, Fbx, warm, Fbz, Fp')
# assume lpiv is in line with lblade
wbre = 2 * 9.81
Fresz = 40
warm_z = -25

# Solve for torque
torqueBlade = 1 * 9550 / 4600 * 1000 # n *mm

# through blade fbd
Fresx = -torqueBlade / 101.6
print('Fresx =',Fresx)

# through blade pulley fbd
T1 = 45
T2 = (torqueBlade + 38.1 * T1) / 38.1 
print('T1 =',T1)
print('T2 =',T2)
        


# assume 
Rpiv_loc = sp.Matrix([0, 0, 0]).T
Ract_loc = sp.Matrix([205, 0, warm_z]).T
Fact_loc = sp.Matrix([205, 0, warm_z]).T
Fres_loc = sp.Matrix([410, 65, warm_z + 25 - 101.6]).T
Fpul_loc = sp.Matrix([410, 100, warm_z + 25]).T
Warm_loc = sp.Matrix([205, 0, warm_z]).T
Wbre_loc = sp.Matrix([410, 0, warm_z + 25]).T

# rpiv_loc2 = sp.Matrix([-410, 0, 0]).T
# ract_loc2 = sp.Matrix([-205, 0, warm_z]).T
# fact_loc2 = sp.Matrix([-205, 0, warm_z]).T
# Fb_loc2 = sp.Matrix([0, 65, warm_z + 25]).T
# Fp_loc2 = sp.Matrix([0, 100, warm_z + 25]).T
# warm_loc2 = sp.Matrix([-205, 0, warm_z]).T
# wbre_loc2 = sp.Matrix([0, 0, warm_z + 25]).T

RpivV = sp.Matrix([rx, ry, rz]).T
RactV = sp.Matrix([rax, ray, raz]).T
FactV = sp.Matrix([fact * sp.cos(theta), 0, -fact * sp.sin(theta)]).T
FresV = sp.Matrix([Fresx, 0, Fresz]).T
FpulV = sp.Matrix([-Fp * np.cos(np.radians(3.545)), 0, Fp * np.sin(np.radians(3.545))]).T
WarmV = sp.Matrix([0, 0, -warm]).T
WbreV = sp.Matrix([0, 0, -wbre]).T

Mract = Ract_loc.cross(RactV)
Mrpiv = Rpiv_loc.cross(RpivV)
Mact = Fact_loc.cross(FactV)
MFres = Fres_loc.cross(FresV)
MFpul = Fpul_loc.cross(FpulV)
Mwarm = Warm_loc.cross(WarmV)
Mwabre = Wbre_loc.cross(WbreV)

# Mract2 = ract_loc2.cross(RactV)
# Mrpiv2 = rpiv_loc2.cross(rpivV)
# Mact2 = fact_loc2.cross(FactV)
# MFb2 = Fb_loc2.cross(FresV)
# MFp2 = Fp_loc2.cross(FpulV)
# Mwarm2 = warm_loc2.cross(WarmV)
# Mwabre2 = wbre_loc2.cross(WbreV)

Ftot = RpivV + FactV + FresV + WarmV + WbreV
Mtot = Mrpiv + Mact + MFres + Mwarm + Mwabre
# Mtot2 = Mract2 + Mrpiv2 + MFp2 + Mwarm2 + Mwabre2 + Mact2
# Ftot = rpivV + FactV + FpulV + WarmV + WbreV
# Mtot = Mrpiv + Mact + MFp + Mwarm + Mwabre


equations = Ftot.row_insert(1, Mtot)
# print(Ftot.col(0))
# print(Ftot.col(1))
# print(Ftot.col(2))
# print(Mtot.col(0))
# print(Mtot.col(1))
# print(Mtot.col(2))

print(Ftot)
print(Mtot)
# print(equations)
print(sp.solve(equations, [rx, ry, rz]))
# print(sp.solve(equations, [rx]))
# print(sp.solve(equations, [ry]))
# print(sp.solve(equations, [rz]))
# print(sp.solve(equations, [rx, ry, rz, rax, ray, raz]))


