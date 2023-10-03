import sympy as sp
import numpy as np
import math as m

rx, ry, rz, warm, fact = sp.symbols('rx, ry, rz, warm, fact')
# givens
wbre = 2 * 9.81
Fresz = 40
L = 410.07
pul_PD = 3 * 25.4 # in mm
blade_dia = 8 * 25.4 # in mm

# assume lpiv is in line with lblade #design
cpiv_off = -25

# Solve for torque (1)
torqueBlade = 1 * 9550 / 4600 * 1000 # n *mm

# through blade fbd
Fresx = -torqueBlade / 101.6
print('Fresx =',Fresx)

# through blade pulley fbd
T1 = 45
T2 = (torqueBlade + 38.1 * T1) / 38.1 
print('T1 =',T1)
print('T2 =',T2)
        
# Locations with respect to the pin 1
Rpiv_loc = sp.Matrix([0, 0, 0]).T
Fact_loc = sp.Matrix([L/2, 0, cpiv_off]).T
Warm_loc = sp.Matrix([L/2, 0, cpiv_off]).T
Ft1_loc = sp.Matrix([L, -100, cpiv_off + 25 + pul_PD/2]).T
Ft2_loc = sp.Matrix([L, -100, cpiv_off + 25 - pul_PD/2]).T
Wbre_loc = sp.Matrix([L, 0, cpiv_off + 25]).T
Fres_loc = sp.Matrix([L, 65, cpiv_off + 25 - blade_dia/2]).T

# angle of pulley forces, assuming pulleys are in line
beta = np.arctan(410.07/(50.8-38.1))
theta = np.arctan((150 + (-1 * cpiv_off)) / (L/2))

# Remember, FactV = Ract! These values are irregardless of location
# RpivV = sp.Matrix([rx, ry, rz]).T
RpivV = sp.Matrix([rx, 0, rz]).T # set ry = 0
FactV = sp.Matrix([fact * sp.cos(theta), 0, -fact * sp.sin(theta)]).T
WarmV = sp.Matrix([0, 0, -warm]).T
Ft1V = sp.Matrix([-T1 * np.sin(beta), 0, -T1 * np.cos(beta)]).T
Ft2V = sp.Matrix([-T2 * np.sin(beta), 0, T2 * np.cos(beta)]).T
WbreV = sp.Matrix([0, 0, -wbre]).T
FresV = sp.Matrix([Fresx, 0, Fresz]).T
print("RpivV =", RpivV)
print("FactV =", FactV)
print("WarmV =", WarmV)
print("Ft1V =", Ft1V)
print("Ft2V =", Ft2V)
print("WbreV =", WbreV)
print("FresV =", FresV)

# Moment calculated with respect to pin1
MRpiv = Rpiv_loc.cross(RpivV)
MFact = Fact_loc.cross(FactV)
MWarm = Warm_loc.cross(WarmV)
MFt1 = Ft1_loc.cross(Ft1V)
MFt2 = Ft2_loc.cross(Ft2V)
MWbre = Wbre_loc.cross(WbreV)
MFres = Fres_loc.cross(FresV)
print("MRpiv =", MRpiv)
print("MFact =", MFact)
print("MWarm =", MWarm)
print("MFt1 =", MFt1)
print("MFt2 =", MFt2)
print("MWbre =", MWbre)
print("MFres =", MFres)

Ftot = RpivV + FactV + WarmV + Ft1V + Ft2V + WbreV + FresV
Mpin1 = MRpiv + MFact + MWarm + MFt1 + MFt2 + MWbre + MFres
equations = Ftot.row_insert(1, Mpin1)

print("Ftot =", Ftot)
print("Mpin1 =", Mpin1)

# print(equations)
print("Solving for Fact =", sp.solve(Mpin1[1], [fact]))
fact_sol = sp.solve(Mpin1[1], [fact])[0] # replace fact with fact_sol

MpinR = sp.Matrix([-1 * Mpin1[0], 0, -1 * Mpin1[2]]).T
print("MpinR =", MpinR)

print("Solving for reactions at pin1 =", sp.solve(Ftot, [rx, rz]))
rx_sol = sp.solve(Ftot, [rx, rz])[rx] # replace rx with rx_sol
rz_sol = sp.solve(Ftot, [rx, rz])[rz] # replace rz with rz_sol
# Remember, FactV = Ract!
RactV = FactV

# Let's look at pin1 here
armd = 2 * 25.4; # design
p1rx1, p1rz1, p1rx2, p1rz2 = sp.symbols('p1rx1, p1rz1, p1rx2, p1rz2')
P1R1 = sp.Matrix([p1rx1, 0, p1rz1]).T
P1R2 = sp.Matrix([p1rx2, 0, p1rz2]).T
P1R1_loc = sp.Matrix([0, -armd/2 + -0.25/2, 0]).T
P1R2_loc = sp.Matrix([0, armd/2 + 0.25/2, 0]).T

FTOT_P1 = P1R1 + P1R2 + -1 * RpivV.subs([(rx, rx_sol), (rz, rz_sol), (fact, fact_sol)])
MTOT_P1 = P1R1_loc.cross(P1R1) + P1R2_loc.cross(P1R2) + -1 * MpinR
print("FTOT_P1 =", FTOT_P1)
print("MTOT_P1 =", MTOT_P1)
equations2 = FTOT_P1.row_insert(1, MTOT_P1)

# print(sp.solve(equations2, [p1rx1, p1rz1, p1rx2, p1rz2]))
p1rx1_sol = sp.solve(equations2, [p1rx1, p1rz1, p1rx2, p1rz2])[p1rx1]
p1rz1_sol = sp.solve(equations2, [p1rx1, p1rz1, p1rx2, p1rz2])[p1rz1]
p1rx2_sol = sp.solve(equations2, [p1rx1, p1rz1, p1rx2, p1rz2])[p1rx2]
p1rz2_sol = sp.solve(equations2, [p1rx1, p1rz1, p1rx2, p1rz2])[p1rz2]

# Since shear stress is function of faverage,
# print((p1rx1_sol + p1rx2_sol)/2)

# let's look at actuator pin here
p2rx1, p2rz1, p2rx2, p2rz2 = sp.symbols('p2rx1, p2rz1, p2rx2, p2rz2')
P2R1 = sp.Matrix([p2rx1, 0, p2rz1]).T
P2R2 = sp.Matrix([p2rx2, 0, p2rz2]).T
P2R1_loc = sp.Matrix([0, -armd/2 + -0.25/2, 0]).T
P2R2_loc = sp.Matrix([0, armd/2 + 0.25/2, 0]).T
MP2 = P2R1_loc.cross(P2R1) + P2R2_loc.cross(P2R2)

FTOT_P2 = P2R1 + P2R2 + -1 * RactV.subs([(fact, fact_sol)])
equations3 = FTOT_P2.row_insert(1, MP2)
p2rx1_sol = sp.solve(equations3, [p2rx1, p2rz1, p2rx2, p2rz2])[p2rx1]
p2rz1_sol = sp.solve(equations3, [p2rx1, p2rz1, p2rx2, p2rz2])[p2rz1]
p2rx2_sol = sp.solve(equations3, [p2rx1, p2rz1, p2rx2, p2rz2])[p2rx2]
p2rz2_sol = sp.solve(equations3, [p2rx1, p2rz1, p2rx2, p2rz2])[p2rz2]


# rotate coords https://www.vcalc.com/wiki/vCalc/Circle%20-%20Radius%20from%20chord%20length%20and%20arc%20height
r = 853.2870245
theta_prime = np.pi/2 - np.arcsin((r - 25) / r)

A = sp.Matrix([[np.cos(theta_prime), 0, np.sin(theta_prime)], [0, 1, 0], [-1 * np.sin(theta_prime), 0, np.cos(theta_prime)]])
# Forces in A coords
RpivV_A = (A * RpivV.T).T
FactV_A = (A * FactV.T).T
WarmV_A = (A * WarmV.T).T
Ft1V_A = (A * Ft1V.T).T
Ft2V_A = (A * Ft2V.T).T
WbreV_A = (A * WbreV.T).T
FresV_A = (A * FresV.T).T
# Moments in A coords
MRpiv_A = (A * MRpiv.T).T
MFact_A = (A * MFact.T).T
MWarm_A = (A * MWarm.T).T
MFt1_A = (A * MFt1.T).T
MFt2_A = (A * MFt2.T).T
MWbre_A = (A * MWbre.T).T
MFres_A = (A * MFres.T).T

print("RpivV_A =", RpivV_A)
print("FactV_A =", FactV_A)
print("WarmV_A =", WarmV_A)
print("Ft1V_A =", Ft1V_A)
print("Ft2V_A =", Ft2V_A)
print("WbreV_A =", WbreV_A)
print("FresV_A =", FresV_A)
print("MRpiv_A =", MRpiv_A)
print("MFact_A =", MFact_A)
print("MWarm_A =", MWarm_A)
print("MFt1_A =", MFt1_A)
print("MFt2_A =", MFt2_A)
print("MWbre_A =", MWbre_A)
