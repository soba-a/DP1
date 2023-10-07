import sympy as sp
import numpy as np
import math as m
# THIS IS THE ORIGINAL, 
# rx, ry, rz, warm, fact = sp.symbols('rx, ry, rz, warm, fact')
rx, ry, rz, fact = sp.symbols('rx, ry, rz, fact')

# givens:
wbre = 2 * 9.81 # N
Fresz = 40 # N

# pulley givens
pul_PD = 3 * 25.4 # in mm
blade_dia = 8 * 25.4 # in mm

# assume lpiv is in line with lblade #design
cpiv_off = -25

# THINGS TO CHANGE:############################################
# arm_density = 2.7 / 10**3 ########################## WAS NOT CHANGED BC WE HAVE UNIT WEIGHT
unit_weight = 26.6
armod = 2 * 25.4; # design of arm diameterm mm #############################
armid = 1.75 * 25.4 # inner diameter of arm, mm #############################
# ys = 47900 / 145.03773773 # in MPa # design
ys = 276

pd = 0.25 * 25.4 # pin diameter in mm #############################
Ep = 190 * 10**9 # modulus of elasticity in Pa, converted from GPa
vp = 0.28 # possion's ratio
Earm = 68.9 * 10**9 # modulus of elasticity in Pa, converted from GPa
varm = 0.33 # possion's ratio
# arm design:

L = 410.07
L_m = L / 1000 # in m

vol = (np.pi/4) * ((armod)**2 - (armid)**2) * 450 # in m^3
# warm = arm_density * vol * 9.81
warm = unit_weight * vol * 1000 / 1000**3 #using unit weight

# Iarm = np.pi * (armod)**4 / 64 # solid
Iarm = np.pi * (armod)**4 / 64 - np.pi  * (armid)**4 / 64 # hollow
Iarm_m = Iarm / (1000**4)
# print("I =", Iarm)
Jarm = np.pi * (armod)**4 / 32 - np.pi  * (armid)**4 / 32 # hollow
Jarm_m = Jarm / (1000**4)
# print("J =", Jarm)
# Aarm = np.pi * (armd)**2 / 4 # solid
Aarm = np.pi / 4 * ((armod)**2 - (armid)**2) # hollow
Aarm_m = Aarm / (1000**2)


ys_safety = ys / 2

# PIN DESIGN
# pin details: 4140 Alloy Steel

# arm details: 6061 Aluminum, 1/4" walls, 1" OD
# this on on the hole! using h8f7 tolerance
armw = (armod - armid) / 2 # wall thickness
ph = pd + 0.022 # pin hole diameter in mm


####### SOLVING HERE
#######
#######
#######
#######
#######
#######
#######
#######
#######
#######

# Solve for torque (1)
torqueBlade = 1 * 9550 / 4600 * 1000 # n *mm

# through blade fbd
Fresx = -torqueBlade / 101.6
# print('Fresx =',Fresx)

# through blade pulley fbd
T1 = 45
T2 = (torqueBlade + 38.1 * T1) / 38.1 
# print('T1 =',T1)
# print('T2 =',T2)
        
# Locations with respect to the pin 1
Rpiv_loc = sp.Matrix([0, 0, 0]).T
Ract_loc = sp.Matrix([L/2, 0, cpiv_off]).T
Warm_loc = sp.Matrix([L/2, 0, cpiv_off]).T
Ft1_loc = sp.Matrix([L, -100, cpiv_off + 25 + pul_PD/2]).T
Ft2_loc = sp.Matrix([L, -100, cpiv_off + 25 - pul_PD/2]).T
Wbre_loc = sp.Matrix([L, 0, cpiv_off + 25]).T
Fres_loc = sp.Matrix([L, 65, cpiv_off + 25 - blade_dia/2]).T

# angles of pulley forces, assuming pulleys are in line & actuator
beta = np.arctan(410.07/(50.8-38.1)) # pulley
theta = np.arctan((150 + (-1 * cpiv_off)) / (L/2)) # actuator

# Remember, FactV = Ract! These values are irregardless of location
# RpivV = sp.Matrix([rx, ry, rz]).T
RpivV = sp.Matrix([rx, 0, rz]).T # set ry = 0, resistance force on arm by pin1
RactV = sp.Matrix([fact * sp.cos(theta), 0, -fact * sp.sin(theta)]).T # resistance force on arm by pin2, we know it = fact
WarmV = sp.Matrix([0, 0, -warm]).T # weight of arm
Ft1V = sp.Matrix([-T1 * np.sin(beta), 0, T1 * np.cos(beta)]).T 
Ft2V = sp.Matrix([-T2 * np.sin(beta), 0, -T2 * np.cos(beta)]).T
WbreV = sp.Matrix([0, 0, -wbre]).T
FresV = sp.Matrix([Fresx, 0, Fresz]).T

# we need this system to solve for: rx, rz, fact, and reaction moments
# print("UNSOLVED FORCES =======================")
# print("RpivV =", RpivV)
# print("RactV =", RactV)
# print("WarmV =", WarmV)
# print("Ft1V =", Ft1V)
# print("Ft2V =", Ft2V)
# print("WbreV =", WbreV)
# print("FresV =", FresV)

# Moment calculated with respect to pin1
MRpiv = Rpiv_loc.cross(RpivV)
MRact = Ract_loc.cross(RactV)
MWarm = Warm_loc.cross(WarmV)
MFt1 = Ft1_loc.cross(Ft1V)
MFt2 = Ft2_loc.cross(Ft2V)
MWbre = Wbre_loc.cross(WbreV)
MFres = Fres_loc.cross(FresV)
Mpin_armx, Mpin_armz = sp.symbols('Mpin_armx, Mpin_armz')
Mpin_arm = sp.Matrix([Mpin_armx, 0, Mpin_armz]).T # Resistant moment on arm by pin
# print("UNSOLVED MOMENTS =======================")
# print("MRpiv =", MRpiv)
# print("MRact =", MRact)
# print("MWarm =", MWarm)
# print("MFt1 =", MFt1)
# print("MFt2 =", MFt2)
# print("MWbre =", MWbre)
# print("MFres =", MFres)

# solve for rx, rz, fact, and reaction moments
Ftot = RpivV + RactV + WarmV + Ft1V + Ft2V + WbreV + FresV
Mpin1 = MRpiv + MRact + MWarm + MFt1 + MFt2 + MWbre + MFres + Mpin_arm
# sysEq1 = Ftot.row_insert(1, Mpin1)

# print("SYSTEM =======================")
# print("Ftot =", Ftot) 
# print("Mpin1 =", Mpin1)

Mpin_armx_sol = sp.solve(Mpin1[0], [Mpin_armx])[0] # pin on arm reaction moment
Mpin_armz_sol = sp.solve(Mpin1[2], [Mpin_armz])[0] # pin on arm reaction moment
Mpin_arm = sp.Matrix([Mpin_armx_sol, 0, Mpin_armz_sol]).T # replace mpin_Arm with Mpin_arm_sol

# print("Solving for Fact =", sp.solve(Mpin1[1], [fact]))
fact_sol = sp.solve(Mpin1[1], [fact])[0] # find fact_sol

# print("Solving for reactions at pin1 =", sp.solve(Ftot, [rx, rz]))
rx_sol = sp.solve(Ftot, [rx, rz])[rx] # replace rx with rx_sol
rz_sol = sp.solve(Ftot, [rx, rz])[rz] # replace rz with rz_sol
rx_sol = rx_sol.subs([(fact, fact_sol)])
rz_sol = rz_sol.subs([(fact, fact_sol)])

RpivV = sp.Matrix([rx_sol, 0, rz_sol]).T # replace RpivV with RpivV_sol
RactV = sp.Matrix([fact_sol * sp.cos(theta), 0, -fact_sol * sp.sin(theta)]).T # replace RactV with RactV_sol
MRact = MRact.subs([(fact, fact_sol)]) # substitute value so fact is not in the equation

# Remember, FactV = Ract!
FactV = RactV

## BASED ON ABOVE, WE HAVE ALL FORCES AND MOMENTS ACTING ON THE ARM!!!
# PRINT THEM HERE:
# print("FORCES AND MOMENTS, SOLVED!=======================")
# print("RpivV =", RpivV)
# print("FactV =", FactV)
# print("WarmV =", WarmV)
# print("Ft1V =", Ft1V)
# print("Ft2V =", Ft2V)
# print("WbreV =", WbreV)
# print("FresV =", FresV)
# print("MRpiv =", MRpiv)
# print("MRact =", MRact)
# print("MWarm =", MWarm)
# print("MFt1 =", MFt1)
# print("MFt2 =", MFt2)
# print("MWbre =", MWbre)
# print("MFres =", MFres)
# print("Mpin_arm =", Mpin_arm)

############################################
# rotate coords https://www.vcalc.com/wiki/vCalc/Circle%20-%20Radius%20from%20chord%20length%20and%20arc%20height
r = 853.2870245 # radius of curvature
theta_prime = np.pi/2 - np.arcsin((r - 25) / r)

A = sp.Matrix([[np.cos(theta_prime), 0, np.sin(theta_prime)], [0, 1, 0], [-1 * np.sin(theta_prime), 0, np.cos(theta_prime)]])
# Forces in A coords
RpivV_A = (A * RpivV.T).T
FactV_A = (A * RactV.T).T
WarmV_A = (A * WarmV.T).T
Ft1V_A = (A * Ft1V.T).T
Ft2V_A = (A * Ft2V.T).T
WbreV_A = (A * WbreV.T).T
FresV_A = (A * FresV.T).T
# Moments in A coords
MRpiv_A = (A * MRpiv.T).T
MRact_A = (A * MRact.T).T
MWarm_A = (A * MWarm.T).T
MFt1_A = (A * MFt1.T).T
MFt2_A = (A * MFt2.T).T
MWbre_A = (A * MWbre.T).T
MFres_A = (A * MFres.T).T
Mpin_arm_A = (A * Mpin_arm.T).T

# print("ROTATED FORCES AND MOMENTS AT ORIGIN, SOLVED!=======================")
# print("RpivV_A =", RpivV_A)
# print("MRpiv_A =", MRpiv_A)
# print("MRact_A =", MRact_A)
# print("MWarm_A =", MWarm_A)
# print("MFt1_A =", MFt1_A)
# print("MFt2_A =", MFt2_A)
# print("MWbre_A =", MWbre_A)
# print("MFres_A =", MFres_A)
# print("Mpin_arm_A =", Mpin_arm_A)

# # compare the rotated force and moment values to the original in one line each
# print("RpivV =", RpivV)
# print("Rpiv_A =", RpivV_A)
# print("MRpiv =", MRpiv)
# print("MRpiv_A =", MRpiv_A)
# print("MRact =", MRact)
# print("MRact_A =", MRact_A)
# print("MWarm =", MWarm)
# print("MWarm_A =", MWarm_A)
# print("MFt1 =", MFt1)
# print("MFt1_A =", MFt1_A)
# print("MFt2 =", MFt2)
# print("MFt2_A =", MFt2_A)
# print("MWbre =", MWbre)
# print("MWbre_A =", MWbre_A)
# print("MFres =", MFres)
# print("MFres_A =", MFres_A)
# print("Mpin_arm =", Mpin_arm)
# print("Mpin_arm_A =", Mpin_arm_A)

############################################
# FINDING CLEVIS ON PIN REATIONS

# Let's look at pin1 here
p1rx1, p1rz1, p1rx2, p1rz2 = sp.symbols('p1rx1, p1rz1, p1rx2, p1rz2')
P1R1 = sp.Matrix([p1rx1, 0, p1rz1]).T
P1R2 = sp.Matrix([p1rx2, 0, p1rz2]).T
# P1R1_loc = sp.Matrix([0, -armd/2 + -0.25/2, 0]).T # assuming in the middle of the clevis !WRONG
# P1R2_loc = sp.Matrix([0, armd/2 + 0.25/2, 0]).T # assuming in the middle of the clevis !WRONG
P1R1_loc = sp.Matrix([0, -(armod)/2, 0]).T # assuming at the edge of the clevis (CORRECT)
P1R2_loc = sp.Matrix([0, (armod)/2, 0]).T # assuming at the edge of the clevis (CORRECT)

FTOT_P1 = P1R1 + P1R2 + -1 * RpivV # RpivV is the reaction force on the arm by pin1, so flip it to get force on pin by arm
MTOT_P1 = P1R1_loc.cross(P1R1) + P1R2_loc.cross(P1R2) + -1 * Mpin_arm # same as above, but with moments
# print("UNSOLVED SYSTEM FOR PIN1======================="")
# print("FTOT_P1 =", FTOT_P1)
# print("MTOT_P1 =", MTOT_P1)
sysEq2 = FTOT_P1.row_insert(1, MTOT_P1)

# SOLVE FOR REACTIONS
# print(sp.solve(equations2, [p1rx1, p1rz1, p1rx2, p1rz2]))
p1rx1_sol = sp.solve(sysEq2, [p1rx1, p1rz1, p1rx2, p1rz2])[p1rx1]
p1rz1_sol = sp.solve(sysEq2, [p1rx1, p1rz1, p1rx2, p1rz2])[p1rz1]
p1rx2_sol = sp.solve(sysEq2, [p1rx1, p1rz1, p1rx2, p1rz2])[p1rx2]
p1rz2_sol = sp.solve(sysEq2, [p1rx1, p1rz1, p1rx2, p1rz2])[p1rz2]
P1R1 = sp.Matrix([p1rx1_sol, 0, p1rz1_sol]).T # replace
P1R2 = sp.Matrix([p1rx2_sol, 0, p1rz2_sol]).T # replace

# let's look at actuator pin here
p2rx1, p2rz1, p2rx2, p2rz2 = sp.symbols('p2rx1, p2rz1, p2rx2, p2rz2')
P2R1 = sp.Matrix([p2rx1, 0, p2rz1]).T
P2R2 = sp.Matrix([p2rx2, 0, p2rz2]).T
# P2R1_loc = sp.Matrix([0, -armd/2 + -0.25/2, 0]).T # wrong because its at middle of clevis
# P2R2_loc = sp.Matrix([0, armd/2 + 0.25/2, 0]).T # wrong because its at middle of clevis
P2R1_loc = sp.Matrix([0, -(armod)/2, 0]).T # assuming at the edge of the clevis (CORRECT)
P2R2_loc = sp.Matrix([0, (armod)/2, 0]).T # assuming at the edge of the clevis (CORRECT)

FTOT_P2 = P2R1 + P2R2 + -1 * RactV.subs([(fact, fact_sol)]) # RactV is the reaction force on the arm by pin2, so flip it to get force on pin by arm
MP2 = P2R1_loc.cross(P2R1) + P2R2_loc.cross(P2R2) # assume no reaction moments on actuator pin because actuator is truss/2 force member. bc truss do not carry moment
sysEq3 = FTOT_P2.row_insert(1, MP2)

# SOLVE FOR REACTIONS
p2rx1_sol = sp.solve(sysEq3, [p2rx1, p2rz1, p2rx2, p2rz2])[p2rx1]
p2rz1_sol = sp.solve(sysEq3, [p2rx1, p2rz1, p2rx2, p2rz2])[p2rz1]
p2rx2_sol = sp.solve(sysEq3, [p2rx1, p2rz1, p2rx2, p2rz2])[p2rx2]
p2rz2_sol = sp.solve(sysEq3, [p2rx1, p2rz1, p2rx2, p2rz2])[p2rz2]

# print pin reaction forces! pin 1 = pivot | pin2 = actuator
# print("PIN REACTION FORCES! SOLVED! =======================")
# print("p1rx1_sol =", p1rx1_sol)
# print("p1rz1_sol =", p1rz1_sol)
# print("p1rx2_sol =", p1rx2_sol)
# print("p1rz2_sol =", p1rz2_sol)

# print("p2rx1_sol =", p2rx1_sol)
# print("p2rz1_sol =", p2rz1_sol)
# print("p2rx2_sol =", p2rx2_sol)
# print("p2rz2_sol =", p2rz2_sol)

###### PIN DESIGN CHECKS
# check pin 1
RpivV_norm = (RpivV[0]**2 + RpivV[1]**2 + RpivV[2]**2)**0.5
p1r1_norm = (p1rx1_sol**2 + p1rz1_sol**2)**0.5
p1r2_norm = (p1rx2_sol**2 + p1rz2_sol**2)**0.5

#check pin2
p2r1_norm = (p2rx1_sol**2 + p2rz1_sol**2)**0.5
p2r2_norm = (p2rx2_sol**2 + p2rz2_sol**2)**0.5

print("================== PIN DESIGN CHECKS ==================")
print("RpivV_norm =", RpivV_norm, "N")
print("p1r1_norm =", p1r1_norm, "N")
print("p1r2_norm =", p1r2_norm, "N")
print("p2r1_norm =", p2r1_norm, "N")
print("p2r2_norm =", p2r2_norm, "N")
tmax_pin = np.max([RpivV_norm, p1r1_norm, p1r2_norm, p2r1_norm, p2r2_norm]) / ((np.pi/4) * pd**2) # MPa
print("Max Shear Stress across both pins =", tmax_pin, "MPa")

###
###
###
### solve for maximum my bending moment
# Calculate My at actuator, from pov of actuator pin in this direction ->
rpivlocation = sp.Matrix([-L/2, 0, -1 * cpiv_off]).T
warmlocation = sp.Matrix([-L/4, 0, -1 * cpiv_off]).T
My_actuator = rpivlocation.cross(RpivV)
Mwarm_actuator = warmlocation.cross(WarmV)
My_actuator_total = My_actuator + Mwarm_actuator
# print("My_actuator_total =", My_actuator_total) # run this value

## Calculate My at blade in this direction ->
# rpivlocation = rpivlocation + sp.Matrix([-L/2, 0, 0]).T
# warmlocation = sp.Matrix([-L/2, 0, -1 * cpiv_off]).T

# My_actuator = rpivlocation.cross(RpivV)
# Mwarm_actuator = warmlocation.cross(WarmV)
# Mract_actuator = warmlocation.cross(RactV)

# My_actuator_total = My_actuator + Mwarm_actuator + Mract_actuator
# # print("My_actuator_total =", My_actuator_total)

# Limiting axial load (6)
peak_axial_load = RpivV_A[0] #N
peak_bend_moment_x = Mpin_arm_A[0] #Nmm
peak_bend_moment_y = My_actuator_total[1] #Nmm
peak_bend_moment_z = Mpin_arm_A[2] #Nmm
peak_torque = peak_bend_moment_x #nmm

print("================== PEAK LOADS ==================")
print("peak_axial_load =", peak_axial_load)
print("peak_bend_moment_x =", peak_bend_moment_x)
print("peak_bend_moment_y =", peak_bend_moment_y)
print("peak_bend_moment_z =", peak_bend_moment_z)

#### can we approximate and use the approximate formula for stress in curved beam?
rc = r # radius of curvature
rn = ((armod)/2)**2 / (2*(rc - (rc**2 - ((armod)/2)**2)**0.5))
# print rn and rc on one line
print("RN = RC ASSUMPTION CHECK ====================================")
print("rn =", rn, "rc =", rc)
# yea we can lol

# Calculate stress on beam
M = peak_bend_moment_z
y = (armod)/2

sigma_axial = peak_axial_load / Aarm
sigma_bending = M*y / Iarm * rc / ((armod)/2)
sigma_x = sigma_axial + sigma_bending
txy = peak_torque * ((armod)/2) / Jarm

vonmise_sigma = (sigma_x**2 + 3 * txy**2)**0.5

# stress at actuator pin
M = peak_bend_moment_y
sigma_bending_y = M*y / Iarm * rn / ((armod)/2)

### Calculate Stress Concentration at Pin
# Bending Stress Concentration, A-16 Shigley
print("================== STRESS CONCENTRATION AT PIN ==================")
# bending
d_D_ratio = armid/armod
print("d/D ratio =", d_D_ratio)
a_D_ratio = ph/armod
print("a/D ratio =", a_D_ratio)
print("===================================")
d_D_interp_vals = [0.6, 0.9]
a_D_interp_vals = [0.125, 0.15]
d_D_interp_A_val_1 = [0.82, 0.82] # top row
d_D_interp_K_val_1 = [2.32, 2.41] # top row
d_D_interp_A_val_2 = [0.79, 0.79] # bottom row
d_D_interp_K_val_2 = [2.29, 2.39] # bottom row
d_D_interp_val_A = [np.interp(d_D_ratio, d_D_interp_vals, d_D_interp_A_val_1), np.interp(d_D_ratio, d_D_interp_vals, d_D_interp_A_val_2)]
d_D_interp_val_K = [np.interp(d_D_ratio, d_D_interp_vals, d_D_interp_K_val_1), np.interp(d_D_ratio, d_D_interp_vals, d_D_interp_K_val_2)]
# print(d_D_interp_val_A, d_D_interp_val_K)
A_val_bending = np.interp(a_D_ratio, a_D_interp_vals, d_D_interp_val_A)
K_val_bending = np.interp(a_D_ratio, a_D_interp_vals, d_D_interp_val_K)
# print("A_val_bending =", A_val_bending)
# print("K_val_bending =", K_val_bending)

Znet = ((np.pi) * A_val_bending) / (32 * armod) * (armod**4 - armid**4)
M = peak_bend_moment_z
sigma_nom = M/Znet
stress_conc_bending = K_val_bending * sigma_nom
print("Stress Concentration Bending =", stress_conc_bending, "MPa")

# Torsion
d_D_interp_vals = [0.6, 0.8]
a_D_interp_vals = [0.125, 0.15]
d_D_interp_A_val_1 = [0.90, 0.91] # top row
d_D_interp_K_val_1 = [1.70, 1.74] # top row
d_D_interp_A_val_2 = [0.87, 0.89] # bottom row
d_D_interp_K_val_2 = [1.69, 1.75] # bottom row
d_D_interp_val_A = [np.interp(d_D_ratio, d_D_interp_vals, d_D_interp_A_val_1), np.interp(d_D_ratio, d_D_interp_vals, d_D_interp_A_val_2)]
d_D_interp_val_K = [np.interp(d_D_ratio, d_D_interp_vals, d_D_interp_K_val_1), np.interp(d_D_ratio, d_D_interp_vals, d_D_interp_K_val_2)]
# print(d_D_interp_val_A, d_D_interp_val_K)
A_val_torsion = np.interp(a_D_ratio, a_D_interp_vals, d_D_interp_val_A)
K_val_torsion = np.interp(a_D_ratio, a_D_interp_vals, d_D_interp_val_K)
# print("A_val_torsion =", A_val_torsion)
# print("K_val_torsion =", K_val_torsion)

Jnet = (np.pi * A_val_torsion * (armod**4 - armid**4)) / (32)
tau_nom = (peak_torque * armod) / (2*Jnet)
tau_conc_torsion = K_val_torsion * tau_nom
print("Stress Concentration Torsion =", tau_conc_torsion, "MPa")

# axial
c1 = 3
c2 = 0.427 - 6.770*(armid/armod) + 22.698*(armid/armod)**2 - 16.670*(armid/armod)**3
c3 = 11.357 + 15.665*(armid/armod) - 60.929*(armid/armod)**2 + 41.501*(armid/armod)**3
Ktg = c1 + c2 * (ph/armod) + c3*(ph/armod)**2
sigma_conc_axial = Ktg * sigma_axial
print("Stress Concentration Axial =", sigma_conc_axial, "MPa")


### Calculate pin arm contact stresses
Fapp = RpivV[0]
# print("Fapp =", Fapp)

l = (armw) #length of contact patch is 1/4 (wall thickness) mm

num = ((1-vp**2)/Ep) + ((1-varm**2)/Earm)
denom = (1/(pd/1000) - 1/(ph/1000))
b = ((2 * Fapp)/(np.pi * (l/1000)) * num/denom)**0.5
# print("CONTACT WIDTH: b =", 2 * b * 1000, "mm")
pmax = 2 * Fapp / (np.pi * (b*1000) * l) # n/mm**2

# print maximum stresses
print("================== ARM DESIGN CHECKS ==================")
print("Yield Safety =", ys_safety, "MPa")
print("===================================")
print("Maximum Axial Stress (Fx) @ Pivot Pin =", sigma_axial, "MPa")
print("Maximum Bending Stress (Mz) @ Pivot Pin =", sigma_bending, "MPa")

print("Maximum Torsional Stress (Mx) @ Pivot Pin =", txy, "MPa")
print("Von Mises Stress @ Pivot Pin =", vonmise_sigma, "MPa")
print("CONTACT STRESS @ Pivot Pin = pmax =", pmax, "MPa")

print("Maximum Bending Stress (My) @ Actuator Pin =", sigma_bending_y, "MPa")

####assume beam is straight for deflections
####
####
#### Deflection in Z direction, superposition of all
# Simple supports - overhanging load
FdefZ_tot = Ft1V[2] + Ft2V[2] + WbreV[2] + FresV[2]
w = warm  / (L_m)
Iarm_m = Iarm / (1000**4)
print("================== DEFLECTION CHECKS ==================")
totdef_z = w / (64 * Earm * Iarm_m) * L_m**4 - FdefZ_tot * L_m**2 / (6 * Earm * Iarm_m)
print("Z Deflection @ end of beam =", totdef_z * 1000, "mm")

### Deflection in y direction, superposition of all
totdef_y = (-1 * Mpin_arm[2]/1000) * L_m**2 / (2 * Earm * Iarm_m)
print("Y Deflection @ end of beam =", totdef_y*1000, "mm")

### Deflection in x direction, superposition of all
P = FactV[0] + Ft1V[0] + Ft2V[0] + FresV[0]
# print("Aarm =", Aarm)
totdef_x = P * L_m / (Aarm_m * Earm)
print("X Deflection @ end of beam =", totdef_x*1000, "mm")

## Torsional Deflection
# T = -1 * Mpin_arm[0] / 1000
T = (MFres[0] + MFt1[0] + MFt2[0]) / 1000
r = armod / 2 / 1000
J = Jarm_m
G = Earm / (2 * (1 + varm))
totdef_tor = T * r / (G * J)
print("TORSIONAL Deflection @ end of beam =", totdef_tor*180/np.pi, "deg")