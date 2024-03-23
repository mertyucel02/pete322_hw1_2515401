from matplotlib import pyplot as plt
import math

print("***ANNULAR PRESSURE DROP AND ECD INPUTS***")

TVD = float(input("Enter TVD (ft): "))
rho = float(input("Enter fluid density (lbs/gal): "))
q = float(input("Enter volumetric flow rate (gal/min): "))
theta_300 = float(input("Enter 300-rpm dial reading: "))
theta_3 = float(input("Enter 3-rpm dial reading: "))
theta_600 = float(input("Enter 600-rpm dial reading: "))
ann_geo_num = int(input("Enter the # of annular geometries: "))
print("-------------------------------------------------------------")
list_d_2 = []
list_d_1 = []
list_L = []

for num in range(0, ann_geo_num):
    if num == 0:
        print("***INFO ABOUT 1ST ANNULAR GEOMETRY***")
    elif num == 1:
        print("***INFO ABOUT 2ND ANNULAR GEOMETRY***")
    elif num == 2:
        print("***INFO ABOUT 3RD ANNULAR GEOMETRY***")
    else:
        print("***INFO ABOUT {}TH ANNULAR GEOMETRY***".format(num + 1))
    d_2 = float(input("Enter outside diameter (in): "))
    list_d_2.append(d_2)
    d_1 = float(input("Enter inside diameter (in): "))
    list_d_1.append(d_1)
    length = float(input("Enter length (ft): "))
    list_L.append(length)
    print("-------------------------------------------------------------")

print("-------------------------------------------------------------")
print("-------------------------------------------------------------")

P_aT = 0
list_P_a = []

for i in range(0, ann_geo_num):
    V_a = (.408 * q) / ((list_d_2[i] ** 2) - (list_d_1[i] ** 2))
    n_a = .5 * math.log10(theta_300 / theta_3)
    K_a = (5.11 * theta_600) / (511 ** n_a)
    mu_ea = 100 * K_a * (((144 * V_a) / (list_d_2[i] - list_d_1[i])) ** (n_a - 1))
    Re_a = (928 * V_a * (list_d_2[i] - list_d_1[i]) * rho) / (
        mu_ea * (((2 * n_a + 1) / (3 * n_a)) ** n_a)
    )
    Re_L = 3470 - 1370 * n_a
    Re_T = 4270 - 1370 * n_a
    if Re_a < Re_L:
        f_a = 24 / Re_a
    elif Re_a > Re_T:
        f_a = ((math.log10(n_a) + 3.93) / 50) / (
            Re_a ** ((1.75 - math.log10(n_a)) / 7)
        )
    else:
        f_a = (((Re_a - Re_L) / 800) * (
                (((math.log10(n_a) + 3.93) / 50) / (Re_T ** (
                        (1.75 - math.log10(n_a)) / 7))) - (24 / Re_L))) + (24 / Re_L)
    P_a = (f_a * (V_a ** 2) * rho * list_L[i]) / (25.81 * (list_d_2[i] - list_d_1[i]))
    list_P_a.append(P_a)
    P_aT += P_a
    V_c = ((Re_L * K_a * (((2 * n_a + 1) / (3 * n_a)) ** n_a)) / (928 * rho * (list_d_2[i] - list_d_1[i]) * ((144 / (list_d_2[i] - list_d_1[i])) ** (1 - n_a)))) ** (1 / (2 - n_a))
    q_c = 2.45 * V_c * (list_d_2[i] ** 2 - list_d_1[i] ** 2)

ECD = (P_aT / (.052 * TVD)) + rho

print("+++++COMPUTED+++++")

print("***DRILL STRING PRESSURE DROP INPUTS***")

pipe_num = int(input("Enter # of pipes: "))
list_d_1_ds = []
list_L_ds = []
for num in range(0, pipe_num):
    if num == 0:
        print("***INFO ABOUT 1ST PIPE***")
    elif num == 1:
        print("***INFO ABOUT 2ND PIPE***")
    elif num == 2:
        print("***INFO ABOUT 3RD PIPE***")
    else:
        print("***INFO ABOUT {}TH PIPE***".format(num + 1))
    d_1_ds = float(input("Enter inside diameter (in): "))
    list_d_1_ds.append(d_1_ds)
    length_ds = float(input("Enter length (ft): "))
    list_L_ds.append(length_ds)
    print("-------------------------------------------------------------")

print("-------------------------------------------------------------")
print("-------------------------------------------------------------")

P_pT = 0
list_P_p = []

for i in range(0, pipe_num):
    V_p = (.408 * q) / (list_d_1_ds[i] ** 2)
    n_p = 3.32 * math.log10(theta_600 / theta_300)
    K_p = (5.11 * theta_600) / (1022 ** n_p)
    mu_ep = 100 * K_p * (((144 * V_p) / (list_d_1_ds[i])) ** (n_p - 1))
    Re_p = (928 * V_p * (list_d_1_ds[i]) * rho) / (
        mu_ep * (((3 * n_p + 1) / (4 * n_p)) ** n_p)
    )
    Re_L = 3470 - 1370 * n_p
    Re_T = 4270 - 1370 * n_p
    if Re_p < Re_L:
        f_p = 16 / Re_p
    elif Re_p > Re_T:
        f_p = ((math.log10(n_p) + 3.93) / 50) / (
            Re_p ** ((1.75 - math.log10(n_p)) / 7)
        )
    else:
        f_p = (((Re_p - Re_L) / 800) * (
                (((math.log10(n_p) + 3.93) / 50) / (Re_T ** (
                        (1.75 - math.log10(n_p)) / 7))) - (16 / Re_L))) + (16 / Re_L)
    P_p = (f_p * (V_p ** 2) * rho * list_L_ds[i]) / (25.81 * (list_d_1_ds[i]))
    list_P_p.append(P_p)
    P_pT += P_p


print("+++++COMPUTED+++++")
print("-------------------------------------------------------------")
print("-------------------------------------------------------------")

print("***BIT HYDRAULICS INPUTS***")
jet_num = int(input("Enter # of jets: "))
D_jet = float(input("Enter diameter of jet: "))
DC_d_2 = float(input("Enter outside diameter of DC: "))
DC_d_1 = float(input("Enter inside diameter of DC: "))
A_t = jet_num * (math.pi * D_jet / 4)
C_d = .95
delta_P_b = (8.311 * math.pow(10, -5) * rho * (q ** 2)) / ((C_d * A_t) ** 2)
HHP_b = q * delta_P_b / 1714
HHP_b_sq_in = HHP_b * 1.27 / (DC_d_1 ** 2)
print("+++++COMPUTED+++++")
print("-------------------------------------------------------------")
print("-------------------------------------------------------------")

print("***JET IMPACT FORCE***")
IF_jet = .01823 * C_d * q * math.sqrt(delta_P_b * delta_P_b)
print("+++++COMPUTED+++++")
print("-------------------------------------------------------------")
print("-------------------------------------------------------------")

print("***OPTIMUM NOZZLE SIZE INPUTS***")
q_1 = float(input("Enter 1st test flow rate: "))
P_1 = float(input("Enter 1st test pressure: "))
q_2 = float(input("Enter 2nd test flow rate: "))
P_2 = float(input("Enter 2nd test pressure: "))
pump_eff = float(input("Enter pump efficiency: "))
max_pump_HP = float(input("Enter maximum pump horsepower: "))
p_max = float(input("Enter maximum pump pressure: "))
q_min = float(input("Enter minimum flow rate to lift cuttings: "))

delta_P_b_1 = (8.311 * math.pow(10, -5) * rho * (q_1 ** 2)) / ((C_d * A_t) ** 2)
delta_P_b_2 = (8.311 * math.pow(10, -5) * rho * (q_2 ** 2)) / ((C_d * A_t) ** 2)
P_para_1 = P_1 - delta_P_b_1
P_para_2 = P_2 - delta_P_b_2
x_axis = [q_1, q_2]
y_axis = [P_para_1, P_para_2]
plt.plot(x_axis, y_axis)
plt.xlabel("Flow Rate ($gal/min$)")
plt.ylabel("Parasitic Pressure Loss ($psi$)")
plt.title("The Plot of Parasitic Pressure Loss ($psi$) versus Flow Rate ($gal/min$)")
for xy in zip(x_axis, y_axis):
    plt.annotate('($%.2f psi, %.2f gal/min$)' % xy, xy=xy)

m = math.log10(P_para_1 / P_para_2) / math.log10(q_1 / q_2)
slope = (P_para_2 - P_para_1) / (q_2 - q_1)
q_max = 1714 * max_pump_HP * pump_eff / p_max
delta_p_d = p_max * (2 / (m + 2))
plt.vlines(q_max, ymin = 0, ymax = delta_p_d, colors = 'r')
plt.hlines(delta_p_d, xmin=q_min, xmax=q_max, colors= 'r')
plt.vlines(q_min, ymin = delta_p_d, colors = 'r')
plt.show()

intersection = delta_P_b_2 + slope * (0 - delta_P_b_2)

delta_p_d_opt = slope * q_max + intersection
delta_P_b_opt = delta_p_d_opt - p_max
A_t_opt = math.sqrt((8.311 * math.pow(10, -5) * rho * (q_max ** 2)) / ((C_d ** 2) * delta_P_b_opt))
d_N_opt = math.sqrt((4 * A_t_opt) / (math.pi * D_jet))
print("+++++COMPUTED+++++")
print("-------------------------------------------------------------")
print("-------------------------------------------------------------")

print("***STANDPIPE PRESSURE***")
standpipe = P_aT + P_pT + delta_P_b
print("+++++COMPUTED+++++")
print("-------------------------------------------------------------")
print("-------------------------------------------------------------")
print("***RESULTS***")
print("Standpipe pressure is {} psi".format(standpipe))
print("Equivalent Circulating Density (ECD) is {} lbs/gal".format(ECD))
print("Bit pressure drop is {} psi".format(delta_P_b))
print("Horsepower per square inch (HSI) is {} HHP/in^2".format(HHP_b_sq_in))
print("Jet Impact Force is {} lbf".format(IF_jet))
print("Optimum Nozzle Size is {} in.".format(d_N_opt))
