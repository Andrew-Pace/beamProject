import math
import matplotlib.pyplot as plt
import numpy as np

##### Dimensions ########
b_in = 6  # Width of the beam in feet
h_in = 12  # Height of the beam in feet
thickness_in = 0.675  # Thickness of the beam in inches

######## Material Properties ############
E = 30000000  # psi
alpha = 1e-5  # /degreeF
nu = 0.3  # Example value for Poisson's ratio

###### Case Variables #######
P_w = 0  # Wind Force, lbf
delta_T = 30  # degreeF
##########################
##########################
##########################



# unchanging parameters
L_in = 102  # Length of Beam, inches
P1 = 357  # P1, lbs
P2 = 234  # P2, lbs
x1_in = 26.68  # P1 position from left, inches
x2_in = x1_in + 42.125  # P2 position from left, inches
yl_ft = 1 # y distance of loads from beam feet
LL_lbsft2 = 75 # live load, lbs/ft2
uL_lbft = LL_lbsft2 * (b_in/12)  # live load->uniform load, lbs/ft

#### DimensionCalculationsForEquations ############
hw = ((h_in * 12) - (2 * thickness_in)) #in
y = b_in / 2  # y, in
z = h_in / 2  # z, in
I_z = 0.015835438449159193  # I_z
I_y = 0.0011864518413922365  # I_y
def calculate_J_I_beam(bf, tf, tw, d):
    """
    Calculate the polar moment of inertia (J) for an I-beam.

    Parameters:
    bf (float): Flange width (inches)
    tf (float): Flange thickness (inches)
    tw (float): Web thickness (inches)
    d (float): Total depth of the I-beam (inches)

    Returns:
    float: Polar moment of inertia (J) (inches^4)
    """
    # Height of the web
    h_web = d - 2 * tf

    # Polar moment of inertia for the flanges
    J_flange = 2 * ((bf * tf**3) / 3)

    # Polar moment of inertia for the web
    J_web = (tw * h_web**3) / 3

    # Total polar moment of inertia
    J_total = J_flange + J_web

    return J_total
J = calculate_J_I_beam(b_in, thickness_in, thickness_in, h_in)
x1_ft = x1_in / 12
x2_ft = x2_in / 12
L_ft = L_in / 12
Pw_lbft = P_w / L_ft


# Calculate vertical reactions at supports
R_cz = ((P1 * x1_ft) + (P2 * x2_ft) + (uL_lbft * L_ft * L_ft / 2)) / L_ft
R_az = ((P1 * (L_ft - x1_ft)) + (P2 * (L_ft - x2_ft)) + (uL_lbft * L_ft * L_ft / 2)) / L_ft
R_ay = P_w / 2

# Initialize lists to store results
xp_values = np.linspace(0, L_ft, num=500)
M_z_values = []
M_y_values = []
M_combined_values = []
V_z_values = []
V_y_values = []
sigma_x_values = []
tau_z_values = []
tau_y_values = []
tau_zy_values = []
tau_t_values = []
tau_combined_values = []
epsilon_t_values = []
sigma_Thermal_values = []
sigma_v_values = []
epsilon_x_values = []


def calculate_torque(xp, P1, P2, x1_ft, x2_ft, yl_ft):
    """
    Calculate the torque at a given point along the beam.

    Parameters:
    xp (float): Position along the beam (feet)
    P1 (float): Load P1 (pounds)
    P2 (float): Load P2 (pounds)
    x1_ft (float): Position of P1 (feet)
    x2_ft (float): Position of P2 (feet)
    yl_ft (float): Leverage distance (feet)

    Returns:
    float: Total torque at the given point (foot-pounds)
    """
    # Calculate torque contribution from P1
    if xp >= x1_ft:
        T1 = P1 * yl_ft
    else:
        T1 = P1 * yl_ft * (xp / x1_ft)  # Linear interpolation before x1_ft

    # Calculate torque contribution from P2
    if xp >= x2_ft:
        T2 = P2 * yl_ft
    else:
        T2 = P2 * yl_ft * (xp / x2_ft)  # Linear interpolation before x2_ft

    # Total torque at point xp
    T_total = T1 + T2
    return T_total

def calc_Q(z1):
    ### From Beam Equations Word doc
    Q = (b_in/8) * (h_in**2 - hw**2) + (thickness_in/8)*(hw**2-(4*(z1**2)))
    return Q

def calculate_Q_flange_horizontal(bf, tf):
    """
    Calculate the first moment of area (Q) for both flanges of an I-beam under horizontal shear force.

    Parameters:
    bf (float): Flange width (inches)
    tf (float): Flange thickness (inches)

    Returns:
    float: First moment of area (Q) (inches^3)
    """
    # Area of one flange
    A_flange = bf * tf
    # Distance to the centroid of one flange
    x_bar_flange = bf / 2
    # First moment of area for one flange
    Q_flange = A_flange * x_bar_flange
    # Since there are two flanges, multiply by 2
    Q_total = 2 * Q_flange
    return Q_total

for xp in xp_values:
    # Calculate bending moment M_z
    if xp < x1_ft:
        M_z = (R_az * xp) - (uL_lbft * (xp * xp) / 2)
    elif x1_ft <= xp < x2_ft:
        M_z = (R_az * xp) - (P1 * (xp - x1_ft)) - (uL_lbft * (xp * xp) / 2)
    else:
        M_z = (R_az * xp) - (P1 * (xp - x1_ft)) - (P2 * (xp - x2_ft)) - (uL_lbft * (xp * xp) / 2)

    M_z_values.append(M_z)

    M_y = R_ay * xp - (Pw_lbft * xp**2) / 2

    M_y_values.append(M_y)

    # Calculate combined bending moment M_combined
    M_combined = math.sqrt(M_z ** 2 + M_y ** 2)
    print(xp, M_combined)
    M_combined_values.append(M_combined)

    # Calculate shear forces V_z and V_y
    if xp < x1_ft:
        V_z = R_az - (uL_lbft * xp)
    elif x1_ft <= xp < x2_ft:
        V_z = R_az - P1 - (uL_lbft * xp)
    else:
        V_z = R_az - P1 - P2 - (uL_lbft * xp)

    V_z_values.append(V_z)

    V_y = R_ay - P_w * xp
    V_y_values.append(V_y)

    # Calculate bending stress σ_x
    sigma_x = -(M_z * y) / I_z + (M_y * z) / I_y
    sigma_x_values.append(sigma_x)

    # Calculate shear stresses τ_z and τ_y
    tau_z = (V_z * calc_Q(0)) / (I_z * thickness_in)
    tau_y = (V_y * calculate_Q_flange_horizontal(thickness_in, thickness_in)) / (I_y * thickness_in)

    tau_z_values.append(tau_z)
    tau_y_values.append(tau_y)

    # Calculate combined shear stress τ_zy and torsional shear stress τ_t
    tau_zy = math.sqrt(tau_z ** 2 + tau_y ** 2)
    tau_t = (calculate_torque(xp, P1, P2, x1_ft, x2_ft, yl_ft) * L_ft / 2) / J

    tau_zy_values.append(tau_zy)
    tau_t_values.append(tau_t)

    # Calculate combined shear stress τ_combined
    tau_combined = math.sqrt(tau_zy ** 2 + tau_t ** 2)
    tau_combined_values.append(tau_combined)

    # Calculate thermal strain ε_t and thermal stress σ_Thermal
    epsilon_t = alpha * delta_T
    sigma_Thermal = E * epsilon_t

    epsilon_t_values.append(epsilon_t)
    sigma_Thermal_values.append(sigma_Thermal)

    # Calculate maximum stress, von Mises stress σ_v
    sigma_zt = sigma_x + sigma_Thermal
    sigma_yt = sigma_x + sigma_Thermal

    sigma_v = math.sqrt(sigma_zt ** 2 + sigma_yt ** 2 - sigma_zt * sigma_yt + 3 * tau_combined ** 2)

    sigma_v_values.append(sigma_v)

    # Calculate strain ε_x using Hooke's Law
    epsilon_x = (sigma_x - nu * (sigma_x + sigma_x)) / E

    epsilon_x_values.append(epsilon_x) ##

# Dictionary to map argument names to their corresponding data arrays and labels
data_dict = {
    'M_z': (M_z_values, 'Bending Moment M_z (foot-pounds)'),
    'M_y': (M_y_values, 'Bending Moment M_y (foot-pounds)'),
    'M_combined': (M_combined_values, 'Combined Bending Moment M_combined (foot-pounds)'),
    'V_z': (V_z_values, 'Shear Force V_z (pounds)'),
    'V_y': (V_y_values, 'Shear Force V_y (pounds)'),
    'sigma_x': (sigma_x_values, 'Bending Stress σ_x (psi)'),
    'tau_z': (tau_z_values, 'Shear Stress τ_z (psi)'),
    'tau_y': (tau_y_values, 'Shear Stress τ_y (psi)'),
    'tau_zy': (tau_zy_values, 'Combined Shear Stress τ_zy (psi)'),
    'tau_t': (tau_t_values, 'Torsional Shear Stress τ_t (psi)'),
    'tau_combined': (tau_combined_values, 'Combined Shear Stress τ_combined (psi)'),
    'epsilon_t': (epsilon_t_values, 'Thermal Strain ε_t'),
    'sigma_Thermal': (sigma_Thermal_values, 'Thermal Stress σ_Thermal (psi)'),
    'sigma_v': (sigma_v_values, 'Von Mises Stress σ_v (psi)'),
    'epsilon_x': (epsilon_x_values, 'Strain ε_x')
}

def plot_value(xp_values, value_name, save_fig = False, mark_max = False):
    """
    Plot the specified value along the beam.

    Parameters:
    xp_values (array): Array of positions along the beam.
    value_name (str): Name of the value to plot.
    """
    if value_name in data_dict:
        values, ylabel = data_dict[value_name]
        plt.figure(figsize=(12, 4))  # Set the figure size to be short and wide
        plt.plot(xp_values, values, label=value_name)
        plt.xlabel('Position along the beam (feet)')
        plt.ylabel(ylabel)
        plt.title(f'{ylabel} along the beam')
        plt.legend()
        plt.grid(True)
        plt.xlim(0, xp_values[-1])
        tick_interval = 0.5  # Adjust this value to change the interval
        plt.xticks(np.arange(0, xp_values[-1] + tick_interval, tick_interval))
        if mark_max:
            # Find the maximum value and its position
            max_value = max(values)
            max_index = values.index(max_value)
            max_x = xp_values[max_index]

            # Annotate and mark the maximum value with just text
            plt.scatter(max_x, max_value, color='red', zorder=2)  # Mark the maximum point
            plt.text(max_x, max_value + 0.1, f'Max: {max_value:.2f}', ha='center', fontweight='bold')
        if save_fig:
            plt.savefig(f'{ylabel}')
        plt.show()
    else:
        print(f"Value name '{value_name}' not recognized. Please choose from: {list(data_dict.keys())}")


# Example usage: Plotting Bending Moment M_z along the beam
plot_value(xp_values, 'M_y', False, False)
