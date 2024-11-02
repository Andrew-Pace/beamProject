import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from scipy.integrate import quad

factor_of_safety = 1
# Beam properties and loads
L = 9.125  # Length of the beam in feet
P1 = factor_of_safety * 357  # Point load 1 in lbs
a1 = 2.48  # Distance of point load 1 from the left support in feet
P2 = factor_of_safety * 234  # Point load 2 in lbs
a2 = L - 3.13  # Distance of point load 2 from the left support in feet
w_psf = factor_of_safety * 75  # Uniform load in lbs/ft^2
F_w = factor_of_safety * 0  # Total horizontal wind force in lbs
wind_force_perpindicular = 1615

##### Dimensions ########
b = 6 / 12  # Width of the beam in feet, using dimensions close to a S12 x 50 I beam from I beam table.
h = 12 / 12  # Height of the beam in feet
web_thickness = 0.675 # in
flange_thickness = 0.675 # in
Ixx = 0.015835438449159193 * 20736 # I ft4 -> I in4
Iyy = 0.0011864518413922365 * 20736 # I ft4 -> I in4
# Sxx = 50.8
# Syy = 5.74
E = 30000000  # psi
alpha = 10.8 * 10**-6 #Farenheight
deltaT = 0 # farenheight, using 0 for potential expansion joints to be used becuase thermal stress is huge for the normal 100F temp differnce
thickness = 0.675  # Thickness of the web&flanges in inches
##########################


binch = b * 12
hinch = h * 12
hw = ((h * 12) - (2 * thickness)) #in
# Convert uniform load from psf to plf to appply it evenly
w_plf = w_psf * b

# Cross-sectional properties of the beam
y = (h * 12) / 2  # Distance from neutral axis to outermost, inches
# Moment of inertia for I-beam
I = (thickness * (((h * 12) - (2 * thickness))**3) / 12) + b * ((h * 12)**3 - (h * 12 - (2 * thickness))**3)
A = 2*binch*thickness + hw*thickness  # Convert to square inches

# Young's modulus for the material (assuming 1095)
E = 30000000  # psi

def calculate_reactions(L, P1, a1, P2, a2, w_plf):
    # reactions at supports for a simply supp beam
    Ra = (P1 * (L - a1) + P2 * (L - a2) + w_plf * L * (L / 2)) / L
    Rb = P1 + P2 + w_plf * L - Ra
    return Ra, Rb

def calculate_shear_force(x, L, P1, a1, P2, a2, w_plf, Ra):
    # shear force at a distance x from the left supp
    V = 0
    if x < a1:
        V = Ra - w_plf * x
    elif x < a2:
        V = Ra - P1 - w_plf * x
    else:
        V = Ra - P1 - P2 - w_plf * x
    return V

def calculate_bending_moment(x, L, P1, a1, P2, a2, w_plf, Ra):
    # bending moment at a distance x from the left supp
    M = 0
    if x < a1:
        M = Ra * x - (w_plf * x**2) / 2
    elif x < a2:
        M = Ra * x - P1 * (x - a1) - (w_plf * x**2) / 2
    else:
        M = Ra * x - P1 * (x - a1) - P2 * (x - a2) - (w_plf * x**2) / 2
    return M

def calculate_bending_stress(Mt, yt, It):
    # Calculate bending stress
    # Mt comes in as ftlbs becasue of calcuilate_bending_moment function uses feet and lbs
    # converting by *12, cause y and I use inches
    return (Mt*12) * yt / It

def calculate_strain(stress, E):
    # Calculate strain using Hooke's Law
    return stress / E

def calc_Q(y1):
    ### From Beam Equations Word doc
    Q = ((b*12)/8) * ((h*12)**2 - hw**2) + (thickness/8)*(hw**2-(4*(y1**2)))
    return Q

def calc_Ic():
    ### From Beam Equations Word doc
    Ic = (binch*(hinch**3) - binch*(hw**3) + thickness*(hw**2))/12
    return Ic

def calc_tau(V, y1):
    Q = calc_Q(y1)
    Ic = calc_Ic()
    tau = V * Q / (Ic * thickness)
    return tau

def calc_tau_hor(V, y1):
    tau = (V*first_moment_of_area_web() / (centroidal_moment_of_inertia() * thickness)) * (binch*(hinch**2 - hw**2) + thickness*(hw**2 - (4*(y1**2))))
    return tau

def first_moment_of_area_web():
    """
    Calculate the first moment of area (Q) of the web for a horizontal I-beam.
    Parameters:
    bw (float): Web width
    tw (float): Web thickness
    tf (float): Flange thickness
    h (float): Total height of the I-beam
    Returns:
    float: First moment of area (Q) of the web
    """
    bf = binch
    tf = thickness
    tw = thickness
    # Height of the web
    h_web = h - 2 * tf
    # Area of the web
    A_web = tw * h_web
    # Distance from the centroid of the web to the neutral axis (centroid) of the I-beam
    y_web = h_web / 2 + tf
    # First moment of area (Q) of the web
    Q_hor = A_web * y_web
    return Q_hor

def centroidal_moment_of_inertia():
    """
    Calculate the centroidal moment of inertia for an I-beam.

    Parameters:
    bf (float): Flange width
    tf (float): Flange thickness
    bw (float): Web width
    tw (float): Web thickness
    h (float): Total height of the I-beam

    Returns:
    float: Centroidal moment of inertia
    """
    bf = binch
    tf = thickness
    tw = thickness
    # Area of the flanges
    A_flange = bf * tf

    # Area of the web
    A_web = tw * (h - 2 * tf)

    # Distance from the centroid of the flanges to the centroid of the I-beam
    d_flange = (h / 2) - (tf / 2)

    # Moment of inertia of the flanges about their own centroidal axis
    I_flange_centroid = (bf * tf ** 3) / 12

    # Moment of inertia of the web about its own centroidal axis
    I_web_centroid = (tw * (h - 2 * tf) ** 3) / 12

    # Parallel axis theorem to find the moment of inertia about the centroidal axis of the I-beam
    I_flange = 2 * (I_flange_centroid + A_flange * d_flange ** 2)
    I_web = I_web_centroid

    # Total centroidal moment of inertia
    Ic_hor = I_flange + I_web

    return Ic_hor

def calculate_horbending_moments(L, F, x):
    if 0 <= x <= L / 2:
        M = F * x / 2
    elif L / 2 < x <= L:
        M = F * (L - x) / 2
    else:
        raise ValueError("x must be within the range of 0 to L")

    return M

def calculate_horbending_stress(I, M, w):
    # Calculate the bending moment at the center of the beam
    c = w / 2
    bending_stress = M*12 * c / I
    return bending_stress

def calculate_von_mises_stress(horizontal_stress, vertical_stress):
    # Calculate the von Mises stress
    von_mises_stress = math.sqrt(horizontal_stress**2 + vertical_stress**2 - horizontal_stress * vertical_stress)
    return von_mises_stress

def calculate_shear_force_due_to_wind(L, w, x):
    """
    Calculate the shear force at any point x along the beam due to a horizontal wind load.

    Parameters:
    L (float): Length of the beam in feet.
    w (float or function): Wind load in pounds per linear foot (plf) or a function defining the load distribution.
    x (float): Position along the beam in feet.

    Returns:
    float: Shear force at point x in pounds.
    """
    if callable(w):
        # If w is a function, calculate the reaction forces by integrating the load distribution

        # Reaction force at the supports
        Ra, _ = quad(w, 0, L)
        Ra /= 2
    else:
        # If w is a constant, calculate the reaction forces directly
        Ra = w * L / 2

    # Shear force at any point x along the beam
    if callable(w):
        # Integrate the load distribution from 0 to x
        V, _ = quad(w, 0, x)
        V = Ra - V
    else:
        # Use the constant load formula
        V = Ra - w * x

    return V

def calculate_horizontal_shear_stress(V, I, t_web, h_web, b_flange, t_flange):
    """
    Calculate the horizontal shear stress at any point x along the beam for an I-beam.

    Parameters:
    V (float): Shear force at the point x in pounds.
    I (float): Moment of inertia of the entire cross-section in inches^4.
    t_web (float): Thickness of the web in inches.
    h_web (float): Height of the web in inches.
    b_flange (float): Width of the flange in inches.
    t_flange (float): Thickness of the flange in inches.

    Returns:
    float: Horizontal shear stress at the point x in psi.
    """
    # Calculate A' and y' for the I-beam
    A_prime = b_flange * t_flange + t_web * h_web / 2
    y_prime = (t_flange / 2 * b_flange * t_flange + h_web / 4 * t_web * h_web) / A_prime

    # Calculate Q
    Q = A_prime * y_prime

    # Calculate shear stress
    tau = V * Q / (I * t_web)

    return tau

def calculate_von_mises_shear_stress(tau_horizontal, tau_vertical):
    """
    Calculate the equivalent shear stress using the von Mises criterion.

    Parameters:
    tau_horizontal (float): Horizontal shear stress in psi.
    tau_vertical (float): Vertical shear stress in psi.

    Returns:
    float: Equivalent shear stress in psi.
    """
    tau_v = math.sqrt(tau_horizontal ** 2 + tau_vertical ** 2 - tau_horizontal * tau_vertical)
    return tau_v

def calculate_thermal_stress(alpha, delta_T, E):
    """
    Calculate the thermal stress given the coefficient of thermal expansion, temperature change, and Young's modulus.

    Parameters:
    alpha (float): Coefficient of thermal expansion (1/°F or 1/°C).
    delta_T (float): Change in temperature (°F or °C).
    E (float): Young's modulus in psi.

    Returns:
    float: Thermal stress in psi.
    """
    epsilon_t = alpha * delta_T
    sigma_t = E * epsilon_t
    return sigma_t

def calculate_combined_stress(sigma_bending, tau_shear, sigma_thermal):
    """
    Calculate the combined stress using the von Mises criterion.

    Parameters:
    sigma_bending (float): Bending stress in psi.
    tau_shear (float): Shear stress in psi.
    sigma_thermal (float): Thermal stress in psi.

    Returns:
    float: Combined stress in psi.
    """
    sigma_v = math.sqrt(sigma_bending ** 2 + 3 * tau_shear ** 2 + sigma_thermal ** 2)
    return sigma_v


Ra, Rb = calculate_reactions(L, P1, a1, P2, a2, w_plf)
x_values = np.linspace(0, L, num=100)

vert_shear_forces = [calculate_shear_force(x, L, P1, a1, P2, a2, w_plf, Ra) for x in x_values]
vert_shear_stresses = [calc_tau(V, 0) for V in vert_shear_forces]

vert_moments = [calculate_bending_moment(x, L, P1, a1, P2, a2, w_plf, Ra) for x in x_values]
vert_bending_stresses = [calculate_bending_stress(M, y, Ixx) for M in vert_moments]

hor_bending_moments = [calculate_horbending_moments(L, wind_force_perpindicular, x) for x in x_values]
hor_bending_stresses = [calculate_horbending_stress(Iyy, M, h) for M in hor_bending_moments]
horizontal_shear_forces = [calculate_shear_force_due_to_wind(L, wind_force_perpindicular/L, x) for x in x_values]
horizontal_shear_stresses = [calculate_horizontal_shear_stress(V, Iyy, web_thickness, hw, b, flange_thickness) for V in horizontal_shear_forces]

thermal_compression_stress = calculate_thermal_stress(alpha, deltaT, E)

resultant_maxbending_stresses = [calculate_von_mises_stress(hor_stress, vert_stress)
                                 for hor_stress, vert_stress in zip(hor_bending_stresses, vert_bending_stresses)]
resulting_maxshear_stresses = [calculate_von_mises_shear_stress(tau_horizontal, tau_vertical)
                               for tau_horizontal, tau_vertical in zip(horizontal_shear_stresses, vert_shear_stresses)]

total_stresses = [calculate_combined_stress(sigma_bending, tau_shear, thermal_compression_stress) for sigma_bending, tau_shear in
                  zip(resultant_maxbending_stresses, resulting_maxshear_stresses)]
strains = [calculate_strain(sigma_bending, E) for sigma_bending in resultant_maxbending_stresses]

# for i, x in enumerate(x_values):
#     print(f"x: {x:.2f}ft, Bending Stress: {vert_bending_stresses[i]:.2f}psi,"
#           f" Hor Bending Stress: {hor_bending_stresses[i]:.2f}psi, Total Bending Stress: {resultant_maxbending_stresses[i]:.2f}psi,"
#           f" HorShearStress: {horizontal_shear_stresses[i]:.2f}psi")
# for i, x in enumerate(x_values):
#     print(f"x: {x:.2f}ft, thermalstresses: {thermal_compression_stress[i]:.2f}")

print(thermal_compression_stress)

# Function to find the closest index
def find_closest_index(array, value):
    return (np.abs(array - value)).argmin()

# Find the indices for the required x values
index_L2 = find_closest_index(x_values, L / 2)
index_second = 1
index_second_last = -2

# Print the results for x = L / 2
print(f"Results at x = L / 2 (x = {x_values[index_L2]}):")
print(f"Vertical Shear Force: {vert_shear_forces[index_L2]}")
print(f"Vertical Shear Stress: {vert_shear_stresses[index_L2]}")
print(f"Vertical Moment: {vert_moments[index_L2]}")
print(f"Vertical Bending Stress: {vert_bending_stresses[index_L2]}")
print(f"Horizontal Bending Moment: {hor_bending_moments[index_L2]}")
print(f"Horizontal Bending Stress: {hor_bending_stresses[index_L2]}")
print(f"Horizontal Shear Force: {horizontal_shear_forces[index_L2]}")
print(f"Horizontal Shear Stress: {horizontal_shear_stresses[index_L2]}")
print(f"Resultant Max Bending Stress: {resultant_maxbending_stresses[index_L2]}")
print(f"Resulting Max Shear Stress: {resulting_maxshear_stresses[index_L2]}")
print(f"Total Stress: {total_stresses[index_L2]}")
print(f"Strain: {strains[index_L2]}")

# Print the results for the second x value
print(f"\nResults at the second x value (x = {x_values[index_second]}):")
print(f"Vertical Shear Force: {vert_shear_forces[index_second]}")
print(f"Vertical Shear Stress: {vert_shear_stresses[index_second]}")
print(f"Vertical Moment: {vert_moments[index_second]}")
print(f"Vertical Bending Stress: {vert_bending_stresses[index_second]}")
print(f"Horizontal Bending Moment: {hor_bending_moments[index_second]}")
print(f"Horizontal Bending Stress: {hor_bending_stresses[index_second]}")
print(f"Horizontal Shear Force: {horizontal_shear_forces[index_second]}")
print(f"Horizontal Shear Stress: {horizontal_shear_stresses[index_second]}")
print(f"Resultant Max Bending Stress: {resultant_maxbending_stresses[index_second]}")
print(f"Resulting Max Shear Stress: {resulting_maxshear_stresses[index_second]}")
print(f"Total Stress: {total_stresses[index_second]}")
print(f"Strain: {strains[index_second]}")

# Print the results for the second-to-last x value
print(f"\nResults at the second-to-last x value (x = {x_values[index_second_last]}):")
print(f"Vertical Shear Force: {vert_shear_forces[index_second_last]}")
print(f"Vertical Shear Stress: {vert_shear_stresses[index_second_last]}")
print(f"Vertical Moment: {vert_moments[index_second_last]}")
print(f"Vertical Bending Stress: {vert_bending_stresses[index_second_last]}")
print(f"Horizontal Bending Moment: {hor_bending_moments[index_second_last]}")
print(f"Horizontal Bending Stress: {hor_bending_stresses[index_second_last]}")
print(f"Horizontal Shear Force: {horizontal_shear_forces[index_second_last]}")
print(f"Horizontal Shear Stress: {horizontal_shear_stresses[index_second_last]}")
print(f"Resultant Max Bending Stress: {resultant_maxbending_stresses[index_second_last]}")
print(f"Resulting Max Shear Stress: {resulting_maxshear_stresses[index_second_last]}")
print(f"Total Stress: {total_stresses[index_second_last]}")
print(f"Strain: {strains[index_second_last]}")


def plot_beam_analysis(x_values, shear_forces, shear_max_stresses, moments, bending_stresses):
    plt.figure(figsize=(14, 10))

    plt.subplot(5, 1, 1)
    plt.plot(x_values, shear_forces)
    plt.title('Shear Force along the Beam')
    plt.xlabel('Position along the beam (ft)')
    plt.ylabel('Shear Force (lbs)')
    plt.grid(True)

    # plt.subplot(5, 1, 2)
    # plt.plot(x_values, shear_max_stresses)
    # plt.title('Shear Max Stress along the Beam')
    # plt.xlabel('Position along the beam (ft)')
    # plt.ylabel('Shear Max Stress (psi)')
    # plt.grid(True)

    plt.subplot(5, 1, 3)
    plt.plot(x_values, moments)
    plt.title('Bending Moment along the Beam')
    plt.xlabel('Position along the beam (ft)')
    plt.ylabel('Bending Moment (lb-ft)')
    plt.grid(True)

    # plt.subplot(5, 1, 4)
    # plt.plot(x_values, bending_stresses)
    # plt.title('Bending Stress along the Beam')
    # plt.xlabel('Position along the beam (ft)')
    # plt.ylabel('Bending Stress (psi)')
    # plt.grid(True)


    plt.tight_layout()
    plt.savefig('beam_analysis.png')
    plt.show()
#
#
# def save_to_excel(x_values, shear_forces, shear_max_stresses, moments, bending_stresses, strains):
#     file_name = input("Enter the name of the Excel file to save results: ")
#
#     data = {
#         'Position along the beam (ft)': x_values,
#         'Shear Force (lbs)': shear_forces,
#         'Shear Max Stress (psi)': shear_max_stresses,
#         'Bending Moment (lb-ft)': moments,
#         'Total Bending Stress (psi)': bending_stresses,
#         'Strains': strains,
#     }
#
#     df = pd.DataFrame(data)
#     df.to_excel(f'{file_name}.xlsx', index=False)
#
#     print(f"Results have been saved to '{file_name}.xlsx'.")
#
# # # Example usage:
plot_beam_analysis(x_values, shear_forces, shear_max_stresses, moments, resultant_maxbending_stresses)
save_to_excel(x_values, shear_forces, shear_max_stresses, moments, resultant_maxbending_stresses, strains)

