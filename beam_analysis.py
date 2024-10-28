import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Beam properties and loads
L = 8.5  # Length of the beam in feet
P1 = 357  # Point load 1 in lbs
a1 = 2.17  # Distance of point load 1 from the left support in feet
P2 = 234  # Point load 2 in lbs
a2 = 8.5 - 2.82  # Distance of point load 2 from the left support in feet
w_psf = 75  # Uniform load in lbs/ft^2
F_w = 0  # Total horizontal wind force in lbs


##### Dimensions ########
b = 6 / 12  # Width of the beam in feet
h = 6 / 12  # Height of the beam in feet
thickness = 1  # Thickness of the beam in inches
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
    R1 = (P1 * (L - a1) + P2 * (L - a2) + w_plf * L * (L / 2)) / L
    R2 = P1 + P2 + w_plf * L - R1
    return R1, R2

def calculate_shear_force(x, L, P1, a1, P2, a2, w_plf, R1):
    # shear force at a distance x from the left supp
    V = 0
    if x < a1:
        V = R1 - w_plf * x
    elif x < a2:
        V = R1 - P1 - w_plf * x
    else:
        V = R1 - P1 - P2 - w_plf * x
    return V

def calculate_bending_moment(x, L, P1, a1, P2, a2, w_plf, R1):
    # bending moment at a distance x from the left supp
    M = 0
    if x < a1:
        M = R1 * x - (w_plf * x**2) / 2
    elif x < a2:
        M = R1 * x - P1 * (x - a1) - (w_plf * x**2) / 2
    else:
        M = R1 * x - P1 * (x - a1) - P2 * (x - a2) - (w_plf * x**2) / 2
    return M

def calculate_bending_stress(M, y, I):
    # Calculate bending stress
    return M * y / I

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
    ### From Beam Equations Word doc
    tau = (V*calc_Q(y1) / (calc_Ic() * thickness)) * (binch*(hinch**2 - hw**2) + thickness*(hw**2 - (4*(y1**2))))
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

def calc_horizontal_shear(F_w, L):
    # Calculate horizontal shear force due to wind
    return F_w / L

R1, R2 = calculate_reactions(L, P1, a1, P2, a2, w_plf)
x_values = np.linspace(0, L, num=100)
shear_forces = [calculate_shear_force(x, L, P1, a1, P2, a2, w_plf, R1) for x in x_values]
shear_max_stresses = [calc_tau(V, 0) for V in shear_forces]
moments = [calculate_bending_moment(x, L, P1, a1, P2, a2, w_plf, R1) for x in x_values]
bending_stresses = [calculate_bending_stress(M, y, I) for M in moments]
strains = [calculate_strain(sigma_bending, E) for sigma_bending in bending_stresses]

horizontal_shear = calc_horizontal_shear(F_w, L)
horizontal_max_stress = calc_tau_hor(horizontal_shear, 0)

# combine vert and hor stresses
resultant_shear_stresses = [np.sqrt(v**2 + horizontal_max_stress**2) for v in shear_max_stresses]

plt.figure(figsize=(14,10))

plt.subplot(5, 1, 1)
plt.plot(x_values, shear_forces)
plt.title('Shear Force along the Beam')
plt.xlabel('Position along the beam (ft)')
plt.ylabel('Shear Force (lbs)')
plt.grid(True)

plt.subplot(5, 1, 2)
plt.plot(x_values, shear_max_stresses)
plt.title('Shear Max Stress along the Beam')
plt.xlabel('Position along the beam (ft)')
plt.ylabel('Shear Max Stress (psi)')
plt.grid(True)

plt.subplot(5, 1, 3)
plt.plot(x_values, moments)
plt.title('Bending Moment along the Beam')
plt.xlabel('Position along the beam (ft)')
plt.ylabel('Bending Moment (lb-ft)')
plt.grid(True)

plt.subplot(5, 1, 4)
plt.plot(x_values, bending_stresses)
plt.title('Bending Stress along the Beam')
plt.xlabel('Position along the beam (ft)')
plt.ylabel('Bending Stress (psi)')
plt.grid(True)

plt.subplot(5, 1, 5)
plt.plot(x_values, resultant_shear_stresses)
plt.title('Resultant Shear Stress along the Beam')
plt.xlabel('Position along the beam (ft)')
plt.ylabel('Resultant Shear Stress (psi)')
plt.grid(True)

plt.tight_layout()
plt.savefig('beam_analysis.png')
plt.show()

data = {
    'Position along the beam (ft)': x_values,
    'Shear Force (lbs)': shear_forces,
    'Shear Max Stress (psi)': shear_max_stresses,
    'Bending Moment (lb-ft)': moments,
    'Bending Stress (psi)': bending_stresses,
    'Horizontal Shear Stress (psi)': horizontal_max_stress,
    'Resultant Shear Stress (psi)': resultant_shear_stresses,
    'Strains': strains,
}
df = pd.DataFrame(data)
df.to_excel('beam_analysis_results_with_wind.xlsx', index=False)
print("Results have been saved to 'beam_analysis_results_with_wind.xlsx'.")
