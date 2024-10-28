import numpy as np

# Currently only tests for load bearing and flexing, the included materials need to be adjusted to include all types
# of steel and other metals. Need to clarify on the "Movement of the stacks" with kolla.

######srryforbeinglatetothepartyguys###################################################################################
#################################inputs_into_system################written by andrewpace###############################
#######################################################################################################################

stack_displacement = 2.625  # inches
stack_frequency = 1 / (13.0*60.0)  # Hz
wind_force = 1615.0  # pounds-force
live_load = 75.0  # pounds-force per square foot (industry standard live load) [^1^][2]
beam_length = (9.0*12.0) + 1.5  # inches (length of the beam)

#limits for beam dimensions in imperial units
min_width = 6.0  # inches
max_width = 8.0  # inches
min_height = 6.0  # inches
max_height = 12.0  # inches
min_diameter = 6.0
max_diameter = 8.0
#safety factors
safety_factor_stress = 1.5
safety_factor_deflection = 1.5
thickness_for_hollow_beams_and_ibeam = 0.1 #multiplier for dimension, so 0.1 is 10% of the width/height whichever is lower

#additional point loads and their positions
point_load1 = 357  # pounds-force
position1 = (3.0*12.0) + 11 + (1/16)  # inches from the left end
point_load2 = 234  # pounds-force
position2 = ((3.0*12.0) + 11 + (1/16)) + ((3.0*12.0) + 6 + (1/8))  # inches from the left end

# basic materials, needs to be filled in with actual material properties, these are currently just generic props
materials = {
                # psi, psi, lb/in^3, $/in^3
    'steel': {'E': 29e6, 'yield_strength': 36e3, 'density': 0.284, 'cost_per_cubic_inch': 0.2},
    'aluminum': {'E': 10e6, 'yield_strength': 35e3, 'density': 0.098, 'cost_per_cubic_inch': 0.8},
    'titanium': {'E': 16.5e6, 'yield_strength': 120e3, 'density': 0.163, 'cost_per_cubic_inch': 10.0},
    'copper': {'E': 17e6, 'yield_strength': 33e3, 'density': 0.321, 'cost_per_cubic_inch': 2.5},
    'beryllium_copper': {'E': 19e6, 'yield_strength': 160e3, 'density': 0.302, 'cost_per_cubic_inch': 1.0},
    'junktestmaterial': {'E': 1e4, 'yield_strength': 21e3, 'density': 0.21, 'cost_per_cubic_inch': 0.3},
    'aisi_4140': {'E': 2.97e7, 'yield_strength': 9.5e4, 'density': 0.284, 'cost_per_cubic_inch': 0.5},
    'aisi_4340': {'E': 3.05e7, 'yield_strength': 1.25e5, 'density': 0.284, 'cost_per_cubic_inch': 0.6},
    's355': {'E': 3.05e7, 'yield_strength': 5.15e4, 'density': 0.284, 'cost_per_cubic_inch': 0.4}

}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

# calculate natural frequency of the beam
def natural_frequency(E, I, L, m):
    return (1 / (2 * np.pi)) * np.sqrt((E * I) / (m * L ** 3))

# calculate stress in the beam
def calculate_stress(F, L, I, y):
    return (F * L * y) / I

# calculate deflection of the beam
def calculate_deflection(F, L, E, I):
    return (F * L ** 3) / (48 * E * I)

# calculate deflection due to point loads
def calculate_point_load_deflection(P, a, L, E, I):
    b = L - a
    return (P * a * b * (L ** 2 - a ** 2 - b ** 2)) / (6 * E * I * L)


# Function to calculate cost pof beam
def calculate_cost(density, volume, cost_per_cubic_inch):
    return volume * cost_per_cubic_inch  # $


# Function to calculate properties for different beam shapes
def calculate_properties(shape, width, height, diameter, thickness):
    if shape == 'rectangular':
        I = (width * height ** 3) / 12
        y = height / 2
        volume = width * height * beam_length  # in^3
    elif shape == 'hollow_rectangular':
        I = (width * height ** 3 - (width - 2 * thickness) * (height - 2 * thickness) ** 3) / 12
        y = height / 2
        volume = (width * height - (width - 2 * thickness) * (height - 2 * thickness)) * beam_length  # in^3
    elif shape == 'I-beam':
        I = (width * height ** 3) / 12 - ((width - thickness) * (height - 2 * thickness) ** 3) / 12
        y = height / 2
        volume = (width * height - (width - thickness) * (height - 2 * thickness)) * beam_length  # in^3
    elif shape == 'round_solid':
        I = (np.pi * diameter ** 4) / 64
        y = diameter / 2
        volume = (np.pi * diameter ** 2 / 4) * beam_length  # in^3
    elif shape == 'round_hollow':
        I = (np.pi * (diameter ** 4 - (diameter - 2 * thickness) ** 4)) / 64
        y = diameter / 2
        volume = (np.pi * (diameter ** 2 - (diameter - 2 * thickness) ** 2) / 4) * beam_length  # in^3
    return I, y, volume

def can_flex_enough_bool(properties, E_calced, stress_calced):
    print(f'E_good?: {properties['E']< E_calced}')
    print(f'bending_stress_ok?: {(properties['yield_strength'] * safety_factor_stress) > stress_calced}')
    if properties['E']< E_calced and properties['yield_strength'] > stress_calced:
        return True
    else:
        return False



def bending_moment_from_deflection(delta, E, I, L):
    """
    Calculate the bending moment from deflection and material properties.

    Parameters:
    delta (float): Deflection in inches
    E (float): Young's modulus in psi
    I (float): Moment of inertia in in^4
    L (float): Length of the beam in inches

    Returns:
    float: Bending moment in lb-in
    """
    return (8 * delta * E * I) / (L ** 2)


def calculate_youngs_modulus_and_stress(moment_of_inertia, height, E_material, length= beam_length,
                                        deflection = stack_displacement, load=wind_force*1000):
    """
    Calculate Young's modulus (E) and bending stress (σ) for a simply supported beam with a central load.

    Parameters:
    moment_of_inertia (float): Moment of inertia of the beam's cross-section.
    length (float): Length of the beam.
    deflection (float): Deflection of the beam under the load.
    load (float): Applied load on the beam.
    height (float): Height of the beam's cross-section.

    Returns:
    tuple: Calculated Young's modulus (E) and bending stress (σ).
    """
    # Calculate Young's modulus (E) needed to allow bending
    E = (load * length ** 3) / (48 * deflection * moment_of_inertia)

    # bending moment (M) at the center of the beam
    M = bending_moment_from_deflection(deflection, E_material, moment_of_inertia, length)

    # distance from the neutral axis to the outermost fiber (c)
    c = height / 2

    # bending stress (σ)
    stress = M * c / moment_of_inertia

    return E, stress

def calculate_bending_stress(E, I_in4, y_in, load_lbf, load_position_in): # load_pos is inches from the left most support
    # Calculate reactions at supports (assuming simply supported beam)
    R1 = load_lbf * (beam_length - load_position_in) / beam_length
    R2 = load_lbf * load_position_in / beam_length

    # Calculate maximum bending moment
    if load_position_in <= beam_length / 2:
        M_max = R1 * load_position_in
    else:
        M_max = R2 * (beam_length - load_position_in)

    # Calculate maximum stress (at the top/bottom surface of the beam)
    stress_max = M_max * y_in / I_in4
    return stress_max #returns psi stress

def calculate_bending_stress_live_load(I, y, live_load_psf):
    # w is the load per unit length (lb/in)
    w = live_load_psf * (1 / 12)  # Convert psf to lb/in
    M = (w * beam_length ** 2) / 8
    # Calculate the bending stress (σ)
    sigma = (M * y) / I

    return sigma


def calculate_bending_stress_from_deflection(width_in, height_in, E_psi,
                                             length_in=2*beam_length, deflection_in=stack_displacement):
        # """Calculate bending stress from the movement of the stacks.
        # :param width_in:
        # :param height_in:
        # :param E_psi:
        # :param length_in: used 2 * beam_length, because these calculation are for the middle of the beam, but the beam
        # is deflecting across the entire beam length and not just down for half and up for half. mught be incorrect but it
        # is giving numbers that are closer to what is expected
        # :param deflection_in:
        # :return:
        # """

    # inverted I for horizontal bendmoment
    I = (height_in * width_in ** 3) / 12

    # Calculate the load (P) from deflection
    P = (deflection_in * 48 * E_psi * I) / (length_in ** 3)

    # Calculate the bending moment (M) at the center
    M = (P * length_in) / 4

    # Calculate the distance from the neutral axis to the outermost fiber (c)
    c = width_in / 2

    # Calculate the bending stress (σ)
    sigma = (M * c) / I

    return sigma

def convert_lbf_to_psf(force_lbf, area_ft2):
    # Calculate pressure in psf
    pressure_psf = force_lbf / area_ft2
    return pressure_psf


results = []

shapes = ['rectangular', 'hollow_rectangular', 'I-beam', 'round_solid', 'round_hollow']

for material, properties in materials.items():
    E = properties['E']
    yield_strength = properties['yield_strength']
    density = properties['density']
    cost_per_cubic_inch = properties['cost_per_cubic_inch']

    for shape in shapes:
        best_design = None
        min_dimensions = float('inf')

        biggest_design = None
        max_dimensions = -float('inf')

        if shape in ['rectangular', 'hollow_rectangular', 'I-beam']:
            for width in np.linspace(min_width, max_width, 10):  # inches
                for height in np.linspace(min_height, max_height, 10):  # inches
                    thickness = thickness_for_hollow_beams_and_ibeam * min(width, height)  # thicknes for hollow beams
                    I, y, volume = calculate_properties(shape, width, height, None, thickness)
                    m = density * volume / 1728  # converting in^3 to ft^3
                    E_needed_for_shape_n_dim = calculate_youngs_modulus_and_stress(I, height, properties['E'])
                    print(material)
                    print(f'For {shape}, width: {width}in, height: {height}in')
                    # print(f'E_needs_to_be_under: {E_needed_for_shape_n_dim[0]}psi')
                    # print(f'stress_on_the_beam_from_flexing: {E_needed_for_shape_n_dim[1]}psi')
                    # print(f'I: {I}')

                    freq = natural_frequency(E, I, beam_length, m)
                    if freq < stack_frequency:
                        continue


                    # New
                    point_load_bending_stress = (calculate_bending_stress(E, I, y, point_load1, position1) +
                                                 calculate_bending_stress(E, I, y, point_load2, position2))

                    live_load_bending_stress = calculate_bending_stress_live_load(I, y, live_load)
                    stack_movement_stress = calculate_bending_stress_from_deflection(width, height, E)
                    wind_converted_psf = convert_lbf_to_psf(wind_force, width*height)
                    wind_stress = calculate_bending_stress_live_load(I, width/2, live_load)
                    #endnew


                    # wind force
                    deflection_horizontal = calculate_deflection(wind_force, beam_length, E, I)

                    # live load and point loads
                    total_vertical_force = live_load * width * beam_length / 144  # Adding live load
                    deflection_vertical = calculate_deflection(total_vertical_force, beam_length, E, I)
                    deflection_vertical += calculate_point_load_deflection(point_load1, position1, beam_length, E, I)
                    deflection_vertical += calculate_point_load_deflection(point_load2, position2, beam_length, E, I)

                    # vertical loads
                    total_force = total_vertical_force + point_load1 + point_load2
                    stress = point_load_bending_stress + live_load_bending_stress + stack_movement_stress + wind_stress
                    it_can_flex_enough = can_flex_enough_bool(properties, E_needed_for_shape_n_dim[0], E_needed_for_shape_n_dim[1])
                    # print(it_can_flex_enough)

                    allowable_stress = yield_strength / safety_factor_stress
                    allowable_deflection_for_load = 1 #inches
                    print(f'Total Force: {total_force} |||||| Deflection: {deflection_vertical}')
                    print(f'Total Stress: {stress}  |||||| Allowable Stress: {allowable_stress}')
                    print(f'stresstest: {stress < allowable_stress}, flextest: {it_can_flex_enough}')
                    print('')




                    if (stress < allowable_stress):
                        dimensions = width * height
                        if dimensions < min_dimensions:
                            min_dimensions = dimensions
                            cost = calculate_cost(density, volume, cost_per_cubic_inch)
                            best_design = (material, shape, width, height, cost, it_can_flex_enough)
                        if dimensions > max_dimensions:
                            max_dimensions = dimensions
                            cost = calculate_cost(density, volume, cost_per_cubic_inch)
                            biggest_design = (material, shape, width, height, cost, it_can_flex_enough)

        elif shape in ['round_solid', 'round_hollow']:
            for diameter in np.linspace(min_diameter, max_diameter, 10):  # inches
                thickness = thickness_for_hollow_beams_and_ibeam * diameter  # thicikness for hollow beams
                I, y, volume = calculate_properties(shape, None, None, diameter, thickness)
                m = density * volume / 1728  # converting in^3 to ft^3

                E_needed_for_shape_n_dim = calculate_youngs_modulus_and_stress(I, diameter, properties['E'])
                print(material)
                print(f'For {shape}, diameter {diameter}in')
                print(f'E_needs_to_be_under: {E_needed_for_shape_n_dim[0]}psi')
                print(f'stress_on_the_beam_from_flexing: {E_needed_for_shape_n_dim[1]}psi')

                freq = natural_frequency(E, I, beam_length, m)
                if freq < stack_frequency:
                    continue

                # wind force
                deflection_horizontal = calculate_deflection(wind_force, beam_length, E, I)

                # live load and point loads
                total_vertical_force = live_load * (np.pi * diameter ** 2 / 4) * beam_length / 144  # Adding live load
                deflection_vertical = calculate_deflection(total_vertical_force, beam_length, E, I)
                deflection_vertical += calculate_point_load_deflection(point_load1, position1, beam_length, E, I)
                deflection_vertical += calculate_point_load_deflection(point_load2, position2, beam_length, E, I)

                # vertical loads
                total_force = total_vertical_force + point_load1 + point_load2
                stress = calculate_stress(total_force, beam_length, I, y)
                it_can_flex_enough = can_flex_enough_bool( properties,E_needed_for_shape_n_dim[0], E_needed_for_shape_n_dim[1])

                allowable_stress = yield_strength / safety_factor_stress
                horizontal_maximum_wind_deflection = 1
                vertical_maximum_deflection = 1
                print(
                    f'stresstest: {stress < allowable_stress}, deflectHtest: {deflection_horizontal < horizontal_maximum_wind_deflection}, deflectVtest: {deflection_vertical < vertical_maximum_deflection}, flextest: {it_can_flex_enough}')
                print('')

                if (stress < allowable_stress and deflection_horizontal < horizontal_maximum_wind_deflection
                        and deflection_vertical < vertical_maximum_deflection and it_can_flex_enough):
                    dimensions = diameter
                    if dimensions < min_dimensions:
                        min_dimensions = dimensions
                        cost = calculate_cost(density, volume, cost_per_cubic_inch)
                        best_design = (material, shape, diameter, None, cost, it_can_flex_enough)

        if best_design:
            results.append(best_design)
        if biggest_design:
            results.append(biggest_design)


# prints "best design" for each shape and material if meets load and flex criteria
for result in results:
    material, shape, width, height, cost, flex_bool = result
    if shape in ['rectangular', 'hollow_rectangular', 'I-beam']:
        print(f"Material: {material}, Shape: {shape}, Width: {width:.2f}in, Height: {height:.2f}in, Cost: ${cost:.2f}")
    else:
        print(f"Material: {material}, Shape: {shape}, Diameter: {width:.2f}in, Cost: ${cost:.2f}")

