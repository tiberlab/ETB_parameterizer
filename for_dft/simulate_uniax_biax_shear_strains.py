import numpy as np
import pickle
import os
import itertools
import sys


def create_energy_range(bandstructure, k_point, FCC_OR_WZ=0):
    bands_array = bandstructure.evaluate().inUnitsOf(eV)
    n_kpoints = bands_array.shape[0]

    if FCC_OR_WZ == 0:
        if k_point == "L":
            k_index = 0
        elif k_point == "G":
            k_index = n_kpoints//2
        elif k_point == "X":
            k_index = -1
        else:
            raise ValueError("Argument k_point must be a string chosen from: ['L', 'G', 'X'].")
    if FCC_OR_WZ == 1:
        if k_point == "G":
            k_index = n_kpoints//2
        if k_point == "Y":
            k_index = 0
        if k_point == "Z":
            k_index = -1


    energies = bands_array[k_index,:]

    return energies* eV


def setup_bulk_FCC(cation, strain_tensor):
    if cation == "Ga":
        elements = [Gallium, Nitrogen]

    elif cation == "Al":
        elements = [Aluminium, Nitrogen]

    elif cation == "In":
        elements = [Indium, Nitrogen]
    
    else:
        raise ValueError("Please insert one of the three possible atoms as cation:\n- Ga\n- Al\n- In")
    

    # Initial GaN lattice. Corresponds to a lattice constant of 4.531
    vector_a = np.array([0.0, 2.2655, 2.2655]).reshape(-1,1) #*Angstrom
    vector_b = np.array([2.2655, 0.0, 2.2655]).reshape(-1,1) #*Angstrom
    vector_c = np.array([2.2655, 2.2655, 0.0]).reshape(-1,1) #*Angstrom

    # Apply strain tensor
    vector_a = (vector_a + np.matmul(strain_tensor, vector_a)).reshape(-1)*Angstrom
    vector_b = (vector_b + np.matmul(strain_tensor, vector_b)).reshape(-1)*Angstrom
    vector_c = (vector_c + np.matmul(strain_tensor, vector_c)).reshape(-1)*Angstrom
    
    # Set up lattice
    lattice = UnitCell(vector_a, vector_b, vector_c)

    # Define coordinates
    fractional_coordinates = [[0., 0., 0.],
                              [0.2500, 0.2500, 0.2500]]

    # Set up configuration
    bulk = BulkConfiguration(
        bravais_lattice=lattice,
        elements=elements,
        fractional_coordinates=fractional_coordinates
    )

    return bulk


def setup_bulk_WZ(cation, strain_tensor):
    if cation == "Ga":
       elements = [Gallium, Gallium, Nitrogen, Nitrogen]

    elif cation == "Al":
        elements = [Aluminium, Aluminium, Nitrogen, Nitrogen]

    elif cation == "In":
        elements = [Indium, Indium, Nitrogen, Nitrogen]
    
    else:
        raise ValueError("Please insert one of the three possible atoms as cation:\n- Ga\n- Al\n- In")
    

    # Initial GaN lattice. 
    vector_a = np.array([1.5945, -2.76176, 0]).reshape(-1,1) #*Angstrom
    vector_b = np.array([1.5945, 2.76176, 0]).reshape(-1,1) #*Angstrom
    vector_c = np.array([0, 0, 5.185]).reshape(-1,1) #*Angstrom

    # Apply strain tensor
    vector_a = (vector_a + np.matmul(strain_tensor, vector_a)).reshape(-1)*Angstrom
    vector_b = (vector_b + np.matmul(strain_tensor, vector_b)).reshape(-1)*Angstrom
    vector_c = (vector_c + np.matmul(strain_tensor, vector_c)).reshape(-1)*Angstrom
    
    # Set up lattice
    lattice = UnitCell(vector_a, vector_b, vector_c)

    # Define coordinates
    fractional_coordinates = [[ 0.333333333333,  0.666666666667,  0.            ],
                            [ 0.666666666667,  0.333333333333,  0.5           ],
                            [ 0.333333333333,  0.666666666667,  0.385         ],
                            [ 0.666666666667,  0.333333333333,  0.885         ]]

    # Set up configuration
    bulk = BulkConfiguration(
        bravais_lattice=lattice,
        elements=elements,
        fractional_coordinates=fractional_coordinates
    )

    return bulk


def run_simulation(alpha, configuration, cation, file_string):
    #----------------------------------------
    # Exchange-Correlation
    #----------------------------------------
    exchange_correlation = HSECustomExchangeCorrelation(
        screening_length=0.11*1/Bohr,
        exx_fraction=alpha,
        number_of_spins=4,
        spin_orbit=True)

    k_point_sampling = KpointDensity(
        density_a=4.0*Angstrom,
        density_b=4.0*Angstrom,
        density_c=4.0*Angstrom
    )

    numerical_accuracy_parameters = NumericalAccuracyParameters(
        k_point_sampling=k_point_sampling
    )

    checkpoint_handler = NoCheckpointHandler

    calculator = PlaneWaveCalculator(
        exchange_correlation=exchange_correlation,
        numerical_accuracy_parameters=numerical_accuracy_parameters,
        checkpoint_handler=checkpoint_handler,
        store_wave_functions=None
    )

    configuration.setCalculator(calculator)

    configuration.update()

    # Save also configuration
    filename = f"{file_string}/{cation}N_PW_configuration_{file_string}"
    nlsave(filename + ".hdf5", configuration)

    return configuration

def compute_and_save_bandstructure_FCC(configuration, cation, file_string):

    # Band structure computation
    bandstructure = Bandstructure(
            configuration=configuration,
            route=['L', 'G', 'B'],  # X point is called B in the UnitCell bulk structure
            points_per_segment=1
        )
    
    filename = f"{file_string}/{cation}N_PW_bands_{file_string}"
    nlsave(filename + ".hdf5", bandstructure)
  
    return bandstructure


def compute_and_save_bandstructure_WZ(configuration, cation, file_string):

    # Band structure computation
    bandstructure = Bandstructure(
            configuration=configuration,
            route=['Y', 'G', 'Z'],  # M and A points are called Y and Z in the UnitCell bulk structure
            points_per_segment=1
        )
    
    filename = f"{file_string}/{cation}N_PW_bands_{file_string}"
    nlsave(filename + ".hdf5", bandstructure)

    return bandstructure


def compute_and_save_pdos_FCC(configuration, bandstructure, cation, file_string):

    # PDOS parameters
    spectrum_method = GaussianBroadening(broadening=0.01*eV)
    G = RegularKpointGrid(ka_range=0.0, kb_range=0.0, kc_range=0.0, na=1, nb=1, nc=1)
    L = RegularKpointGrid(ka_range=0.5, kb_range=0.5, kc_range=0.5, na=1, nb=1, nc=1)
    X = RegularKpointGrid(ka_range=0.5, kb_range=0.0, kc_range=0.5, na=1, nb=1, nc=1)
    energies_G = create_energy_range(bandstructure, "G")
    energies_L = create_energy_range(bandstructure, "L")
    energies_X = create_energy_range(bandstructure, "X")

    # PDOS for Gamma point
    pdos = ProjectedDensityOfStates(configuration, kpoints=G, projections=ProjectOnOrbitalsBySite, 
                                    energies=energies_G, energy_zero_parameter=FermiLevel,
                                    spectrum_method=spectrum_method)
    
    filename = f"{file_string}/{cation}N_PW_PDOS_{file_string}"
    nlsave(filename + ".hdf5", pdos, object_id="G")

    # PDOS for L point
    pdos = ProjectedDensityOfStates(configuration, kpoints=L, projections=ProjectOnOrbitalsBySite, 
                                    energies=energies_L, energy_zero_parameter=FermiLevel,
                                    spectrum_method=spectrum_method)
    
    nlsave(filename + ".hdf5", pdos, object_id="L")

    # PDOS for X point
    pdos = ProjectedDensityOfStates(configuration, kpoints=X, projections=ProjectOnOrbitalsBySite, 
                                    energies=energies_X, energy_zero_parameter=FermiLevel,
                                    spectrum_method=spectrum_method)
    
    nlsave(filename + ".hdf5", pdos, object_id="X")


def compute_and_save_pdos_WZ(configuration, bandstructure, cation, file_string):

    # PDOS parameters
    spectrum_method = GaussianBroadening(broadening=0.001*eV)
    G = RegularKpointGrid(ka_range=0.0, kb_range=0.0, kc_range=0.0, na=1, nb=1, nc=1)
    Y = RegularKpointGrid(ka_range=0.0, kb_range=0.5, kc_range=0.0, na=1, nb=1, nc=1)
    Z = RegularKpointGrid(ka_range=0.0, kb_range=0.0, kc_range=0.5, na=1, nb=1, nc=1)
    
    energies_G = create_energy_range(bandstructure, "G", FCC_OR_WZ=1)
    energies_Z = create_energy_range(bandstructure, "Z", FCC_OR_WZ=1)
    energies_Y = create_energy_range(bandstructure, "Y", FCC_OR_WZ=1)
        
    # PDOS for Gamma point
    pdos = ProjectedDensityOfStates(configuration, kpoints=G, projections=ProjectOnOrbitalsBySite, 
                                    energies=energies_G, energy_zero_parameter=FermiLevel,
                                    spectrum_method=spectrum_method)
    
    filename = f"{file_string}/{cation}N_PW_PDOS_{file_string}"
    nlsave(filename + ".hdf5", pdos, object_id="G")


    # PDOS for Z point
    pdos = ProjectedDensityOfStates(configuration, kpoints=Z, projections=ProjectOnOrbitalsBySite, 
                                    energies=energies_Z, energy_zero_parameter=FermiLevel,
                                    spectrum_method=spectrum_method)
    
    filename = f"{file_string}/{cation}N_PW_PDOS_{file_string}"
    nlsave(filename + ".hdf5", pdos, object_id="Z")


    # PDOS for Y point
    pdos = ProjectedDensityOfStates(configuration, kpoints=Y, projections=ProjectOnOrbitalsBySite, 
                                    energies=energies_Y, energy_zero_parameter=FermiLevel,
                                    spectrum_method=spectrum_method)
    
    filename = f"{file_string}/{cation}N_PW_PDOS_{file_string}"
    nlsave(filename + ".hdf5", pdos, object_id="Y")




def simulate_and_save_objects_FCC(optimal_alpha, bulk, cation, strain_type, eps_values):

    if strain_type == "uniax":
        file_string = f"{strain_type}_exx_{eps_values}"
    
    if strain_type == "biax":
        file_string = f"{strain_type}_exx_{eps_values[0]}_eyy_{eps_values[1]}"

    if strain_type == "shear":
        file_string = f"{strain_type}_exy_{eps_values}"

    os.makedirs(f"{file_string}", exist_ok = True)

    configuration = run_simulation(optimal_alpha, bulk, cation, file_string)

    bandstructure = compute_and_save_bandstructure_FCC(configuration, cation, file_string)

    compute_and_save_pdos_FCC(configuration, bandstructure, cation, file_string)


def simulate_and_save_objects_WZ(optimal_alpha, bulk, cation, strain_type, eps_values):

    if strain_type == "uniax":
        file_string = f"{strain_type}_exx_{eps_values}"
    
    if strain_type == "biax":
        file_string = f"{strain_type}_exx_{eps_values[0]}_eyy_{eps_values[1]}"

    if strain_type == "shear":
        file_string = f"{strain_type}_exy_{eps_values}"

    os.makedirs(f"{file_string}", exist_ok = True)

    configuration = run_simulation(optimal_alpha, bulk, cation, file_string)

    bandstructure = compute_and_save_bandstructure_WZ(configuration, cation, file_string)

    compute_and_save_pdos_WZ(configuration, bandstructure, cation, file_string)
    

def main(cation, alpha, FCC_OR_WZ):
    
    # Uniaxial strain
    eps_xx = [-0.04, -0.02, 0.02, 0.04]

    if isMainProcess():
        print(f"\nUniaxial strain: {cation}N, "
              f"eps_xx = {eps_xx}.\n")

    for eps in eps_xx:
        strain_tensor = np.zeros((3,3))
        strain_tensor[0,0] = eps
	
        if FCC_OR_WZ == 0:
            bulk = setup_bulk_FCC(cation, strain_tensor)
            simulate_and_save_objects_FCC(alpha, bulk, cation, "uniax", eps)
        elif FCC_OR_WZ == 1:
            bulk = setup_bulk_WZ(cation, strain_tensor)
            simulate_and_save_objects_WZ(alpha, bulk, cation, "uniax", eps)


    # Biaxial strain
    eps_vals = [-0.04, -0.02, 0.02, 0.04]
    
    if isMainProcess():
        print(f"\nBiaxial strain: {cation}N, "
            f"(eps_xx, eps_yy) =  {list(itertools.combinations_with_replacement(eps_vals, 2))}.\n")
	
    for eps_couple in itertools.combinations_with_replacement(eps_vals, 2):
        eps_xx = eps_couple[0]
        eps_yy = eps_couple[1]

        strain_tensor = np.zeros((3,3))
        strain_tensor[0,0] = eps_xx
        strain_tensor[1,1] = eps_yy

        if FCC_OR_WZ == 0:
            bulk = setup_bulk_FCC(cation, strain_tensor)
            simulate_and_save_objects_FCC(alpha, bulk, cation, "biax", eps_couple)
        elif FCC_OR_WZ == 1:
            bulk = setup_bulk_WZ(cation, strain_tensor)
            simulate_and_save_objects_WZ(alpha, bulk, cation, "biax", eps_couple)

	
    # Shear strain
    eps_xy = [-0.04, -0.02, 0.02, 0.04]
    
    if isMainProcess():
        print(f"\nShear strain: {cation}N, "
            f"eps_xy = {eps_xy}.\n\n")
	
    for eps in eps_xy:
        strain_tensor = np.zeros((3,3))
        strain_tensor[0,1] = eps
        strain_tensor[1,0] = eps

        if FCC_OR_WZ == 0:
            bulk = setup_bulk_FCC(cation, strain_tensor)
            simulate_and_save_objects_FCC(alpha, bulk, cation, "shear", eps)
        elif FCC_OR_WZ == 1:
            bulk = setup_bulk_WZ(cation, strain_tensor)
            simulate_and_save_objects_WZ(alpha, bulk, cation, "shear", eps)

	

if __name__ == '__main__':
    
    # args = sys.argv[1:]
    # args = ["Ga", 0.39, 4.531]
    args = ["Ga", 0.399]

    FCC_OR_WZ = 1  # 0 for FCC, 1 for WZ    

    if len(args) < 2:
        raise ValueError("Please insert the follwing quantities:\n1) Cation element (Ga, Al, In);\n2) Alpha mixing parameter.")

    cation, alpha = args

    cation = str(cation)
    try:
        alpha = float(alpha)
    except:
        raise ValueError("Please provide a float value for the alpha parameter.")

    main(cation, alpha, FCC_OR_WZ)







