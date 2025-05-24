import numpy as np
import pickle
import os
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
        if k_point == "M":
            k_index = 0
        if k_point == "A":
            k_index = -1


    energies = bands_array[k_index,:]

    return energies* eV


def setup_bulk_FCC(cation, lattice_constant):
    if cation == "Ga":
        elements = [Gallium, Nitrogen]

    elif cation == "Al":
        elements = [Aluminium, Nitrogen]

    elif cation == "In":
        elements = [Indium, Nitrogen]
    
    else:
        raise ValueError("Please insert one of the three possible atoms as cation:\n- Ga\n- Al\n- In")

    # Set up lattice
    lattice = FaceCenteredCubic(lattice_constant * Angstrom)

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


def setup_bulk_WZ(cation, a, c):
    if cation == "Ga":
       elements = [Gallium, Gallium, Nitrogen, Nitrogen]

    elif cation == "Al":
        elements = [Aluminium, Aluminium, Nitrogen, Nitrogen]

    elif cation == "In":
        elements = [Indium, Indium, Nitrogen, Nitrogen]
    
    else:
        raise ValueError("Please insert one of the three possible atoms as cation:\n- Ga\n- Al\n- In")

    # Set up lattice
    lattice = Hexagonal(a*Angstrom, c*Angstrom)  

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


def run_simulation(alpha, configuration, lc_string, cation, strain_string):
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
    filename = f"{lc_string}_{strain_string}/{cation}N_PW_configuration_{lc_string}_{strain_string}"
    nlsave(filename + ".hdf5", configuration)

    return configuration

def compute_and_save_bandstructure_FCC(configuration, lc_string, cation, strain_string):

    # Band structure computation
    bandstructure = Bandstructure(
            configuration=configuration,
            route=['L', 'G', 'X'],
            points_per_segment=20
        )
    
    filename = f"{lc_string}_{strain_string}/{cation}N_PW_bands_{lc_string}_{strain_string}"
    nlsave(filename + ".hdf5", bandstructure)

    # with open(filename + ".pickle", "wb") as f:
    #    pickle.dump(bandstructure.evaluate().inUnitsOf(eV), f)
    
    return bandstructure


def compute_and_save_bandstructure_WZ(configuration, lc_string, cation, strain_string):

    # Band structure computation
    bandstructure = Bandstructure(
            configuration=configuration,
            route=['M', 'G', 'A'],
            points_per_segment=20
        )
    
    filename = f"{lc_string}_{strain_string}/{cation}N_PW_bands_{lc_string}_{strain_string}"
    nlsave(filename + ".hdf5", bandstructure)

    # with open(filename + ".pickle", "wb") as f:
    #    pickle.dump(bandstructure.evaluate().inUnitsOf(eV), f)
    
    return bandstructure


def compute_and_save_pdos_FCC(configuration, bandstructure, lc_string, cation, strain_string):

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
    
    filename = f"{lc_string}_{strain_string}/{cation}N_PW_PDOS_{lc_string}_{strain_string}"
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


def compute_and_save_pdos_WZ(configuration, bandstructure, lc_string, cation, strain_string):

    # PDOS parameters
    spectrum_method = GaussianBroadening(broadening=0.001*eV)
    G = RegularKpointGrid(ka_range=0.0, kb_range=0.0, kc_range=0.0, na=1, nb=1, nc=1)
    M = RegularKpointGrid(ka_range=0.0, kb_range=0.5, kc_range=0.0, na=1, nb=1, nc=1)
    A = RegularKpointGrid(ka_range=0.0, kb_range=0.0, kc_range=0.5, na=1, nb=1, nc=1)
    
    energies_G = create_energy_range(bandstructure, "G", FCC_OR_WZ=1)
    energies_A = create_energy_range(bandstructure, "A", FCC_OR_WZ=1)
    energies_M = create_energy_range(bandstructure, "M", FCC_OR_WZ=1)
        
    # PDOS for Gamma point
    pdos = ProjectedDensityOfStates(configuration, kpoints=G, projections=ProjectOnOrbitalsBySite, 
                                    energies=energies_G, energy_zero_parameter=FermiLevel,
                                    spectrum_method=spectrum_method)
    
    filename = f"{lc_string}_{strain_string}/{cation}N_PW_PDOS_{lc_string}_{strain_string}"
    nlsave(filename + ".hdf5", pdos, object_id="G")


    # PDOS for A point
    pdos = ProjectedDensityOfStates(configuration, kpoints=A, projections=ProjectOnOrbitalsBySite, 
                                    energies=energies_A, energy_zero_parameter=FermiLevel,
                                    spectrum_method=spectrum_method)
    
    filename = f"{lc_string}_{strain_string}/{cation}N_PW_PDOS_{lc_string}_{strain_string}"
    nlsave(filename + ".hdf5", pdos, object_id="A")


    # PDOS for M point
    pdos = ProjectedDensityOfStates(configuration, kpoints=M, projections=ProjectOnOrbitalsBySite, 
                                    energies=energies_M, energy_zero_parameter=FermiLevel,
                                    spectrum_method=spectrum_method)
    
    filename = f"{lc_string}_{strain_string}/{cation}N_PW_PDOS_{lc_string}_{strain_string}"
    nlsave(filename + ".hdf5", pdos, object_id="M")



def simulate_and_save_objects_FCC(optimal_alpha, bulk, cation, lattice_constant, strain):

    lc_string = "lc_{:.3f}".format(lattice_constant)
    if strain < 0:
        strain_string = "s_m{:d}".format(round(100*abs(strain)))
    else:
        strain_string = "s_{:d}".format(round(100*strain))

    os.makedirs(f"{lc_string}_{strain_string}", exist_ok = True)

    configuration = run_simulation(optimal_alpha, bulk, lc_string, cation, strain_string)

    bandstructure = compute_and_save_bandstructure_FCC(configuration, lc_string, cation, strain_string)

    compute_and_save_pdos_FCC(configuration, bandstructure, lc_string, cation, strain_string)


def simulate_and_save_objects_WZ(optimal_alpha, bulk, cation, a, c, strain):

    lc_string = "lc_{:.3f}_{:.3f}".format(a, c)
    if strain < 0:
        strain_string = "s_m{:d}".format(round(100*abs(strain)))
    else:
        strain_string = "s_{:d}".format(round(100*strain))

    os.makedirs(f"{lc_string}_{strain_string}", exist_ok = True)

    configuration = run_simulation(optimal_alpha, bulk, lc_string, cation, strain_string)

    bandstructure = compute_and_save_bandstructure_WZ(configuration, lc_string, cation, strain_string)

    compute_and_save_pdos_WZ(configuration, bandstructure, lc_string, cation, strain_string)
    

 
def main_FCC(cation, alpha, strained_constants, strains):
    
    for const, strain in zip(strained_constants, strains):
        bulk = setup_bulk_FCC(cation, const)
        simulate_and_save_objects_FCC(alpha, bulk, cation, const, strain)


def main_WZ(cation, alpha, strained_a, strained_c, strains):
    
    for a, c, strain in zip(strained_a, strained_c, strains):
        bulk = setup_bulk_WZ(cation, a, c)
        simulate_and_save_objects_WZ(alpha, bulk, cation, a, c, strain)
    

if __name__ == '__main__':
    
    #args = sys.argv[1:]
    FCC_OR_WZ = 1  # FCC=0, WZ=1
    args = ["Ga", 0.39, 4.531] #FOR FCC
    #args = ["Ga", 0.399, 3.189, 5.185]  #FOR WZ
    
    if FCC_OR_WZ == 0:
        if len(args) < 3:
            raise ValueError("Please insert the follwing quantities:\n1) Cation element (Ga, Al, In);\n2) Alpha mixing parameter;"
                            "\n3) Lattice constant (Ang).")
        cation, alpha, lattice_constant = args

        try:
            lattice_constant = float(lattice_constant)
        except:
            raise ValueError("Please provide a float value for the lattice constant (Ang).")
        
    if FCC_OR_WZ == 1:
        if len(args) < 4:
            raise ValueError("Please insert the follwing quantities:\n1) Cation element (Ga, Al, In);\n2) Alpha mixing parameter;"
                            "\n3) Lattice constant a (Ang)\n4) Lattice constant c (Ang).")
        cation, alpha, a, c = args

        try:
            a = float(a)
            c = float(c)
        except:
            raise ValueError("Please provide a float value for the lattice constants a and c (Ang).")

    cation = str(cation)
    try:
        alpha = float(alpha)
    except:
        raise ValueError("Please provide a float value for the alpha parameter.")
    
    strain_min = -0.06
    strain_max =  0.06
    step = 0.01
    strains = np.arange(strain_min, strain_max + step, step)
    
    # Remove strain = 0
    # strains = strains[np.where(np.abs(strains)>1e-15)]

    if FCC_OR_WZ == 0:
        lcs = lattice_constant + strains * lattice_constant
        
        if isMainProcess():
            print(f"Hydrostatic strain: {cation}N Zb cell, lattice constant of {lattice_constant} Ang,"
                f"strained from {strain_min} to {strain_max}\n\n")
        
        main_FCC(cation, alpha, lcs, strains)

    if FCC_OR_WZ == 1:
        a_strained = a + strains * a
        c_strained = c + strains * c

        if isMainProcess():
            print(f"Hydrostatic strain: {cation}N Wz cell, lattice constant a of {a} Ang, lattice constant c of {c} Ang"
                f"strained from {strain_min} to {strain_max}\n\n")
        
        main_WZ(cation, alpha, a_strained, c_strained, strains)

    

        
