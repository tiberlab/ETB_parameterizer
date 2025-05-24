import numpy as np
import scipy
import sys
import os
import fnmatch


def create_header(filename, kpoint, lc, strain):
    with open(filename, "a+") as f:
        f.write(f"# k-point: {kpoint}\n")
        f.write(f"# lattice constant: {lc}, strain: {strain}\n")
        f.write(f"# units energies: [eV]; units pdos: [1/eV]\n")
        f.write(f"# energies sc pc dc sa pa da s*\n\n")

    return



def create_coefficient_report(out_path, filename, kpoint, coeffs_for_report):
        
    with open(f"{out_path}/Coefficients_normalization_report.dat", "a+") as f:
        f.write(f"Sum of coefficients in {filename}_{kpoint}.csv:")
        sum_vector = np.sum(np.array(coeffs_for_report), 0)
        f.write(f"\n{sum_vector.shape[0]} energies\n")
        string_row = ""
        for s in sum_vector:
            string_row += f"{s}, "
        f.write(string_row)
        f.write("\n\n")
    return



def gaussian(energy1, energy2, sigma):

    return 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-0.5 * ((energy1 - energy2)/(sigma))**2 )



def main(in_path, filename, out_path, k_list, sigma, FCC_OR_WZ):

    for kpoint in k_list:
        # create_header(f"{filename}_{kpoint}.csv", kpoint, lc, strain)

        pdos = nlread(f"{in_path}/{filename}.hdf5", object_id=kpoint)[0]

        pdos_projections = pdos.evaluate().inUnitsOf(eV**-1)
        energies = pdos.energies().inUnitsOf(eV)

        if FCC_OR_WZ == 0:
            orbital_names = ["sc", "pc-1", "pc0", "pc1" ,"dc-2", "dc-1", "dc0", "dc1", "dc2", \
                             "sa", "pa-1", "pa0", "pa1" ,"da-2", "da-1", "da0", "da1", "da2"]
        elif FCC_OR_WZ == 1:
            orbital_names = ["s1c", "p1c-1", "p1c0", "p1c1" ,"d1c-2", "d1c-1", "d1c0", "d1c1", "d1c2", \
                             "s2c", "p2c-1", "p2c0", "p2c1" ,"d2c-2", "d2c-1", "d2c0", "d2c1", "d2c2", \
                             "s1a", "p1a-1", "p1a0", "p1a1" ,"d1a-2", "d1a-1", "d1a0", "d1a1", "d1a2"\
                             "s2a", "p2a-1", "p2a0", "p2a1" ,"d2a-2", "d2a-1", "d2a0", "d2a1", "d2a2"]
 
        coeffs_for_report = []

        with open(f"{out_path}/{filename}_{kpoint}.csv", "w+") as f:
            
            print(f"Solving for {filename}/{kpoint} ")

            for orbital in range(len(orbital_names)):
                linsys_A = []
                linsys_b = []
                orbital_pdos = pdos_projections[orbital,:]

                for i, energy1 in enumerate(energies):
                    b_term = orbital_pdos[i]
                    linsys_b.append(b_term)
                    
                    row = []
                    for energy2 in energies:
                        g = gaussian(energy1, energy2, sigma)
                        row.append(g)
                    row = np.array(row)
                    linsys_A.append(row)

                linsys_A = np.array(linsys_A)
                linsys_b = np.array(linsys_b)

                orbital_coefficients = scipy.linalg.solve(linsys_A, linsys_b)

                coeffs_for_report.append(orbital_coefficients)


                f.write(f"{orbital_names[orbital]}\n")
                string_row = ""
                for coeff in orbital_coefficients:
                    if coeff < 0:
                        string_row += f"{0.0}, "
                    else:   
                        string_row += f"{coeff}, "

                f.write(string_row + "\n\n")

        create_coefficient_report(out_path, filename, kpoint, coeffs_for_report)
    
    return


if __name__ == "__main__":

    try:
        in_path = sys.argv[1]
        out_path = sys.argv[2]

    except:
        raise ValueError("Please provide the input path (location of biax/shear/uniax folders) and the output path.")

    FCC_OR_WZ = 0  # 0 for FCC, 1 for WZ
    sigma = 0.01  # Has to match the Gaussian broadening parameter in the DFT PDOS calculation

    if FCC_OR_WZ == 0:
        k_list = ["L", "G", "X"]
    elif FCC_OR_WZ == 1:
        k_list = ["Y", "G", "Z"]

    
    if os.path.exists(f"{out_path}/Coefficients_normalization_report.dat"):
        os.remove(f"{out_path}/Coefficients_normalization_report.dat")

    folders = [name for name in os.listdir(in_path) if ("biax" in name and os.path.isdir(f"{in_path}/{name}"))]
    for fold in folders:
        pdos_files = fnmatch.filter(os.listdir(f"{in_path}/{fold}"), "*PDOS*.hdf5")

        for file in pdos_files:
            filename = file.split(".hdf5")[0]
            main(f"{in_path}/{fold}", filename, out_path, k_list, sigma, FCC_OR_WZ)

        
    folders = [name for name in os.listdir(in_path) if ("uniax" in name and os.path.isdir(f"{in_path}/{name}"))]
    for fold in folders:
        pdos_files = fnmatch.filter(os.listdir(f"{in_path}/{fold}"), "*PDOS*.hdf5")

        for file in pdos_files:
            filename = file.split(".hdf5")[0]
            main(f"{in_path}/{fold}", filename, out_path, k_list, sigma, FCC_OR_WZ)

        
    folders = [name for name in os.listdir(in_path) if ("shear" in name and os.path.isdir(f"{in_path}/{name}"))]
    for fold in folders:
        pdos_files = fnmatch.filter(os.listdir(f"{in_path}/{fold}"), "*PDOS*.hdf5")

        for file in pdos_files:
            filename = file.split(".hdf5")[0]
            main(f"{in_path}/{fold}", filename, out_path, k_list, sigma, FCC_OR_WZ)


