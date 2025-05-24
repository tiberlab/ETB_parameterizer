import numpy as np
import sys
import os
import fnmatch
import warnings


def read_bands(filename):
    bands = nlread(filename + ".hdf5", Bandstructure)[0]
    bands_array = bands.evaluate().inUnitsOf(eV)
    bands_array -= bands.valenceBandEdge().inUnitsOf(eV)

    return bands, bands_array


def count_bands(b_array, k_idx_Gamma):

    pos_counts = 0
    for eigval in b_array[k_idx_Gamma, :]:
        if eigval > 0:
            pos_counts += 1

    return pos_counts, b_array.shape[1] - pos_counts


def find_gamma_index(bands):
    route = bands.route()
    if any(["G" in path for path in route]):
        path_with_gamma = 0
        for path in route:
            if path[0] == "G":
                break
            else:
                path_with_gamma += 1

        points_per_segment = (len(bands.kpoints()) -1) // len(route)
        gamma_index = int(points_per_segment * path_with_gamma)

        return gamma_index

    else:
        warnings.warn("No Gamma point in K-route. Info in file *info.dat file will be wrong.")
        return 0


def create_k_path(bands):

    frac_k_points = bands.kpoints()
    route = bands.route()
    lattice =  bands._lattice()

    num_segments = len(route)
    # Points per segment
    pps = (len(frac_k_points) -1) // num_segments

    k_points = [lattice.convertFractionalKPoint(kp).inUnitsOf(Ang**-1) for kp in frac_k_points]
    
    # Get the symmetry points
    symmetry_points = []
    for i in range(num_segments+1):
        symmetry_points.append(k_points[i*pps])

    # Calculate the length of each segment and total length
    total_length = 0
    lengths = []
    for i in range(1, len(symmetry_points)):
        l = np.linalg.norm(symmetry_points[i] - symmetry_points[i-1])
        total_length += l
        lengths.append(l)
    
    unnormalized_route = np.cumsum(lengths)
    unnormalized_route = np.insert(unnormalized_route, 0, 0.0)
    
    symmetry_pts_k_route = unnormalized_route / total_length

    # Now fill the k-route with equally spaced points between the symmetry points
    N_k = len(k_points)
    k_path = np.zeros(N_k)
    num_sym = len(symmetry_pts_k_route)    # How many symmetry points in the path

    for i in range(1,num_sym):
        segment_path = np.linspace(symmetry_pts_k_route[i-1], symmetry_pts_k_route[i], pps+1)
        if i != num_sym-1:
            k_path[(i-1)*pps:i*pps] = segment_path[:pps]
        else:
            k_path[(i-1)*pps:] = segment_path

    return k_path


def main(in_path, filename, output_path):
    bands, b_array = read_bands(f"{in_path}/{filename}")

    N_bands = b_array.shape[1]
    N_k = b_array.shape[0]
    
    # Create k-points with the QATK convention
    kpoints = create_k_path(bands)
    
    with open(f"{output_path}/{filename}.csv", "w+") as f:
        
        for i in range(N_bands):
            for j in range(N_k):
                f.write(f"{kpoints[j]}, {b_array[j, i]}\n")
            f.write("\n")

    Gamma_idx = find_gamma_index(bands)
    positive_bands, negative_bands = count_bands(b_array, Gamma_idx)

    with open(f"{output_path}/Other_strains_bands_info.dat", "a+") as f:
        f.write(f"File: {filename}\n")
        f.write(f"Band gap: {bands.directBandGap().inUnitsOf(eV)} eV\n")
        f.write(f"Valence band edge: {bands.valenceBandEdge().inUnitsOf(eV)} eV\n")
        f.write(f"Conduction band edge: {bands.conductionBandEdge().inUnitsOf(eV)} eV\n")
        f.write(f"Number of conduction bands: {positive_bands}\n")
        f.write(f"Number of valence bands: {negative_bands}\n")
        f.write(f"Path: {bands.route()}\n")
        f.write(f"Symmetry points and associated coordinates: {bands._symmetryPoints()}\n\n")


if __name__ == "__main__":

    try:
        in_path = sys.argv[1]
        output_path = sys.argv[2]

    except:
        raise ValueError("Please provide the input path (location of biax/shear/uniax folders) and the output path.")

    # Biaxial
    with open(f"{output_path}/Other_strains_bands_info.dat", "a+") as f:
        f.write("Biaxial strain:\n")

    folders = [name for name in os.listdir(in_path) if ("biax" in name and os.path.isdir(f"{in_path}/{name}"))]
    for fold in folders:
        filename = fnmatch.filter(os.listdir(f"{in_path}/{fold}"), '*bands*.hdf5')[0].split(".hdf5")[0]
        main(f"{in_path}/{fold}", filename, output_path)

    # Uniaxial
    with open(f"{output_path}/Other_strains_bands_info.dat", "a+") as f:
        f.write("\nUniaxial strain:\n")
        
    folders = [name for name in os.listdir(in_path) if ("uniax" in name and os.path.isdir(f"{in_path}/{name}"))]
    for fold in folders:
        filename = fnmatch.filter(os.listdir(f"{in_path}/{fold}"), '*bands*.hdf5')[0].split(".hdf5")[0]
        main(f"{in_path}/{fold}", filename, output_path)

    # Shear
    with open(f"{output_path}/Other_strains_bands_info.dat", "a+") as f:
        f.write("\nShear strain:\n")
        
    folders = [name for name in os.listdir(in_path) if ("shear" in name and os.path.isdir(f"{in_path}/{name}"))]
    for fold in folders:
        filename = fnmatch.filter(os.listdir(f"{in_path}/{fold}"), '*bands*.hdf5')[0].split(".hdf5")[0]
        main(f"{in_path}/{fold}", filename, output_path)
