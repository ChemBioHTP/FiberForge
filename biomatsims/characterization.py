




def calculate_cross_sectional_area(fibril_axis, pdb_file, probe_size=0.6, n_cross_sections=10):

    def project_onto_surface(vector_to_project, surface_normal):
        surface_normal = surface_normal / np.linalg.norm(surface_normal)
        scalar_projection = np.dot(vector_to_project, surface_normal)
        vector_projection = scalar_projection * surface_normal
        projected_vector = vector_to_project - vector_projection
        return projected_vector
    
    # Load the GRO file using MDTraj
    traj = mdtraj.load(pdb_file)

    # Select the protein atoms
    protein_traj = traj.atom_slice(traj.top.select('protein'))

    # Calculate the van der Waals surface using MDTraj for the protein
    vdw_surface = mdtraj.shrake_rupley(protein_traj, probe_radius=probe_size)

    # Identify the indices of atoms that are part of the fibril
    fibril_atoms_indices = np.array([atom.index for atom in traj.topology.atoms if atom.residue.is_protein])

    # Calculate the cross-sectional area for each frame along the fibril axis
    cross_section_areas = []
    for frame_idx in range(traj.n_frames):
        # fibril_frame = traj.xyz[frame_idx][fibril_atoms_indices]
        # axis_frame = traj.xyz[frame_idx][fibril_axis]

        # Project the van der Waals surface onto the fibril axis
        projected_surface = np.array([project_onto_surface(v, fibril_axis) for v in vdw_surface[frame_idx]])

        # Calculate cross-sectional area by summing the projected surface
        cross_section_area = np.sum(projected_surface)

        cross_section_areas.append(cross_section_area)

    # Calculate the average cross-sectional area
    average_cross_section_area = np.mean(cross_section_areas)

    return average_cross_section_area


def calculate_var_over_time(xvg_file):
    with open(xvg_file, 'r') as f:
        lines = f.readlines()
        time_data = []
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                time_data.append(list(map(float, line.split())))
    return time_data