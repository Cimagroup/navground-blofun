from perdiver.distances import *

import matplotlib as mpl

import iblofunmatch.inter as ibfm
output_dir="output"


def plot_matching_diagram_trajectories(Dist_X, Dist_Y, Dist_Z, ax, color="blue", print_barcode_n_reps=False):
    # Compute X, Z barcodes and matching
    # half_shift=int(shift_time/2)
    # Dist_X = trajectory_distance_weighted_velocities_2Dtorus(positions[:half_shift], velocities[:half_shift], weight, side)
    # Dist_Y = trajectory_distance_weighted_velocities_2Dtorus(positions, velocities, weight, side)
    # Dist_Z = np.minimum(Dist_X, Dist_Y)
    idx_S = list(range(Dist_X.shape[0]))
    # Compute matching from X to Z
    ibfm_out = ibfm.get_IBloFunMatch_subset(Dist_X, Dist_Z, idx_S, output_dir, max_rad=-1, num_it=1, store_0_pm=True, points=False, max_dim=1)
    # Plot 0 persistence diagram of matching 
    match_diagram = []
    for idx, bar_X in enumerate(ibfm_out["S_barcode_0"]):
        idx_match = ibfm_out["induced_matching_0"][idx]
        bar_Z = ibfm_out["X_barcode_0"][idx_match]
        match_diagram.append([bar_X[1], bar_Z[1]])
    # end for
    match_diagram = np.array(match_diagram)
    # Plot matching diagram
    ax.scatter(match_diagram[:,0], match_diagram[:,1], color=color)
    ax.plot([0,np.max(match_diagram)*1.1], [0,np.max(match_diagram)*1.1], color="gray")
    if print_barcode_n_reps:
        match_diagram = []
        for idx, bar_X in enumerate(ibfm_out["S_barcode_0"]):
            idx_match = ibfm_out["induced_matching_0"][idx]
            bar_Z = ibfm_out["X_barcode_0"][idx_match]
            match_diagram.append([bar_X[1], bar_Z[1]])
        # end for
        match_diagram = np.array(match_diagram)
        
        print(np.array(match_diagram))
        print(ibfm_out["S_reps_0"])
    # end printing barcode and representatives
# end persistence diagram function


def plot_sequence(X_list, ax, mark_points=[], color="red"):
    # Plot figure
    X_old = X_list[0]
    ax.scatter(X_old[:,0], X_old[:,1], s=20, marker="o", color=mpl.colormaps["RdBu"](1/(len(X_list)+1)), zorder=2)
    for idx, X in enumerate(X_list[1:]):
        ax.scatter(X[:,0], X[:,1], s=20, marker="o", color=mpl.colormaps["RdBu"]((idx+1)/(len(X_list)+1)), zorder=2, label="X")
        if len(mark_points)>0:
            mark_X = X[mark_points]
            ax.scatter(mark_X[:,0], mark_X[:,1], s=20, marker="+", color=color, zorder=3)
        X_old = X
    #end for 

def compute_divergence_vector(Dist_X, Dist_Y, Dist_Z):
    assert(Dist_X.shape[0]==Dist_Y.shape[0])
    assert(Dist_X.shape[0]==Dist_Z.shape[0])
    idx_S = list(range(int(Dist_X.shape[0])))
    # Compute matching
    ibfm_out = [
        ibfm.get_IBloFunMatch_subset(Dist_X, Dist_Z, idx_S, output_dir, max_rad=-1, num_it=1, store_0_pm=True, points=False, max_dim=1),
        ibfm.get_IBloFunMatch_subset(Dist_Y, Dist_Z, idx_S, output_dir, max_rad=-1, num_it=1, store_0_pm=True, points=False, max_dim=1)
    ]
    # Compute induced matchings
    matching_XZ = ibfm_out[0]["induced_matching_0"]
    matching_YZ = ibfm_out[1]["induced_matching_0"]
    composition_XY = [matching_YZ.index(i) for i in matching_XZ]
    endpoints_0 = np.array(ibfm_out[0]["S_barcode_0"][:,1])
    endpoints_1 = np.array(ibfm_out[1]["S_barcode_0"][:,1])
    endpoints_1 = endpoints_1[composition_XY]
    return  endpoints_1-endpoints_0, ibfm_out[0]["X_barcode_0"][:,1]

def compute_divergence_Z_arrays(positions, velocities, steps_list, shift_step, weight, side):
    divergence_list = []
    Z_barcodes_list = []
    for start_step in steps_list:
        Dist_X, Dist_Y, Dist_Z = compute_distance_matrices_trajectories_2D_torus(positions, velocities, start_step, shift_step, weight, side)
        divergence, Z_barcode = compute_divergence_vector(Dist_X, Dist_Y, Dist_Z)
        divergence_list.append(divergence)
        Z_barcodes_list.append(Z_barcode)
    # want x-coordinate: steps, and y-coordinate: interval idx
    divergence_arr = np.array(divergence_list).transpose()
    Z_barcodes_arr = np.array(Z_barcodes_list).transpose()
    return divergence_arr, Z_barcodes_arr

def compute_cumulative_array(div_arr):
    cumulative_list = []
    for j, divergence in enumerate(div_arr.transpose()):
        cumulative = []
        for i in range(div_arr.shape[0]):
            cumulative.append(np.sum(divergence[:i+1]))
        cumulative_list.append(cumulative)
    return np.array(cumulative_list).transpose()