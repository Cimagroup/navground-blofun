import scipy.spatial.distance as dist
import numpy as np

def distances_weighted_velocities(points, velocities, weight):
    points_velocities = np.hstack((points, velocities*weight))
    return dist.squareform(dist.pdist(points_velocities, "minkowski", p=2))

def trajectory_distance_weighted_velocities(positions, velocities, weight):
    assert(len(positions)>0)
    assert(len(positions)==len(velocities))
    positions_velocities_list = []
    for idx, points in enumerate(positions):
        positions_velocities_list.append(np.hstack((points, velocities[idx]*weight)))
    distances_list = []
    for idx, points_vel in enumerate(positions_velocities_list):
        # Compare trajectories at same time
        points_vel_compare = positions_velocities_list[idx]
        distances_list.append(dist.cdist(points_vel, points_vel_compare, "minkowski", p=2))
    # end for
    distances_arr = np.array(distances_list)
    return np.min(distances_arr, axis=0)

def trajectory_distance_weighted_velocities_2Dtorus(positions, velocities, weight, side):
    assert(len(positions)>0)
    assert(len(positions)==len(velocities))
    positions_velocities_list = []
    # Horizontal shift [side, 0] 
    for idx, points in enumerate(positions):
        points_shift = np.array(points)
        points_shift[points_shift[:,0] < (side/2)] += [side, 0]
        positions_velocities_list.append(np.hstack((points_shift, velocities[idx]*weight)))
    # Vertical shift [0, side]
    for idx, points in enumerate(positions):
        points_shift = np.array(points)
        points_shift[points_shift[:,1] < (side/2)] += [0, side]
        positions_velocities_list.append(np.hstack((points_shift, velocities[idx]*weight)))
    # Diagonal shift [side, side]
    for idx, points in enumerate(positions):
        points_shift = np.array(points)
        # Take into account torus periodicity
        points_shift[np.logical_and(points_shift[:,0] <  (side/2), points_shift[:,1] < (side/2))]  += [side, side]
        points_shift[np.logical_and(points_shift[:,0] <  (side/2), points_shift[:,1] >= (side/2))] += [side, 0]
        points_shift[np.logical_and(points_shift[:,0] >= (side/2), points_shift[:,1] < (side/2))]  += [0, side]
        positions_velocities_list.append(np.hstack((points_shift, velocities[idx]*weight)))
    distances_list = []
    for idx, points_vel in enumerate(positions_velocities_list):
        # Compare trajectories at same time
        points_vel_compare = positions_velocities_list[idx]
        distances_list.append(dist.cdist(points_vel, points_vel_compare, "minkowski", p=2))
    # end for
    distances_arr = np.array(distances_list)
    return np.min(distances_arr, axis=0)

def compute_distance_matrices_trajectories_2D_torus(positions, velocities, start_step, shift_time, weight, side):
    half_shift=int(shift_time/2)
    start_X, start_Y = start_step, start_step + half_shift
    end_X, end_Y = start_step+half_shift, start_step + shift_time
    Dist_X = trajectory_distance_weighted_velocities_2Dtorus(positions[start_X:end_X], velocities[start_X:end_X], weight, side)
    Dist_Y = trajectory_distance_weighted_velocities_2Dtorus(positions[start_Y:end_Y], velocities[start_Y:end_Y], weight, side)
    Dist_Z = np.minimum(Dist_X, Dist_Y)
    return Dist_X, Dist_Y, Dist_Z