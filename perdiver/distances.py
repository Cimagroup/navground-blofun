import scipy.spatial.distance as dist
import numpy as np

def check_points_velocities(points, velocities):
    if (not np.all(points.shape==velocities.shape)) or (points.shape[1]!=2):
        raise(ValueError)
        
def distances_weighted_velocities(points, velocities, weight):
    """ This function takes positions and velocities of a moving point cloud at a given time and returns 
    the distance matrix where velocities are weighted by the parameter `weight`. The metric space 
    corresponds to points given as four coordinates (x,y,wv_x, wv_y) where `w`is the weight, `x`, `y` are the coordinates
    and `v_x`, `v_y` are the vector components.
    """
    check_points_velocities(points, velocities)
    points_velocities = np.hstack((points, velocities*weight))
    return dist.squareform(dist.pdist(points_velocities, "minkowski", p=2))

def distances_corridor_weighted_velocities(points, velocities, weight, length):
    """ Returns the distance matrix of points in a corridor. Velocities are weighted by the parameter `weight`.
    The metric space corresponds to points given as four coordinates (x, y, wv_x, wv_y) where `w`is the weight, 
    `x`, `y` are the coordinates and `v_x`, `v_y` are the vector components.
    This requieres a specific function because of the periodic boundary conditions.
    `length` is the corridor length.
    """
    check_points_velocities(points, velocities)
    points_velocities = np.hstack((points, velocities*weight))
    dist_0 = dist.squareform(dist.pdist(points_velocities, "minkowski", p=2))
    shift_points_vels = np.array(points_velocities) # make a copy
    left_pts_idx = shift_points_vels[:,0] < length/2
    shift_points_vels[left_pts_idx] += [length,0,0,0]
    dist_1 = dist.squareform(dist.pdist(shift_points_vels, "minkowski", p=2))
    return np.minimum(dist_0, dist_1)

def distances_2Dtorus_weighted_velocities(points, velocities, weight, side):
    """ Returns the distance matrix of points in a square with side `side` and with periodic
    boundary conditions, so that it is a torus. Velocities are weighted by the parameter `weight`.
    The metric space corresponds to points given as four coordinates (x, y, wv_x, wv_y) where `w`is the weight, 
    `x`, `y` are the coordinates and `v_x`, `v_y` are the vector components.
    This requieres a specific function because of the periodic boundary conditions.
    """
    check_points_velocities(points, velocities)
    ####
    positions_velocities_list = []
    # Horizontal shift [side, 0] 
    points_shift = np.array(points)
    points_shift[points_shift[:,0] < (side/2)] += [side, 0]
    positions_velocities_list.append(np.hstack((points_shift, velocities*weight)))
    # Vertical shift [0, side]
    points_shift = np.array(points)
    points_shift[points_shift[:,1] < (side/2)] += [0, side]
    positions_velocities_list.append(np.hstack((points_shift, velocities*weight)))
    # Diagonal shift [side, side]
    points_shift = np.array(points)
    # Take into account torus periodicity
    points_shift[np.logical_and(points_shift[:,0] <  (side/2), points_shift[:,1] < (side/2))]  += [side, side]
    points_shift[np.logical_and(points_shift[:,0] <  (side/2), points_shift[:,1] >= (side/2))] += [side, 0]
    points_shift[np.logical_and(points_shift[:,0] >= (side/2), points_shift[:,1] < (side/2))]  += [0, side]
    positions_velocities_list.append(np.hstack((points_shift, velocities*weight)))
    ### Now, compute the distances for each instance of points with applied periodic conditions
    distances_list = []
    for idx, points_vel in enumerate(positions_velocities_list):
        # Compare trajectories at same time
        points_vel_compare = positions_velocities_list[idx]
        distances_list.append(dist.cdist(points_vel, points_vel_compare, "minkowski", p=2))
    # end for
    #### Take minimum distances and return
    distances_arr = np.array(distances_list)
    return np.min(distances_arr, axis=0)

### Functions that return distance matrices corresponding to a pair of timesteps

def compute_distance_matrices_timesteps_corridor(X, Y, vel_X, vel_Y, weight, length):
    Dist_X = distances_corridor_weighted_velocities(X, vel_X, weight, length)
    Dist_Y = distances_corridor_weighted_velocities(Y, vel_Y, weight, length)
    Dist_Z = np.minimum(Dist_X, Dist_Y)
    return Dist_X, Dist_Y, Dist_Z

def compute_distance_matrices_timesteps_2Dtorus(X, Y, vel_X, vel_Y, weight, side):
    Dist_X = distances_2Dtorus_weighted_velocities(X, vel_X, weight, side)
    Dist_Y = distances_2Dtorus_weighted_velocities(Y, vel_Y, weight, side)
    Dist_Z = np.minimum(Dist_X, Dist_Y)
    return Dist_X, Dist_Y, Dist_Z

### Functions that return distance matrix corresponding to trajectories

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

def trajectory_corridor_distance_weighted_velocities(positions, velocities, weight, length):
    assert(len(positions)>0)
    assert(len(positions)==len(velocities))
    positions_velocities_list = []
    for idx, points in enumerate(positions):
        positions_velocities_list.append(np.hstack((points, velocities[idx]*weight)))
    distances_list = []
    for idx, points_vel in enumerate(positions_velocities_list):
        # Compare trajectories at same time
        points_vel_compare = positions_velocities_list[idx]
        dist_0 = dist.cdist(points_vel, points_vel_compare, "minkowski", p=2)
        shift_points_vel_compare = np.array(points_vel_compare)
        shift_points_vel_compare[shift_points_vel_compare[:,0]<length/2]+=[length, 0, 0, 0]
        shift_points_vel = np.array(points_vel)
        shift_points_vel[shift_points_vel[:,0]<length/2]+=[length, 0, 0, 0]
        dist_1 = dist.cdist(shift_points_vel, shift_points_vel_compare, "minkowski", p=2)
        distances_list.append(np.minimum(dist_0, dist_1))
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


def compute_distance_matrices_trajectories_corridor(positions, velocities, start_step, shift_time, weight, length):
    half_shift=int(shift_time/2)
    start_X, start_Y = start_step, start_step + half_shift
    end_X, end_Y = start_step+half_shift, start_step + shift_time
    Dist_X = trajectory_corridor_distance_weighted_velocities(positions[start_X:end_X], velocities[start_X:end_X], weight, length)
    Dist_Y = trajectory_corridor_distance_weighted_velocities(positions[start_Y:end_Y], velocities[start_Y:end_Y], weight, length)
    Dist_Z = np.minimum(Dist_X, Dist_Y)
    return Dist_X, Dist_Y, Dist_Z

def compute_distance_matrices_trajectories_2D_torus(positions, velocities, start_step, shift_time, weight, side):
    half_shift=int(shift_time/2)
    start_X, start_Y = start_step, start_step + half_shift
    end_X, end_Y = start_step+half_shift, start_step + shift_time
    Dist_X = trajectory_distance_weighted_velocities_2Dtorus(positions[start_X:end_X], velocities[start_X:end_X], weight, side)
    Dist_Y = trajectory_distance_weighted_velocities_2Dtorus(positions[start_Y:end_Y], velocities[start_Y:end_Y], weight, side)
    Dist_Z = np.minimum(Dist_X, Dist_Y)
    return Dist_X, Dist_Y, Dist_Z