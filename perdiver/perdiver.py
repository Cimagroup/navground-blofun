from perdiver.distances import *

import gudhi

import matplotlib as mpl

import tdqual.topological_data_quality_0 as tdqual

from navground.sim.ui.render import image_for_world

################################################################################
### Computing Matching diagram and divergence vectors
################################################################################

def get_matching_diagram(Dist_X, Dist_Y):
    """ Returns the matching diagram given two square matrices which are positive hollow symmetric."""
    # Compute matching from X to Y
    filt_X, filt_Y, matching = tdqual.compute_Mf_0(Dist_X, Dist_Y, is_dist=True)
    # Plot 0 persistence diagram of matching 
    match_diagram = []
    for idx, idx_match in enumerate(matching):
        match_diagram.append([filt_X[idx], filt_Y[idx_match]])
    # end for
    return np.array(match_diagram)
# end def

def compute_divergence_vector(Dist_X, Dist_Y):
    """ Returns the divergence vector given two square matrices which are positive hollow symmetric."""
    assert(Dist_X.shape[0]==Dist_Y.shape[0])
    # Compute matching
    matching_diagram = get_matching_diagram(Dist_X, Dist_Y)
    return  matching_diagram[:,1] - matching_diagram[:,0]
# end def

def compute_divergence_vector(matching_diagram):
    """ Given a matching diagram, returns the divergence vector."""
    return  matching_diagram[:,1] - matching_diagram[:,0]
# end def

################################################################################
### Other persistence divergence computations
################################################################################

def absolute_perdiver_signal(matching_diagram_list):
    absolute_perdiver = []
    for matching_diagram in matching_diagram_list:
        absolute_perdiver.append(np.sum(np.abs(matching_diagram[:,0]- matching_diagram[:,1])))
    return absolute_perdiver
# end def

def compute_cumulative_array(div_arr):
    cumulative_list = []
    for j, divergence in enumerate(div_arr.transpose()):
        cumulative = []
        for i in range(div_arr.shape[0]):
            cumulative.append(np.sum(divergence[:i+1]))
        cumulative_list.append(cumulative)
    return np.array(cumulative_list).transpose()
# end def

def bottleneck_matching_diagrams(diag_1, diag_2):
    """ Returns the bottleneck distance between two matching diagrams with the same size.
    It first computes that both diagrams have the same number of points. Then it computes the 
    bottleneck distance by first performing a translation of points and then using the bottleneck 
    distance from the GUDHI package. 
    """
    # Check that both matching diagrams have the same size.
    if diag_1.shape[0]!= diag_2.shape[0]:
        print("We assume that both diagrams have the same number of points!")
        raise ValueError
    # Translate points by a vertical shift (depending on the matching diagrams)
    max_val = max(np.max(diag_1), np.max(diag_2))
    diag_1_shift = diag_1 + [0, 2*max_val]
    diag_2_shift = diag_2 + [0, 2*max_val]
    # Use GUDHI to compute bottleneck distance
    return gudhi.bottleneck_distance(diag_1_shift, diag_2_shift)
# end bottleneck distance between diagrams

################################################################################
#### PLOTTING FUNCTIONS
################################################################################

def plot_matching_diagram(match_diagram, ax, color="blue", print_barcode_n_reps=False, max_val_diag=-1):
    """ Plots the matching diagram"""
    ax.scatter(match_diagram[:,0], match_diagram[:,1], color=color)
    if max_val_diag<0:
        max_val_diag = np.max(match_diagram)
    ax.plot([0,max_val_diag*1.1], [0,max_val_diag*1.1], color="gray")
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
# end plotting sequence 

def same_diagram_scale(ax_list):
    """ ax_list is a list of axes to scale to same dimensions"""
    # Scale both diagrams to same scale 
    xmax = np.max(np.vstack([ax.get_xlim() for ax in ax_list]))
    for ax in ax_list:
        plot_matching_diagram(np.array([[xmax, xmax]]), ax, color="white")

def plot_two_timesteps(X, Y, ax, X_col="blue", Y_col="red"):
    # Plot figure
    ax.scatter(X[:,0], X[:,1], s=20, marker="s", c=X_col, zorder=2, label="X")
    ax.scatter(Y[:,0], Y[:,1], s=23, marker="x", c=Y_col, zorder=2, label="Y")
    for edge in zip(X, Y):
        edge_pts = np.array(edge)
        ax.plot(edge_pts[:,0], edge_pts[:,1], c="gray", zorder=1)

def plot_two_timesteps_with_velocities(X, Y, vel_X, vel_Y, ax, X_col="blue", Y_col="red", arrow_width=0.05):
    # Plot two timesteps and connecting segment
    plot_two_timesteps(X, Y, ax, X_col=X_col, Y_col=Y_col)
    # Draw velocities 
    for pos, vel in zip(X, vel_X):
        ax.arrow(pos[0], pos[1], vel[0], vel[1], color=X_col, zorder=2, width=arrow_width)
    
    for pos, vel in zip(Y, vel_Y):
        ax.arrow(pos[0], pos[1], vel[0], vel[1], color=Y_col, zorder=2, width=arrow_width)

def decorate_red(agent):
    return {'fill': 'red'}

def decorate_blue(agent):
    return {'fill': 'blue'}

def image_run_timestep_list(run, timestep_list, to_decorate=False):
    poses = run.poses
    image_list = []
    # Collect all images for world on timesteps   
    for i, timestep in enumerate(timestep_list):
        run.go_to_step(timestep)
        world = run.world
        if to_decorate:
            if i==0:
                image_list.append(
                    image_for_world(world, background_color="black", relative_margin=0, width=1500, decorate=decorate_red)
                )
            else:
                image_list.append(
                    image_for_world(world, background_color="black", relative_margin=0, width=1500, decorate=decorate_blue)
                )
            # end if-else
        else:
            image_list.append(
                image_for_world(world, background_color="black", relative_margin=0, width=1500)
            )
    # Combine images 
    image_sum = np.max(np.array(image_list), axis=0)
    # Set up background to cream
    np.any(image_sum==0)
    zero_positions = np.where(image_sum==0)
    for x,y in zip(zero_positions[0], zero_positions[1]):
        image_sum[x,y]=[255,253,208]
    # end for
    return image_sum



#################################################################################################

#### Corridor depiction #########################################################################

def plot_image_corridor(run, image, timestep_list, length, width, ax):
    # Load poses and number of agents
    poses = run.poses
    num_agents = poses.shape[1]
    # Plot image and adjust bounds
    ax.imshow(image, extent=[0,length, 0, width])
    ax.set_xlim([0, length])
    ax.set_ylim([-0.05*width, 1.05*width])
    # Plot trajectories
    for i in range(num_agents):
        for f, t in zip(timestep_list[:-1], timestep_list[1:]):
            edge = np.array([poses[f,i], poses[t,i]])
            edge = edge[np.argsort(edge[:,0])]
            if np.abs(edge[1,0]-edge[0,0]) < length/2:
                ax.plot(edge[:,0], edge[:,1], color=mpl.colormaps['Set1'](i/(num_agents+3)))
            else:
                ax.plot(edge[:,0]+[length,0], edge[:,1], color=mpl.colormaps['Set1'](i/(num_agents+3)))
                ax.plot(edge[:,0]-[0,length], edge[:,1], color=mpl.colormaps['Set1'](i/(num_agents+3)))
            # end if-else
        # end for
    # end for
# end def

def plot_timesteps_corridor(run, timestep_list, length, width, ax):
    image = image_run_timestep_list(run, timestep_list)
    plot_image_corridor(run, image, timestep_list, length, width, ax)

#################################################################################################

#### Cross depiction     ########################################################################

def plot_image_cross(run, image, timestep_list, side, ax):
    # Load poses and number of agents
    poses = run.poses
    num_agents = poses.shape[1]
    # Plot image and adjust bounds
    ax.imshow(image, extent=[-side/2, side/2, -side/2, side/2])
    ax.set_xlim([-side/2, side/2])
    ax.set_ylim([-side/2, side/2])
    # Plot trajectories
    for i in range(num_agents):
        for f, t in zip(timestep_list[:-1], timestep_list[1:]):
            edge = np.array([poses[f,i], poses[t,i]])
            edge = edge[np.argsort(edge[:,0])]
            ax.plot(edge[:,0], edge[:,1], color=mpl.colormaps['Set1'](i/(num_agents+3)))
        # end for
    # end for
# end def

def plot_timesteps_cross(run, timestep_list, side, ax):
    image = image_run_timestep_list(run, timestep_list)
    plot_image_cross(run, image, timestep_list, side, ax)


#################################################################################################

#### Corridor depiction #########################################################################
def plot_image_cross_torus(run, image, timestep_list, side, ax):
    # Load poses and number of agents
    poses = run.poses
    num_agents = poses.shape[1]
    # Plot image and adjust bounds
    ax.imshow(image, extent=[0, side, 0, side])
    ax.set_xlim([0, side])
    ax.set_ylim([0, side])
    # Plot trajectories
    for i in range(num_agents):
        for f, t in zip(timestep_list[:-1], timestep_list[1:]):
            edge = np.array([poses[f,i], poses[t,i]])
            edge = edge[np.argsort(edge[:,0])]
            if np.abs(edge[1,0]-edge[0,0]) < side/2:
                if np.abs(edge[1,1]-edge[0,1]) < side/2:
                    ax.plot(edge[:,0], edge[:,1], color=mpl.colormaps['Set1'](i/(num_agents+3)))
                else: # y periodic
                    up_or_down = 1 if edge[1,1] < edge[0,1] else -1 # 1 is up, -1 is going downards
                    ax.plot(edge[:,0], edge[:,1] - [up_or_down * side, 0], color=mpl.colormaps['Set1'](i/(num_agents+3)))
                    ax.plot(edge[:,0], edge[:,1] + [0, up_or_down * side], color=mpl.colormaps['Set1'](i/(num_agents+3)))
                # end else
            else: # x periodic
                left_or_right = 1 if edge[1,0] < edge[0,0] else -1 # 1 is up, -1 is going downards
                ax.plot(edge[:,0] - [left_or_right * side, 0], edge[:,1], color=mpl.colormaps['Set1'](i/(num_agents+3)))
                ax.plot(edge[:,0] + [0, left_or_right * side], edge[:,1], color=mpl.colormaps['Set1'](i/(num_agents+3)))
            # end if-else
        # end for
    # end for
# end def

def plot_timesteps_cross_torus(run, timestep_list, side, ax, to_decorate=False):
    image = image_run_timestep_list(run, timestep_list, to_decorate)
    plot_image_cross_torus(run, image, timestep_list, side, ax)




