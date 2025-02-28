### This file contains a few convenient functions and objects to interact with Navground

import argparse
from navground import sim, core

import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Simulation Parameters')
parser.add_argument('--scenario', type=str, 
                    choices=['Corridor', 'Cross', 'CrossTorus'],
                    help='Choose a scenario from the list: Corridor, Cross, CrossTorus')
parser.add_argument('--side', type=float, default=10.0, help='Side of the environment square')
parser.add_argument('--length', type=float, default=10.0, help='Length of the environment')
parser.add_argument('--width', type=float, default=10.0, help='Width of the environment')
parser.add_argument('--num_runs', type=int, default=1, help='Number of simulation runs')
parser.add_argument('--num_steps', type=int, default=100, help='Number of steps in the simulation')
parser.add_argument('--time_step', type=float, default=0.1, help='Time step for the simulation')
parser.add_argument('--num_agents', type=int, default=10, help='Number of agents in the simulation')
parser.add_argument('--max_speed', type=float, default=1.0, help='Maximum speed of agents')
parser.add_argument('--optimal_speed_min', type=float, default=1.0, help='Minimal optimal speed')
parser.add_argument('--optimal_speed_max', type=float, default=1.0, help='Maximum optimal speed')
parser.add_argument('--radius', type=float, default=0.25, help='Radius of agents')
parser.add_argument('--safety_margin', type=float, default=0.1, help='Safety margin for agents')
parser.add_argument('--behavior', type=str, default='HL', help='Behavior type')
parser.add_argument('--max_edge_length', type=float, default=100.0, help='Maximum edge length in the simplicial complex')
parser.add_argument('--time_delay', type=int, default=1, help='Time delay to analise simulation intervals')
parser.add_argument('--embedding_length', type=int, default=10, help='Length of the simulation intervals')
parser.add_argument('--epsilon', type=int, default=50, help='time differences for matching and bottleneck distance computation')


def run_navground(args):
    yaml = f"""
    runs: {args.num_runs}
    steps: {args.num_steps}
    time_step: {args.time_step}
    save_directory: ''
    record_pose: true
    record_twist: true
    record_collisions: true
    record_deadlocks: true
    record_safety_violation: true
    record_efficacy: true
    scenario:
      type: {args.scenario}""" 
    # Give special parameters for each scenario type
    if args.scenario == "Corridor":
        yaml+= f"""
      width: {args.width}
      length: {args.length}"""
    else:
        yaml+= f"""
      side: {args.side}"""
    # Continue with other parameters 
    yaml+= f"""
      tolerance: 1
      groups:
        -
          type: thymio
          number: {args.num_agents}
          radius: {args.radius}
          control_period: 0.1
          speed_tolerance: 0.02
          kinematics:
            type: 2WDiff
            wheel_axis: 0.094
            max_speed: {args.max_speed}
          behavior:
            type: {args.behavior}
            optimal_speed:
                sampler: uniform
                from: {args.optimal_speed_min}
                to: {args.optimal_speed_max}
            horizon: 5.0
            safety_margin: {args.safety_margin}
          state_estimation:
            type: Bounded
            range: 5.0
    """
    experiment = sim.load_experiment(yaml)
    experiment.run()
    return experiment.runs


### The following functions are for extracting collisions, efficacy and deadlocks from Navground
### These were taken from the Navground tutorial notebook from the following URL: 
### https://idsia-robotics.github.io/navground/tutorials/deadlocks_and_collisions.html

def count_deadlocks(deadlock_time, initial_time, final_time):
    is_deadlocked = np.logical_and(deadlock_time > initial_time, deadlock_time < (final_time - 5.0))
    return sum(is_deadlocked)

def extract_data(runs, initial_step, final_step):
    collisions = []
    deadlocks = []
    efficacy = []
    sms = []
    seeds = []
    for i in runs.keys():
        run = runs[i]
        world = run.world
        initial_time, final_time = initial_step*run.time_step, final_step*run.time_step
        deadlocks.append(count_deadlocks(np.array(run.deadlocks), initial_time, final_time))
        collisions.append(np.sum(np.logical_and(
            initial_step < run.collisions[:,0], run.collisions[:,0] < final_step
        )))
        efficacy.append(np.array(run.efficacy[initial_step:final_step]).mean())

    df = pd.DataFrame({
        'deadlocks': deadlocks,
        'collisions': collisions,
        'efficacy': efficacy})
    df['safe'] = (df.collisions == 0).astype(int)
    df['fluid'] = (df.deadlocks == 0).astype(int)
    df['ok'] = ((df.deadlocks == 0) & (df.collisions == 0)).astype(int)
    return df
    