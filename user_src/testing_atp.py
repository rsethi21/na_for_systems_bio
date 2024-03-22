import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pdb
import json
import os
from itertools import product as P
from tqdm import tqdm
import sys
sys.path.append(f"{sys.path[0]}/..")

from src.substrate import Substrate
from src.interaction import Interaction
from src.rate import Rate
from src.network import Network
from src.parse import parse_rates, parse_substrates, parse_interactions
from src.score import error

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--substrates", help="substrates csv", required=True)
parser.add_argument("-r", "--rates", help="rates csv", required=True)
parser.add_argument("-i", "--interactions", help="interactions csv", required=True)
parser.add_argument("-o", "--output", help="path to the output directory to save output files", required=False, default=".")
parser.add_argument("-p", "--parameters", help="pretrained parameters json file", required=False, default=None)
parser.add_argument("-n", "--number", help="number of random samples for area plot", type=int, required=False, default=500)
parser.add_argument("-t", "--toplot", help="list of substrates to plot", nargs="*", required=False, default=[])

if __name__ == "__main__":
    args = parser.parse_args()
   
    # parse and create all necessary objects for creating a network
    rates = parse_rates(args.rates)
    substrates = parse_substrates(args.substrates, rates)
    interactions = parse_interactions(args.interactions, rates, substrates)
    
    # create a new instance of a network
    network = Network("example", rates, interactions, substrates)
    
    # visualize network ordinary differential equations
    time = np.linspace(0, 500, num=501)
    derivatives = network.get_representation_dydt()
    for l, t in derivatives.items():
        print(f"{l} = {t}")
    print()

    # load in pretrained parameters if any
    if args.parameters != None:
        with open(args.parameters, "r") as fitted_params:
            parameters = json.load(fitted_params)
        network.set_parameters(list(parameters.values()), list(parameters.keys()))
 
    amts = [1.0]
    ranges = [[60,120]]
    combos = list(P(amts, ranges))
    combo_figure = plt.figure()
    for i, combo in tqdm(enumerate(combos), total=len(amts)*len(ranges)):
        network.substrates["ATP"].max_value = combo[0]
        network.substrates["ATP"].time_ranges = [combo[1]]

    # check post-training training error
        mean_y_regular = network.graph_distributions(time, args.number, substrates_to_plot=["ATP"], normalize=False)
        print(f"Originally r4 was: {network.parameters['r4'].value}")
        network.parameters["r4"].value = 0.15
        mean_y_decrP2Y12 = network.graph_distributions(time, args.number, substrates_to_plot=["ATP"], normalize=False)
        print(f"Adjusted r4 was: {network.parameters['r4'].value}")
        index = list(network.substrates.keys()).index("pAKT")
        plt.plot(range(len(mean_y_regular[:,index])), mean_y_regular[:,index], label="ctrl")
        plt.plot(range(len(mean_y_decrP2Y12[:,index])), mean_y_decrP2Y12[:,index], label="cond")
    plt.legend()
    combo_figure.savefig("./combo_fig.png")
        # temp_df = pd.DataFrame(mean_y, columns=list(network.substrates.keys()))
        # temp_df.to_csv(os.path.join(args.output, f"atp_{combo[0]}_{combo[1]}"))
        
