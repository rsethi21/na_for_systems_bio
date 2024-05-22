import numpy as np
import matplotlib.pyplot as plt
import argparse
import pdb
import json
import os
from tqdm import tqdm
import sys
sys.path.append(f"{sys.path[0]}/..")

from src.substrate import Substrate
from src.interaction import Interaction
from src.rate import Rate
from src.network import Network
from src.parse import parse_rates, parse_substrates, parse_interactions

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--substrates", help="substrates csv", required=True)
parser.add_argument("-r", "--rates", help="rates csv", required=True)
parser.add_argument("-i", "--interactions", help="interactions csv", required=True)
parser.add_argument("-o", "--output", help="path to the output directory to save output files", required=False, default=".")
parser.add_argument("-n", "--number", help="number of random samples for area plot", type=int, required=False, default=500)

if __name__ == "__main__":
    args = parser.parse_args()
   
    # parse and create all necessary objects for creating a network
    rates = parse_rates(args.rates)
    substrates = parse_substrates(args.substrates, rates)
    interactions = parse_interactions(args.interactions, rates, substrates)
    
    # create a new instance of a network
    network = Network("example", rates, interactions, substrates)
    
    # visualize network ordinary differential equations
    time = np.linspace(0, 1000, num=1001)
    derivatives = network.get_representation_dydt()
    for l, t in derivatives.items():
        print(f"{l} = {t}")
    print()

    

    # graph
    mean_y_regular = network.graph_distributions(time, args.number, normalize=True, initials=initials)
    indices = [list(network.substrates.keys()).index("PIP3"), list(network.substrates.keys()).index("pAKT"), list(network.substrates.keys()).index("pPTEN"), list(network.substrates.keys()).index("LPS")]
    plt.plot(range(len(mean_y_regular[:,indices[0]])), mean_y_regular[:,indices])
    plt.legend(["ctrl-PIP3", "ctrl-pAKT", "ctrl-pPTEN", "ctrl-LPS"])
    combo_figure.savefig(os.path.join(args.output, "./combo_fig.png"))