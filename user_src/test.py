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
parser.add_argument("-d", "--data", help="fit data json", required=True)
parser.add_argument("-a", "--arguments", help="fitting algorithm arguments", required=True)
parser.add_argument("-o", "--output", help="path to the output directory to save output files", required=False, default=".")

if __name__ == "__main__":
    args = parser.parse_args()
   
    # parse and create all necessary objects for creating a network
    rates = parse_rates(args.rates)
    substrates = parse_substrates(args.substrates, rates)
    interactions = parse_interactions(args.interactions, rates, substrates)
    
    # create a new instance of a network
    network = Network("example", rates, interactions, substrates)
    
    # visualize network ordinary differential equations
    time = np.linspace(0, 250, num=250)
    derivatives = network.get_representation_dydt()
    for l, t in derivatives.items():
        print(f"{l} = {t}")
    print()

    # original fit
    network.graph_distributions(time, 10, normalize=False, substrates_to_plot=["A", "B", "C", "S1", "S2"], path=os.path.join(args.output, "original.png"))

    # fit model to data
    network.fit(json.load(open(args.data, "r")), time, json.load(open(args.arguments, "r")), number=10, normalize=True, mlp=12, plots_path=args.output)