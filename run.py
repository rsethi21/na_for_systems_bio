import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pdb
import json
import os

from substrate import Substrate
from interaction import Interaction
from rate import Rate
from network import Network
from parse import parse_rates, parse_substrates, parse_interactions
from score import error

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--substrates", help="substrates csv", required=True)
parser.add_argument("-r", "--rates", help="rates csv", required=True)
parser.add_argument("-i", "--interactions", help="interactions csv", required=True)
parser.add_argument("-f", "--fitting", help="fitting parameters json file", required=False, default=None)
parser.add_argument("-a", "--arguments", help="fitting algorithm arguments json file", required=False, default=None)
parser.add_argument("-o", "--output", help="path to the output directory to save output files", required=False, default=".")
parser.add_argument("-p", "--parameters", help="pretrained parameters json file", required=False, default=None)
parser.add_argument("-m", "--multi", help="number of processers to use for fitting", type=int, required=False, default=1)
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
    time = np.array([i for i in range(500)])
    derivatives = network.get_representation_dydt()
    print(derivatives)
    print()

    # open fitting data
    with open(args.fitting, "r") as file:
        fit_dictionary = json.load(file)

    # check unfitted rates error
    print("Randomly Generated Error")
    network.set_parameters(np.random.rand(len([p for p in network.parameters.values() if p.fixed==False])), [p.identifier for p in network.parameters.values() if p.fixed==False])
    original = network.graph_distributions(time, args.number, normalize=True, path=None)
    print(error(original, fit_dictionary, list(network.substrates.keys())))
    print()

    # load in pretrained parameters if any
    if args.parameters != None:
        with open(args.parameters, "r") as fitted_params:
            parameters = json.load(fitted_params)
        network.set_parameters(list(parameters.values()), list(parameters.keys()))

    # train model if needed
    if args.fitting != None and args.arguments != None:
        with open(args.fitting, "r") as file:
            fit_dictionary = json.load(file)
        with open(args.arguments, "r") as argument_file:
            argues = json.load(argument_file)
        network.fit(fit_dictionary, time, argues, normalize=True, mlp=args.multi)
        with open(os.path.join(args.output, "fitted_params.json"), "w") as out_file:
            json.dump({i: r.value for i, r in network.parameters.items()}, out_file)
    
    # check post-training training error
    print("LPS, HDACi, LY294-002 Error")
    mean_y, f = network.graph_distributions(time, args.number, substrates_to_plot=["PI3K", "pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi", "GSK3B", "LY294-002"], normalize=True, path=os.path.join(args.output, "figure.png"), output_figure=True)
    print(error(mean_y, fit_dictionary, list(network.substrates.keys())))
    print()

    # check post-training validation error
    print("LPS and HDACi Error")
    network.substrates["LY294-002"].time_ranges = None
    lps_hdaci = {
        "pAKT": {"410": 0.8, "420": 0.8,},
        "pPTEN": {"410": 1.4, "420": 1.4,},
        "PIP3": {"410": 0.8, "420": 0.8,},
        "GSK3B": {"410": 1.1, "420": 1.1,}
        }
    mean_lps_hdaci_y = network.graph_distributions(time, args.number, normalize=True, path=None)
    print(error(mean_lps_hdaci_y, lps_hdaci, list(network.substrates.keys())))
    print()

    print("LPS Error")
    network.substrates["HDACi"].time_ranges = None
    lps = {
        "pAKT": {"410": 0.2, "420": 0.2,},
        "pPTEN": {"410": 0.6, "420": 0.6,},
        "PIP3": {"410": 0.2, "420": 0.2,},
        "GSK3B": {"410": 0.7, "420": 0.7,},
        "Phagocytosis": {"360": 0.80413556, "420": 1.04238109}
        }
    mean_lps_y = network.graph_distributions(time, args.number, normalize=True, path=None)
    print(error(mean_lps_y, lps, list(network.substrates.keys())))
    print()
