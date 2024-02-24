import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pdb
import json
import os
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
    for l, t in derivatives.items():
        print(f"{l} = {t}")
    print()

    # open fitting data
    with open(args.fitting, "r") as file:
        fit_dictionary = json.load(file)

    '''
    # check unfitted rates error
    print("Randomly Generated Error")
    network.set_parameters(np.random.rand(len([p for p in network.parameters.values() if p.fixed==False])), [p.identifier for p in network.parameters.values() if p.fixed==False])
    original, original_fig = network.graph_distributions(time, args.number, substrates_to_plot=["PI3K", "pAKT", "pPTEN", "Phagocytosis", "TNFa", "LPS", "HDACi", "GSK3B", "LY294-002"], normalize=True, path=os.path.join(args.output, "figure_0.png"), output_figure=True, axis=None)
    print(error(original, fit_dictionary, list(network.substrates.keys())))
    plt.figure(original_fig)
    plt.plot(420, 3.17, "ko")
    original_fig.savefig(os.path.join(args.output, "figure_0.png"))
    # plt.close(original_fig)
    print()
    '''
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
        print("Training...")
        network.fit(fit_dictionary, time, argues, number=10, normalize=True, mlp=args.multi)
        print("Saving learned parameters...")
        with open(os.path.join(args.output, "fitted_params_tnfa.json"), "w") as out_file:
            json.dump({i: r.value for i, r in network.parameters.items()}, out_file)
    
    # check post-training training error
    print("TNF-alpha Error")
    mean_y, f = network.graph_distributions(time, args.number, substrates_to_plot=["PI3K", "pAKT", "pPTEN", "Phagocytosis", "TNFa", "LPS", "HDACi", "GSK3B", "LY294-002"], normalize=True, path=os.path.join(args.output, "figure_3.png"), output_figure=True)
    plt.figure(f)
    plt.plot(420, 3.17, "ko")
    f.savefig(os.path.join(args.output, "figure_3.png"))
    print(error(mean_y, fit_dictionary, list(network.substrates.keys())))
    # plt.figure(f)
    # plt.plot(420, 0.6, "ko")
    # f.savefig(os.path.join(args.output, "figure_3.png"))
    plt.close(f)
    print()
