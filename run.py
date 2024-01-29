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

    # load in pretrained parameters if any
    if args.parameters != None:
        with open(args.parameters, "r") as fitted_params:
            parameters = json.load(fitted_params)
        network.set_parameters(list(parameters.values()), list(parameters.keys()))

    # train model if needed
    if args.fitting != None or args.arguments != None:
        with open(args.fitting, "r") as file:
            fit_dictionary = json.load(file)
        with open(args.arguments, "r") as argument_file:
            argues = json.load(argument_file)
        network.fit(fit_dictionary, time, argues, normalize=True, mlp=args.multi)
        with open(os.path.join(args.output, "fitted_params.json"), "w") as out_file:
            json.dump({i: r.value for i, r in network.parameters.items()}, out_file)
    
    # create area plots for the network
    y, f = network.graph_distributions(time, args.number, substrates_to_plot=["PI3K", "pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi"], normalize=True, path=os.path.join(args.output, "figure.png"), output_figure=True)

    # plot fitting data to visually evaluate the fit
    plt.figure(f)
    times = [60, 240, 245, 250, 260, 270, 300]
    times2 = [420, 420, 420, 420]
    values = [0.4, 0.4, 0.6, 1.1]
    names = ["pAKT", "PIP3", "pPTEN", "GSK3B"]
    phagocytosis_ys = [1.0, 0.3235695, 0.35537492, 0.33414908, 0.34072799, 0.46634439, 0.44325203]
    
    plt.scatter(times, phagocytosis_ys, s=2, color=network.colors[list(network.substrates.keys()).index("Phagocytosis")])
    plt.scatter(times2, values, s=2, color=[network.colors[list(network.substrates.keys()).index(n)] for n in names])
    f.savefig(os.path.join(args.output, "figure_w_data.png"))
