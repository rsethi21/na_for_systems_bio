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
    print("LPS, HDACi, LY294-002 Error")
    mean_y3, min_y3, max_y3, f3 = network.graph_distributions(time, args.number, substrates_to_plot=["LPS", "HDACi", "LY294-002"], path=os.path.join(args.output, "figure_3.png"), output_figure=True, return_min_max=True)
    plt.figure(f3)
    plt.plot(420, 0.938, "ko")
    # plt.fill_between(time, np.log10(min_y3[:,list(network.substrates.keys()).index("TNFa")]), np.log10(max_y3[:,list(network.substrates.keys()).index("TNFa")]), "b--", alpha=0.2)
    plt.plot(time, mean_y3[:,list(network.substrates.keys()).index("TNFa")], "b-", label="TNFa", linewidth=1.0)
    f3.savefig(os.path.join(args.output, "figure_3.png"))
    print(error(mean_y3, fit_dictionary, list(network.substrates.keys())))
    plt.legend()
    plt.close(f3)
    print()

    network.substrates["LY294-002"].time_ranges = None
    network.substrates["LY294-002"].max_value = 0.0
    print("LPS, HDACi Error")
    mean_y2, min_y2, max_y2, f2 = network.graph_distributions(time, args.number, substrates_to_plot=["LPS", "HDACi", "LY294-002"], path=os.path.join(args.output, "figure_2.png"), output_figure=True, return_min_max=True)
    plt.figure(f2)
    plt.plot(420, 0.469, "ko")
    # plt.fill_between(time, np.log10(min_y2[:,list(network.substrates.keys()).index("TNFa")]), np.log10(max_y2[:,list(network.substrates.keys()).index("TNFa")]), "b--", alpha=0.2)
    plt.plot(time, mean_y2[:,list(network.substrates.keys()).index("TNFa")], "b-", label="TNFa", linewidth=1.0)
    f2.savefig(os.path.join(args.output, "figure_2.png"))
    print(error(mean_y2, {"TNFa": {"0": 0.0, "30": 0.0, "410": 0.469, "420": 0.469}}, list(network.substrates.keys())))
    plt.legend()
    plt.close(f2)
    print()

    network.substrates["HDACi"].time_ranges = None
    network.substrates["HDACi"].max_value = 0.0
    print("LPS Error")
    mean_y, min_y, max_y, f = network.graph_distributions(time, args.number, substrates_to_plot=["LPS", "HDACi", "LY294-002"], path=os.path.join(args.output, "figure.png"), output_figure=True, return_min_max=True)
    plt.figure(f)
    plt.plot(420, 1.0, "ko")
    # plt.fill_between(time, np.log10(min_y[:,list(network.substrates.keys()).index("TNFa")]), np.log10(max_y[:,list(network.substrates.keys()).index("TNFa")]), "b--", alpha=0.2)
    plt.plot(time, mean_y[:,list(network.substrates.keys()).index("TNFa")], "b-", label="TNFa", linewidth=1.0)
    f.savefig(os.path.join(args.output, "figure.png"))
    print(error(mean_y, {"TNFa": {"0": 1.0, "30": 1.0, "410": 1.0, "420": 1.0}}, list(network.substrates.keys())))
    plt.legend()
    plt.close(f)
    print()
