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
   
    main_figure, axs = plt.subplots(1, 4, figsize=(24, 9), sharex=True, sharey=True)

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

    # check unfitted rates error
    print("Randomly Generated Error")
    network.set_parameters(np.random.rand(len([p for p in network.parameters.values() if p.fixed==False])), [p.identifier for p in network.parameters.values() if p.fixed==False])
    original, original_fig = network.graph_distributions(time, args.number, substrates_to_plot=["PI3K", "pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi", "GSK3B", "LY294-002"], normalize=True, path=os.path.join(args.output, "figure_0.png"), output_figure=True, axis=None)
    print(error(original, fit_dictionary, list(network.substrates.keys())))
    print("BV2 --> Primary")
    print(error(original, {"Phagocytosis": {420: 0.6}}, list(network.substrates.keys())))
    # plt.figure(original_fig)
    # plt.plot(420, 0.6, "ko")
    # original_fig.savefig(os.path.join(args.output, "figure_0.png"))
    # plt.close(original_fig)
    #print()

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
        network.fit(fit_dictionary, time, argues, number=5, normalize=True, mlp=args.multi)
        with open(os.path.join(args.output, "fitted_params.json"), "w") as out_file:
            json.dump({i: r.value for i, r in network.parameters.items()}, out_file)
    
    # check post-training training error
    print("LPS, HDACi, LY294-002 Error")
    mean_y = network.graph_distributions(time, args.number, substrates_to_plot=["PI3K", "pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi", "GSK3B", "LY294-002"], normalize=True, path=os.path.join(args.output, "figure_3.png"), output_figure=True, axis=axs[0])
    print(error(mean_y, fit_dictionary, list(network.substrates.keys())))
    print("BV2 --> Primary")
    print(error(mean_y, {"Phagocytosis": {420: 0.6}}, list(network.substrates.keys())))
    # plt.figure(f)
    # plt.plot(420, 0.6, "ko")
    # f.savefig(os.path.join(args.output, "figure_3.png"))
    # plt.close(f)
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
    mean_lps_hdaci_y = network.graph_distributions(time, args.number, substrates_to_plot=["PI3K", "pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi", "GSK3B", "LY294-002"], normalize=True, path=os.path.join(args.output, "figure_2.png"), output_figure=True, axis=axs[1])
    print(error(mean_lps_hdaci_y, lps_hdaci, list(network.substrates.keys())))
    print("BV2 --> Primary")
    print(error(mean_lps_hdaci_y, {"Phagocytosis": {420: 2.7}}, list(network.substrates.keys())))
    # plt.figure(lh_fig)
    # plt.plot(420, 2.7, "ko")
    # axs[1].set_figure(lh_fig)
    # lh_fig.savefig(os.path.join(args.output, "figure_2.png"))
    # plt.close(lh_fig)
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
    mean_lps_y = network.graph_distributions(time, args.number, substrates_to_plot=["PI3K", "pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi", "GSK3B", "LY294-002", "ATP"], normalize=True, path=os.path.join(args.output, "figure_1.png"), output_figure = True, axis=axs[3])
    print(error(mean_lps_y, lps, list(network.substrates.keys())))
    print("BV2 --> Primary")
    print(error(mean_lps_y, {"Phagocytosis": {240: 0.5}}, list(network.substrates.keys())))
    # plt.figure(l_f)
    # plt.plot(420, 0.5, "ko")
    # axs[2].set_figure(l_f)
    # l_f.savefig(os.path.join(args.output, "figure_1.png"))
    # plt.close(l_f)
    print()

    network.substrates["ATP"].time_ranges = [[300, 420]]
    network.substrates["ATP"].max_value = 2.0
    mean_atp_y = network.graph_distributions(time, args.number, substrates_to_plot=["PI3K", "pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi", "GSK3B", "LY294-002", "ATP"], normalize=True, path=os.path.join(args.output, "figure_11.png"), output_figure = True, axis=axs[2])

    axs[0].set_title("LPS,HDACi,LY294-002")
    axs[1].set_title("LPS,HDACi")
    axs[3].set_title("LPS")
    axs[2].set_title("LPS,ATP")
    axs[1].set_xlabel("Time")
    axs[2].set_xlabel("(mins)")
    axs[0].set_ylabel("Concentration (AU)")
    axs[3].legend(loc="upper right", fontsize=8)
    main_figure.savefig(os.path.join(args.output, "birdseye.png"))
