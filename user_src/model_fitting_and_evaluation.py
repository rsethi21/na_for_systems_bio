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
# parser.add_argument("-h", "--hyperparameters", help="path to fitting hyperparameters json file", required=True, default=None)
parser.add_argument("-a", "--arguments", help="fitting algorithm arguments json file", required=False, default=None)
parser.add_argument("-o", "--output", help="path to the output directory to save output files", required=False, default=".")
parser.add_argument("-p", "--parameters", help="pretrained parameters json file", required=False, default=None)
parser.add_argument("-m", "--multi", help="number of processers to use for fitting", type=int, required=False, default=1)
parser.add_argument("-n", "--number", help="number of random samples for area plot", type=int, required=False, default=10)
parser.add_argument("-l", "--less", help="scaling factors to apply to parameters", required=False, default=None)
parser.add_argument("-e", "--scale", help="scaling factors to apply to parameters", required=False, default=None)

if __name__ == "__main__":
    args = parser.parse_args()
   
    # parse and create all necessary objects for creating a network
    rates = parse_rates(args.rates)
    substrates = parse_substrates(args.substrates, rates)
    interactions = parse_interactions(args.interactions, rates, substrates)
    
    # create a new instance of a network
    network = Network("example", rates, interactions, substrates)
    
    # # fitting hyperparamters
    # with open(args.hyperparameters, "r") as hp_file:
    #     hps = json.load(hp_file)

    # visualize network ordinary differential equations
    # time = np.array([i for i in range(hps["end_time"])])
    time = np.array([i for i in range(1001)])
    derivatives = network.get_representation_dydt()
    for l, t in derivatives.items():
        print(f"{l} = {t}")
    print()

    # open fitting data
    if args.fitting != None:
        with open(args.fitting, "r") as file:
            fit_dictionary = json.load(file)

    # check unfitted rates error
    # print("Randomly Generated Error")
    
    # random_parameters = np.random.rand(len([p for p in network.parameters.values() if p.fixed==False]))
    # network.set_parameters(random_parameters, [p.identifier for p in network.parameters.values() if p.fixed==False])
    original = network.graph_distributions(time, args.number, normalize=True, substrates_to_plot=["pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi", "GSK3B", "LY294-002"], path=os.path.join(args.output, "figure_0_fit.png"), output_figure=False)
    
    # print(error(original, fit_dictionary, list(network.substrates.keys())))
    # print("BV2 --> Primary")
    # print(error(original, {"Phagocytosis": {420: 0.6}}, list(network.substrates.keys())))
    # # plt.figure(original_fig)
    # plt.plot(420, 0.6, "ko")
    # original_fig.savefig(os.path.join(args.output, "figure_0.png"))
    # plt.close(original_fig)
    #print()

    # load in pretrained parameters if any
    if args.parameters != None:
        with open(args.parameters, "r") as fitted_params:
            parameters = json.load(fitted_params)
        network.set_parameters(list(parameters.values()), list(parameters.keys()))

    if args.less != None:
        with open(args.less) as less_params:
            scales = json.load(less_params)
        network.apply_scaling(list(scales.values(), list(scales.names())))

    # train model if needed
    if args.fitting != None and args.arguments != None:
        with open(args.fitting, "r") as file:
            fit_dictionary = json.load(file)
        with open(args.arguments, "r") as argument_file:
            argues = json.load(argument_file)
        if args.scale != None:
            network.scale(fit_dictionary, time, argues, scaling_bounds=[1,3], number=args.number, normalize=True, mlp=args.multi)
        else:
            network.fit(fit_dictionary, time, argues, number=args.number, normalize=True, mlp=args.multi)
        with open(os.path.join(args.output, "fitted_params_scaled.json"), "w") as out_file:
            json.dump({i: r.value for i, r in network.parameters.items()}, out_file)

    # check post-training training error
    print("LPS, HDACi, LY294-002 Error")
    mean_y = network.graph_distributions(time, args.number, substrates_to_plot=["pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi", "GSK3B", "LY294-002"], normalize=True, path=os.path.join(args.output, "figure_3_scaled.png"), output_figure=False)
    # print(error(mean_y, fit_dictionary, list(network.substrates.keys())))
    # print(error(mean_y, {"Phagocytosis": {420: 0.6}}, list(network.substrates.keys())))
    # plt.figure(f)
    # plt.plot(420, 0.6, "ko")
    # f.savefig(os.path.join(args.output, "figure_3.png"))
    # plt.close(f)

    # check post-training validation error
    print("LPS and HDACi Error")
    network.substrates["LY294-002"].time_ranges = None
    network.substrates["LY294-002"].max_value = 0
    mean_lps_hdaci_y = network.graph_distributions(time, args.number, substrates_to_plot=["pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi", "GSK3B", "LY294-002"], normalize=True, path=os.path.join(args.output, "figure_2_scaled.png"), output_figure=False)
    # print(error(mean_lps_hdaci_y, lps_hdaci, list(network.substrates.keys())))
    # print("BV2 --> Primary")
    # print(error(mean_lps_hdaci_y, {"Phagocytosis": {420: 2.7}}, list(network.substrates.keys())))
    # plt.figure(lh_fig)
    # plt.plot(420, 2.7, "ko")
    # axs[1].set_figure(lh_fig)
    # lh_fig.savefig(os.path.join(args.output, "figure_2.png"))
    # plt.close(lh_fig)

    print("LPS, HDACi, shRNA Error")
    network.substrates["LY294-002"].time_ranges = None
    network.substrates["LY294-002"].max_value = 0
    network.substrates["shRNA_GSK3B"].max_value = 1.0
    network.substrates["shRNA_GSK3B"].time_ranges = [[60,240]]
    mean_lps_hdaci_y = network.graph_distributions(time, args.number, substrates_to_plot=["pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi", "GSK3B", "LY294-002"], normalize=True, path=os.path.join(args.output, "figure_4_scaled.png"), output_figure=False)

    print("LPS Error")
    network.substrates["HDACi"].time_ranges = None
    network.substrates["HDACi"].max_value = 0
    network.substrates["shRNA_GSK3B"].max_value = 0
    network.substrates["shRNA_GSK3B"].time_ranges = None
    mean_lps_y = network.graph_distributions(time, args.number, substrates_to_plot=["pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi", "GSK3B", "LY294-002"], normalize=True, path=os.path.join(args.output, "figure_1_scaled.png"), output_figure = False)
    # print(error(mean_lps_y, lps, list(network.substrates.keys())))
    # print("BV2 --> Primary")
    # print(error(mean_lps_y, {"Phagocytosis": {240: 0.5}}, list(network.substrates.keys())))
    # plt.figure(l_f)
    # plt.plot(420, 0.5, "ko")
    # axs[2].set_figure(l_f)
    # l_f.savefig(os.path.join(args.output, "figure_1.png"))
    # plt.close(l_f)