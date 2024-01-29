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

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--substrates", help="substrates csv", required=True)
parser.add_argument("-r", "--rates", help="rates csv", required=True)
parser.add_argument("-i", "--interactions", help="interactions csv", required=True)
parser.add_argument("-f", "--fitting", help="fitting parameters json file", required=False, default=None)
parser.add_argument("-a", "--arguments", help="fitting algorithm arguments json file", required=False, default=None)
parser.add_argument("-o", "--output", help="path to the output directory to save output files", required=False, default=".")
parser.add_argument("-p", "--parameters", help="pretrained parameters json file", required=False, default=None)
parser.add_argument("-m", "--multi", help="number of processers to use for fitting", required=False, default=1)

if __name__ == "__main__":
    args = parser.parse_args()
    
    rates_df = pd.read_csv(args.rates)
    substrates_df = pd.read_csv(args.substrates)
    interactions_df = pd.read_csv(args.interactions)

    rates = []
    substrates = []
    interactions = []
    
    for rate_index, rate in rates_df.iterrows():
        parameters = rate.to_dict()
        if pd.isna(parameters["bounds"]):
            del parameters["bounds"]
            del parameters["bounds_type"]
        else:
            if not pd.isna(parameters["bounds_type"]):
                convert_type = parameters["bounds_type"]
                if convert_type == "int":
                    convert_type = int
                elif convert_type == "real":
                    convert_type = float
            else:
                convert_type = float
            parameters["bounds"] = [convert_type(b) for b in parameters["bounds"].split(",")]
        rate_obj = Rate(**parameters)
        rates.append(rate_obj)

    for substrate_index, substrate in substrates_df.iterrows():
        parameters = substrate.to_dict()
        if not pd.isna(parameters["k"]):
            parameters["k"] = [r for r in rates if r.identifier == parameters["k"]][0]
        else:
            del parameters["k"]
        if not pd.isna(parameters["r"]):
            parameters["r"] = [r for r in rates if r.identifier == parameters["r"]][0]
        else:
            del parameters["r"]
        if not pd.isna(parameters["trs"]):
            range_strings = parameters["trs"].split(";")
            time_ranges = []
            for range_string in range_strings:
                time_ranges.append([int(t) for t in range_string.split(",")])
            parameters["trs"] = time_ranges
        else:
            del parameters["trs"]
        if not pd.isna(parameters["other"]):
            pass
        else:
            del parameters["other"]
        substrate_obj = Substrate(**parameters)
        substrates.append(substrate_obj)

    for interaction_index, interaction in interactions_df.iterrows():
        parameters = interaction.to_dict()
        parameters["affected"] = [s for s in substrates if s.identifier == parameters["affected"]][0]
        parameters["effector"] = [s for s in substrates if s.identifier == parameters["effector"]][0]
        parameters["rate"] = [r for r in rates if r.identifier == parameters["rate"]][0]
        if not pd.isna(parameters["dissociation"]):
            parameters["dissociation"] = [r for r in rates if r.identifier == parameters["dissociation"]][0]
        else:
            del parameters["dissociation"]
        if not pd.isna(parameters["hill_coefficient"]):
            parameters["hill_coefficient"] = [r for r in rates if r.identifier == parameters["hill_coefficient"]][0]
        else:
            del parameters["hill_coefficient"]
        if not pd.isna(parameters["effect"]):
            pass
        else:
            del parameters["effect"]
        interaction_obj = Interaction(**parameters)
        interactions.append(interaction_obj)
    
    network = Network("example", rates, interactions, substrates)
    time = np.array([i for i in range(500)])
    derivatives = network.get_representation_dydt()
    print(derivatives)

    if args.parameters == None and args.fitting != None and args.arguments != None:
        with open(args.fitting, "r") as file:
            fit_dictionary = json.load(file)
        with open(args.arguments, "r") as argument_file:
            argues = json.load(argument_file)
        network.fit(fit_dictionary, time, argues, normalize=True, mlp=args.multi)
        with open(os.path.join(args.output, "fitted_params.json"), "w") as out_file:
            json.dump({i: r.value for i, r in network.parameters.items()}, out_file)
    else:
        with open(args.parameters, "r") as fitted_params:
            parameters = json.load(fitted_params)
        network.set_parameters(list(parameters.values()), list(parameters.keys()))
    
    y, f = network.graph_distributions(time, 500, substrates_to_plot=["PI3K", "pAKT", "pPTEN", "Phagocytosis", "LPS", "HDACi"], normalize=True, path=os.path.join(args.output, "figure.png"), output_figure=True)

    plt.figure(f)
    times = [60, 240, 245, 250, 260, 270, 300, 360, 420]
    phagocytosis_ys = [1.0, 0.3235695, 0.35537492, 0.33414908, 0.34072799, 0.46634439, 0.44325203, 0.80413556, 1.04238109]
    plt.scatter([60, 240, 245, 250, 260, 270, 300, 360, 420], phagocytosis_ys, s=1)
    f.savefig(os.path.join(args.output, "figure_w_data.png"))
