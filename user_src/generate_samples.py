import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pdb
import json
import os
from tqdm import tqdm
import concurrent.futures as cf
from itertools import repeat, product
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
parser.add_argument("-m", "--multiprocess", help="number of threads", required=False, type=int, default=1)

def generate_output(random_input, labels, y0s, time, network, time_to_extract = 0, substrates_to_track=["Phagocytosis"]):
    for i, l in zip(random_input, labels):
        if network.substrates[l].time_ranges != None:
            network.substrates[l].max_value = i
        else:
            network.substrates[l].max_value = 0.0

    mean_y = network.graph_distributions(time, None, normalize=True, path=None, initials=y0s, verbose=False)
    indices_to_extract = [list(network.substrates.keys()).index(s) for s in substrates_to_track]
    output = []
    for t in time_to_extract:
        response = list(mean_y[t, indices_to_extract])
        output.extend(response)
    return output

if __name__ == "__main__":
    args = parser.parse_args()

    # parse and create all necesfsary objects for creating a network
    rates = parse_rates(args.rates)
    substrates = parse_substrates(args.substrates, rates)
    interactions = parse_interactions(args.interactions, rates, substrates)

    # create a new instance of a network
    network = Network("example", rates, interactions, substrates)

    # visualize network ordinary differential equations
    time = np.array([i for i in range(501)])
    derivatives = network.get_representation_dydt()

    # load in pretrained parameters if any
    if args.parameters != None:
        with open(args.parameters, "r") as fitted_params:
            parameters = json.load(fitted_params)
        network.set_parameters(list(parameters.values()), list(parameters.keys()))

    # random initial states
    y0s = []
    samples = 10
    for _ in range(samples):
        y0 = []
        for i, s in enumerate(network.substrates.values()):
            if s.type == "stimulus":
                y0.append(0.0)
            else:
                y0.append(2**np.random.randn())
        y0s.append(y0)
    
    # set appropriate conditions
    ranges = [[60, 120], [120, 180], [90, 210], [60, 240], None]
    conditions = list(product(ranges, ranges, ranges, ranges))
    external = ["LPS", "ATP", "LY294-002", "HDACi"]
    data = []
    for i_c, c in tqdm(enumerate(conditions), desc="Different Ranges", total=len(conditions)):
        for i, e in enumerate(external):
            if c[i] is not None:
                network.substrates[e].time_ranges = [c[i]]
            else:
                network.substrates[e].time_ranges = c[i]
        possible = []
        for p in c:
            if p == None:
                pass
            else:
                possible.extend(p)
        if len(possible) != 0:
            measure_times = [max(possible), max(possible) + 60, max(possible) + 120, max(possible)+180]
        else:
            measure_times = [240, 300, 360, 420]
        # random external stimuli
        random_inputs = 10**np.random.normal(0, 0.33, size=(args.number, 4))
        # generate samples
        labels = ["LPS", "ATP", "LY294-002", "HDACi"]
        subs = ["AKT", "pAKT", "PI3Ks", "PI3K", "GSK3B", "pGSK3B", "PTEN", "pPTEN", "PIP2", "PIP3", "P2Y12act", "P2Y12s", "Gio", "Phagocytosis"]
        subs2 = ["AKT_1hr", "pAKT_1hr", "PI3Ks_1hr", "PI3K_1hr", "GSK3B_1hr", "pGSK3B_1hr", "PTEN_1hr", "pPTEN_1hr", "PIP2_1hr", "PIP3_1hr", "P2Y12act_1hr", "P2Y12s_1hr", "Gio_1hr", "Phagocytosis_1hr"]
        subs3 = ["AKT_3hr", "pAKT_3hr", "PI3Ks_3hr", "PI3K_3hr", "GSK3B_3hr", "pGSK3B_3hr", "PTEN_3hr", "pPTEN_3hr", "PIP2_3hr", "PIP3_3hr", "P2Y12act_3hr", "P2Y12s_3hr", "Gio_3hr", "Phagocytosis_3hr"]
        subs4 = ["AKT_2hr", "pAKT_2hr", "PI3Ks_2hr", "PI3K_2hr", "GSK3B_2hr", "pGSK3B_2hr", "PTEN_2hr", "pPTEN_2hr", "PIP2_2hr", "PIP3_2hr", "P2Y12act_2hr", "P2Y12s_2hr", "Gio_2hr", "Phagocytosis_2hr"]
        for i, l in enumerate(labels):
            if network.substrates[l].time_ranges == None:
                random_inputs[:,i] = np.array([0.0 for _ in list(range(int(args.number)))])

        cols = subs.copy()
        cols.extend(subs2)
        cols.extend(subs4)
        cols.extend(subs3)
        with cf.ProcessPoolExecutor(max_workers=args.multiprocess) as executor:
            output = list(tqdm(executor.map(generate_output, random_inputs, repeat(labels), repeat(y0s), repeat(time), repeat(network), repeat(measure_times), repeat(subs)), total=len(random_inputs), desc=f"Generating Data with Range Pairs {c}"))
        output = np.array(output)
        output_df = pd.DataFrame(random_inputs, columns=labels)
        output_df[[f"{ext}_start" for ext in external]] = np.array([[c_i[0] if c_i != None else 0 for c_i in c]]*len(random_inputs))
        output_df[[f"{ext}_end" for ext in external]] = np.array([[c_i[1] if c_i != None else 0 for c_i in c]]*len(random_inputs))
        output_df["measure_time1"] = [measure_times[0] for _ in range(len(random_inputs))]
        output_df["measure_time2"] = [measure_times[1] for _ in range(len(random_inputs))]
        output_df["measure_time3"] = [measure_times[2] for _ in range(len(random_inputs))]
        output_df["measure_time4"] = [measure_times[3] for _ in range(len(random_inputs))]

        output_df[cols] = output
        data.append(output_df)
    final_df = data[0]
    for df in data[1:]:
        try:
            final_df = final_df._append(df, ignore_index=True)
        except:
            final_df = final_df.append(df, ignore_index=True)
    final_df.to_csv(os.path.join(args.output, "data.csv"))
