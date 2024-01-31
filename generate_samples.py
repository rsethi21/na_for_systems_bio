import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pdb
import json
import os
from tqdm import tqdm
import concurrent.futures as cf
from itertools import repeat

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
        network.substrates[l].max_value = i
    mean_y = network.graph_distributions(time, None, normalize=True, path=None, initials=y0s, verbose=False)
    indices_to_extract = [list(network.substrates.keys()).index(s) for s in substrates_to_track]
    response = mean_y[time_to_extract, indices_to_extract]
    return response

if __name__ == "__main__":
    args = parser.parse_args()

    # parse and create all necesfsary objects for creating a network
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

    # load in pretrained parameters if any
    if args.parameters != None:
        with open(args.parameters, "r") as fitted_params:
            parameters = json.load(fitted_params)
        network.set_parameters(list(parameters.values()), list(parameters.keys()))

    # set appropriate conditions
    conditions = ["LPS", "HDACi", "LY294-002"] 
    for c in conditions:
        network.substrates[c].time_ranges = [[60, 240]]

    # generate data
    y0s = []
    samples = 50
    for _ in tqdm(range(samples),desc="Generating Random Initial",total=samples):
        y0 = []
        for i, s in enumerate(network.substrates.values()):
            if s.type == "stimulus":
                y0.append(0.0)
            else:
                y0.append(2**np.random.randn())
        y0s.append(y0)
    labels = ["LPS", "HDACi", "LY294-002"]
    random_inputs = 2**np.random.randn(args.number, 3)
    subs = ["pAKT", "PI3K", "GSK3B", "PTEN"]
    with cf.ProcessPoolExecutor(max_workers=args.multiprocess) as executor:
        output = list(tqdm(executor.map(generate_output, random_inputs, repeat(labels), repeat(y0s), repeat(time), repeat(network), repeat(420), repeat(subs)), total=len(random_inputs), desc="Generating Data"))
    output = np.array(output)
    output_df = pd.DataFrame(random_inputs, columns=labels)
    output_df[subs] = output
    output_df.to_csv(os.path.join(args.output, "data.csv"))
