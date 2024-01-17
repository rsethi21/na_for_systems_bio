import pandas as pd
import argparse

from substrate import Substrate
from interaction import Interaction
from rate import Rate
from network import Network

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--substrates", help="substrates csv", required=True)
parser.add_argument("-r", "--rates", help="rates csv", required=True)
parser.add_argument("-i", "--interactions", help="interactions csv", required=True)
parser.add_argument("-f", "--fitting", help="fitting parameters json file", required=False,default=None)


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
        if parameters["bounds"] == None:
            del parameters["bounds"]
            del parameters["bounds_type"]
        else:
            if parameters["bounds_type"] != None:
                convert_type = parameters["bounds_type"]
            else:
                convert_type = int
            parameters["bounds"] = [convert_type(b) for b in parameters["bounds"].split(",")]
        rate_obj = Rate(**parameters)
        rates.append(rate_obj)
    
    for substrate_index, substrate in substrates_df.iterrows():
        parameters = substrate.to_dict()
        if parameters["k"] != None:
            parameters["k"] = [r for r in rates if r.identifier == parameters["k"]][0]
        if parameters["r"] != None:
            parameters["r"] = [r for r in rates if r.identifier == parameters["r"]][0]
        if parameters["trs"] != None:
            range_strings = parameters["trs"].split(";")
            time_ranges = []
            for range_string in range_strings:
                time_ranges.append([int(t) for t in range_string.split(",")])
            parameters["trs"] = time_ranges
        substrate_obj = Substrate(**parameters)
        substrates.append(substrate_obj)

    for interaction_index, interaction in interactions_df.iterrows():
        parameters = interaction.to_dict()
        parameters["affected"] = [s for s in substrates if s.identifier == parameters["affected"]][0]
        parameters["effector"] = [s for s in substrates if s.identifier == parameters["effector"]][0]
        parameters["rate"] = [r for r in rates if r.identifier == parameters["rate"]][0]
        if parameters["dissociation"] != None:
            parameters["dissociation"] = [r for r in rates if r.identifier == parameters["dissociation"]][0]
        if parameters["hill_coefficient"] != None:
            parameters["hill_coefficient"] = [r for r in rates if r.identifier == parameters["hill_coefficient"]][0]
        interaction_obj = Interaction(**parameters)
        interactions.append(interaction_obj)
    
    network = Network("example", rates, interactions, substrates)
    print(network.identifier)