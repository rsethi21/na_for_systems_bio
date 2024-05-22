import sys
sys.path.append(f"{sys.path[0]}/..")

from src.network import Network
from src.parse import *
import unittest
import numpy as np
import argparse
import json
from itertools import product as P
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--substrates", help="substrates csv", required=True)
parser.add_argument("-r", "--rates", help="rates csv", required=True)
parser.add_argument("-i", "--interactions", help="interactions csv", required=True)
parser.add_argument("-p", "--parameters", help="pretrained parameters json file", required=False, default=None)
parser.add_argument("-t", "--test_configs", help="test configuration json file", required=True)

def setUp(self, network, substrate_of_interest, stimuli_of_interest, initial_conditions, stimuli_range, stimuli_amts, time, expectations):
    self.network = network
    self.substrate_index = substrate_of_interest
    self.stimuli_index = stimuli_of_interest
    self.initial_conditions = initial_conditions
    self.stimuli_range = stimuli_range
    self.stimuli_amts = stimuli_amts
    self.time = time
    self.expectations = [True if expectation == "T" else False for expectation in expectations]

class TestStimuliEffectOnSubstrates(unittest.TestCase):
    def __init__(self, testName, network, substrate_of_interests, stimuli_of_interest, initial_conditions, stimuli_range, stimuli_amts, time, expectations):
        super().__init__(testName)
        self.network = network
        self.substrate_of_interests = substrate_of_interests
        self.substrate_indices = [list(network.substrates.keys()).index(i) for i in substrate_of_interests]
        self.stimuli_index = stimuli_of_interest
        self.initial_conditions = initial_conditions
        self.stimuli_range = stimuli_range
        self.stimuli_amts = stimuli_amts
        self.time = time
        self.expectations = expectations

    def StimuliEffectOnSubstrates(self):
        # [AKT, pAKT, PTEN, pPTEN, PIP2, PIP3, PI3K, PI3Ks, GSK3B, pGSK3B, TNFa, Phagocytosis, P2Y12act, P2Y12s, Gio]
        self.network.substrates[self.stimuli_index].time_ranges = [self.stimuli_range]
        self.network.set_initials(self.initial_conditions)
        max_values = []
        min_values = []
        for amt in self.stimuli_amts:
            self.network.substrates[self.stimuli_index].max_value = amt
            trajectory = self.network.graph(np.linspace(self.time[0], self.time[1], num=self.time[1]+1))
            max_values.append(max(trajectory[self.stimuli_range[0]-5:self.stimuli_range[1]+5,self.substrate_indices[0]]))
            min_values.append(min(trajectory[self.stimuli_range[0]-5:self.stimuli_range[1]+5,self.substrate_indices[0]]))

    def ConstatntTotalSubstrateUponStimuli(self):
        # [AKT, pAKT, PTEN, pPTEN, PIP2, PIP3, PI3K, PI3Ks, GSK3B, pGSK3B, TNFa, Phagocytosis, P2Y12act, P2Y12s, Gio]
        self.network.substrates[self.stimuli_index].max_value = self.stimuli_amts[0]
        self.network.substrates[self.stimuli_index].time_ranges = [self.stimuli_range]
        self.network.set_initials(self.initial_conditions)
        trajectory = self.network.graph(np.linspace(self.time[0], self.time[1], num=self.time[1]+1), substrates_to_plot=self.substrate_of_interests, path="./test.png")
        sums = trajectory[:,self.substrate_indices[0]] + trajectory[:,self.substrate_indices[1]]
        check = sum(sums[self.stimuli_range[0]-5:self.stimuli_range[1]+5])/len(sums[self.stimuli_range[0]-5:self.stimuli_range[1]+5])
        print(f"Checking if {check} ~= {sums[0]}...")
        self.assertAlmostEqual(sums[0], check, f"Total substrate upon stimulation by {self.stimuli_index} is not constant.")
if __name__ == "__main__":

    # open arguments
    args = parser.parse_args()
    # parse and create all necessary objects for creating a network
    rates = parse_rates(args.rates)
    substrates = parse_substrates(args.substrates, rates)
    interactions = parse_interactions(args.interactions, rates, substrates)
    # create a new instance of a network
    network = Network("example", rates, interactions, substrates)
    if args.parameters != None:
        with open(args.parameters, "r") as fitted_params:
            parameters = json.load(fitted_params)
        network.set_parameters(list(parameters.values()), list(parameters.keys()))
    with open(args.test_configs, "r") as configurations:
        tc = json.load(configurations)

    suite = unittest.TestSuite()
    suite.addTest(TestStimuliEffectOnSubstrates("ConstatntTotalSubstrateUponStimuli", network, tc["substrate_of_interests"], tc["stimuli_of_interest"], tc["initial_conditions"], tc["stimuli_range"], tc["stimuli_amts"], tc["time"], tc["expectations"]))
    unittest.TextTestRunner(verbosity=2).run(suite)