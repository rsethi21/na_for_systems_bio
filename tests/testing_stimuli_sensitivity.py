import sys
sys.path.append(f"{sys.path[0]}/..")

from src.network import Network
import unittest
import argparse
import pyyaml

class TestLPSEffectOnPTEN(unittest.TestCase):
    def runTest(self):
        # [AKT, pAKT, PTEN, pPTEN, PIP2, PIP3, PI3K, PI3Ks, GSK3B, pGSK3B, TNFa, Phagocytosis, P2Y12act, P2Y12s, Gio]
        initial_conditions = [[1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0]]
        stimuli_range = [[60, 120]]
        stimuli_amt = 1.0
        time = None

class TestLPSEffectOnTotalPTEN(unittest.TestCase):
    def runTest(self):
        # [AKT, pAKT, PTEN, pPTEN, PIP2, PIP3, PI3K, PI3Ks, GSK3B, pGSK3B, TNFa, Phagocytosis, P2Y12act, P2Y12s, Gio]
        initial_conditions = [[1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0]]
        stimuli_range = [[60, 120]]
        stimuli_amt = 1.0
        time = None