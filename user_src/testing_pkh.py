import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pdb
import json
import os
from itertools import product as P
from tqdm import tqdm

from src.substrate import Substrate
from src.interaction import Interaction
from src.rate import Rate
from src.network import Network
from src.parse import parse_rates, parse_substrates, parse_interactions
from src.score import error

def runSim(
#    ranges = [[60,120], [60,180], [60,240]],
    amtsAtp = [1.0,1.0],
    ranges = [[60,120]],
    amtsLps = [0.0,1.0],    # need to be same length as amtsAtp
    params = None,
    tag = None
):

    # parse and create all necessary objects for creating a network
    #path="na_for_systems_bio"
    path="pten/input/input_atp_data/"
    rates = parse_rates(path+"/rates.csv")
    substrates = parse_substrates(path+"/substrates.csv", rates)
    interactions = parse_interactions(path+"/interactions.csv", rates, substrates)

    # create a new instance of a network
    network = Network("example", rates, interactions, substrates)
    # needed for enforcing parameters 
    for r in network.rates:
      r.fixed=False

    # visualize network ordinary differential equations
    time = np.linspace(0, 500, num=501)
    derivatives = network.get_representation_dydt()
    for l, t in derivatives.items():
        print(f"{l} = {t}")
    print()

    if params is not None:
      print("Applying parameters", params.keys())
      network.set_parameters(list(params.values()), list(params.keys()))
    
    combos = list(P(amtsAtp, ranges))
    fnames=[]
    for i, combo in tqdm(enumerate(combos), total=len(amtsAtp)*len(ranges)):
        network.substrates["ATP"].max_value = combo[0]
        network.substrates["ATP"].time_ranges = [combo[1]]
        # RS - generalize 
        print(amtsLps[i])
        network.substrates["LPS"].max_value = amtsLps[i]
        network.substrates["LPS"].time_ranges = [combo[1]]

    # check post-training training error
    # args.number
        n=5; # hopefully this is correct
        mean_y, f = network.graph_distributions(time,n , substrates_to_plot=["PI3K", "pAKT", "pPTEN", "Phagocytosis","ATP"], normalize=True, path="./", output_figure=True)
        temp_df = pd.DataFrame(mean_y,columns=list(network.substrates.keys()))
        # args.output
        output="./"
        fname = f"atp_{combo[0]}_{combo[1]}_lps_{amtsLps[i]}.csv"
        if tag is not None:
          fname = tag + "_" + fname
        temp_df.to_csv(os.path.join(output, fname))
        fnames.append(fname)
    return fnames 

def Plotter(df,substrate='ATP',tag="",**kwargs):

  pAkt = np.asarray(df['pAKT'])
  npAkt = pAkt - np.min(pAkt)
  npAkt /= np.max(npAkt)

  ts = np.arange( np.shape(pAkt) [0])

  #plt.plot(ts,npAkt,'r-',label='pAkt (normalized)') # normalized
  try:
    linestyle = kwargs['linestyle']
  except:
    linestyle = 'r-'
  try:
    substratePos = kwargs['substratePos']
    substrateCol = 'b-'
  except:
    substratePos = 0.0175
    substrateCol = 'k-'
  plt.plot(ts,pAkt,linestyle,label='pAkt '+tag)
  plt.ylabel("Conc [AU]")
  plt.xlabel('time [min]')

  #ax2 = plt.twinx()
  #ax2.plot(ts,df['ATP'],label='ATP')
  #ax2.set_ylabel("Conc [AU]")
  daMax = np.max( df[substrate] )
  idx = np.where(np.asarray( df[substrate] )>0.5*daMax)
  idx = idx[0]
  conc = np.ones(np.shape(idx)[0])+substratePos
  plt.plot(idx,conc,substrateCol,linewidth=5,label=substrate )
  plt.legend( loc=0 )


