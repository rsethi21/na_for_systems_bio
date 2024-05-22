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
    ranges = [[60,120],[60,120]],
    amtsLps = [0.0,1.0],    # need to be same length as amtsAtp
    path="pten/input/input_atp_data/",
    params = None,
    tag = None,
    parametersFile = None,
    norm=True   
):
    """
    This function runs the PTEN/PI3K simulation engine. 
    - amtsATP: ATP concentration(s) in AU. Accepts a list (minimum len of 1)
    - ranges: Time range (min) over which ATP concentration is applied. Accepts a list of len 2 lists 
    - amtsLps: LPS concentrations (in same format as amtsATP)
    - params: dictionary of parameters and their values (optional) 
    - parametersFile: json file of pretraining parameters 
    - norm: normalize results based on their value at t=0
    """ 

    # parse and create all necessary objects for creating a network
    #path="na_for_systems_bio"
    rates = parse_rates(path+"/rates.csv")
    substrates = parse_substrates(path+"/substrates.csv", rates)
    interactions = parse_interactions(path+"/interactions.csv", rates, substrates)
    # also need to load fitted_params.json

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

    # load in pretrained parameters if any
    if parametersFile != None:
        with open(parametersFile, "r") as fitted_params:
            parameters = json.load(fitted_params)
        network.set_parameters(list(parameters.values()), list(parameters.keys()))

    

    if params is not None:
      print("Applying parameters", params.keys())
      network.set_parameters(list(params.values()), list(params.keys()))
    
    combos = list(P(amtsAtp, ranges))
    fnames=[]
    for i, combo in tqdm(enumerate(combos), total=len(amtsAtp)*len(ranges)):
        network.substrates["ATP"].max_value = combo[0]
        network.substrates["ATP"].time_ranges = [combo[1]]
    for i, combo in tqdm(enumerate(combos), total=len(amtsLps)*len(ranges)):
        network.substrates["LPS"].max_value =  combo[0]
        network.substrates["LPS"].time_ranges = [combo[1]]

    # check post-training training error
    # args.number
        n=5; # hopefully this is correct
        mean_y, f = network.graph_distributions(
          time,n , 
          substrates_to_plot=["PI3K", "pAKT", "pPTEN", "Phagocytosis","ATP"], 
          normalize=norm, path="./", output_figure=True)
        temp_df = pd.DataFrame(mean_y,columns=list(network.substrates.keys()))
        # args.output
        output="./"
        a=network.substrates["ATP"].max_value
        t=network.substrates["ATP"].time_ranges
        l=network.substrates["LPS"].max_value
        fname = f"atp_{a}_{t}_lps_{l}.csv"
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


#!/usr/bin/env python
import sys
##################################
#
# Revisions
#       10.08.10 inception
#
##################################

#
# ROUTINE  
#
def doit(fileIn):
  1


#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -validation" % (scriptName)
  msg+="""
  
 
Notes:

"""
  return msg

#
# MAIN routine executed when launching this script from command line 
#
if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-validation"):
      #arg1=sys.argv[i+1] 
      #doit(arg1)
      path="input/input_atp_data/"
      runSim(parametersFile="fitted_params/fitted_params.json")
  





  raise RuntimeError("Arguments not understood")




