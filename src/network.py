from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from geneticalgorithm2 import geneticalgorithm2 as ga
import numpy as np
from tqdm import tqdm
import pdb
import os
import joblib

class Network:
    def __init__(self, identifier: str, rates: list, interactions: list, substrates: list):
        self.identifier = identifier
        self.interactions = {interaction.identifier: interaction for interaction in interactions}
        self.substrates = {substrate.identifier: substrate for substrate in substrates}
        self.parameters = {rate.identifier: rate for rate in rates}
        self.initial_parameter_values = {parameter_id: parameter.__getattribute__("value") for parameter_id, parameter in self.parameters.items()}
        self.parsed_interactions = self.parse_interactions()
        self.rates = rates
        self.colors = list(mcolors.CSS4_COLORS.keys())
        self.colors = sorted(self.colors)

    def get_all_parameters(self):
        parameters = {}
        for interaction in list(self.interactions.values()):
            parameters.update(interaction.get_interaction_parameters())
        return parameters
    
    def find_rate_pairs(self):
        rate_pairs = {s.k.identifier+","+s.r.identifier: [s.k, s.r] for s in list(self.substrates.values()) if s.type == "non-stimulus"}
        parameters_accounted = []
        rate_pairs_filtered = {}
        for pair in list(rate_pairs.keys()):
            p1, p2 = pair.split(",")
            if p2 in parameters_accounted or p1 in parameters_accounted:
                pass
            else:
                parameters_accounted.extend([p1, p2])
                rate_pairs_filtered[pair] = rate_pairs[pair]
        return rate_pairs_filtered

    def apply_scaling(self, scales, names):
        pairs = self.find_rate_pairs()
        for i in names:
            if i in list(self.parameters.keys()):
                parameter = self.parameters[i]
                parameter.update_rate(parameter.value * scales[names.index(i)])
            else:
                v1, v2 = pairs[i]
                p1, p2 = i.split(",")
                v1.update_rate(v1.value * scales[names.index(i)])
                v2.update_rate(v2.value * scales[names.index(i)])

    def get_bounds_information(self):
        bounds = []
        bounds_types = []
        rate_names = []
        for parameter in list(self.parameters.values()):
            if not parameter.fixed:
                bounds.append(parameter.__getattribute__("bounds"))
                bounds_types.append(parameter.__getattribute__("bounds_type"))
                rate_names.append(parameter.__getattribute__("identifier"))
        return bounds, bounds_types, rate_names

    def set_parameters(self, X, rate_names):
        for i, parameter in self.parameters.items():
            if i in rate_names:
                parameter.update_rate(X[rate_names.index(i)])

    def get_initials(self):
        y0 = []
        for substrate in list(self.substrates.values()):
            y0.append(substrate.__getattribute__("initial_value"))
        return y0
    
    def set_initials(self, new_initials):
        for substrate, new_initial in zip(list(self.substrates.values()), new_initials):
            substrate.__setattr__("initial_value", new_initial)

    def set_currents(self, new_initials):
        for substrate, new_initial in zip(list(self.substrates.values()), new_initials):
            substrate.__setattr__("current_value", new_initial)

    def reset_intials(self):
        for sub_id, substrate in self.substrates.items():
            substrate.__setattr__("current_value", substrate.initial_value)
        for rate_id, rate in self.parameters.items():
            rate.update_rate(self.initial_parameter_values[rate_id])
    
    def reset_stimuli(self):
        for sub_id, substrate in self.substrates.items():
            if substrate.__getattribute__("type") == "stimulus":
                substrate.__setattr__("max_value", 0.0)
                substrate.__setattr__("time_ranges", None)

    def unfreeze_paramters(self, parameter_ids):
        for parameter_id in parameter_ids:
            self.parameters[parameter_id].__setattr__("fixed", False)

    def freeze_parameters(self, parameter_ids):
        for parameter_id in parameter_ids:
            self.parameters[parameter_id].__setattr__("fixed", True)

    def parse_interactions(self):
        parsed_interactions = {name: {} for name in self.substrates.keys()}
        for id, interaction in self.interactions.items():
            substrate_id = interaction.__getattribute__("s1").__getattribute__("identifier")
            rate_id = interaction.__getattribute__("rate").__getattribute__("identifier")
            substrate2 = interaction.__getattribute__("s2").__getattribute__("identifier")
            Km = interaction.__getattribute__("Km")
            n = interaction.__getattribute__("n")
            if  Km != None and n != None:
                if rate_id not in parsed_interactions[substrate_id].keys():
                    parsed_interactions[substrate_id][rate_id] = [(substrate2, interaction.__getattribute__("effect"), Km.__getattribute__("identifier"), n.__getattribute__("identifier"))]
                else:
                    parsed_interactions[substrate_id][rate_id].append((substrate2, interaction.__getattribute__("effect"), Km.__getattribute__("identifier"), n.__getattribute__("identifier")))
            else:
                if rate_id not in parsed_interactions[substrate_id].keys():
                    parsed_interactions[substrate_id][rate_id] = [(substrate2, interaction.__getattribute__("effect"))]
                else:
                    parsed_interactions[substrate_id][rate_id].append((substrate2, interaction.__getattribute__("effect")))
        return parsed_interactions

    # piece together derivatives for plotting purposes
    def get_dydt(self, y, time, parameter_sets=None):
        '''
        Inputs: current state of substrates, current time
        Outputs: predicted change in substrates
        '''
        self.set_currents(y) # set the current state of system substrates
        if parameter_sets != None:
            for pi, config in parameter_sets.items():
                if time >= config["times"][0] and time < config["times"][1]:
                    self.parameters[pi].value = config["value"]
                else:
                    self.parameters[pi].value = self.initial_parameter_values[pi]
        dydt = {substrate_id: 0 for substrate_id in list(self.substrates.keys())} # instantiate dictionary to hold differential equations
        # iterate through each substrate
        for substrate_id in list(dydt.keys()):
            substrate_of_interest = self.substrates[substrate_id] # extract substrate object
            # check if substrate is not a stimulus
            if substrate_of_interest.__getattribute__("type") == "non-stimulus":
                k = substrate_of_interest.__getattribute__("k") # extract phophorylation/activation rate
                r = substrate_of_interest.__getattribute__("r") # extract dephosphorylation/deactivation rate
                if substrate_of_interest.other != None:
                    upreg_term = 1*k.__getattribute__("value")*self.substrates[substrate_of_interest.other].__getattribute__("current_value") # by default apply opposite form of the substrate to positive
                else:
                    upreg_term = 1*k.__getattribute__("value") # otherwise just instantiate with activation rate
                downreg_term = -1*r.__getattribute__("value")*substrate_of_interest.__getattribute__("current_value") #  by default substrate will be deacticated by itself
                terms = {k.__getattribute__("identifier"): upreg_term, r.__getattribute__("identifier"): downreg_term} # store terms for ease of parsing interactions
                for parameter_id, associated_interactions in self.parsed_interactions[substrate_id].items(): # iterate throat each interaction parameter related to a specific substrate
                    if parameter_id not in list(terms.keys()):
                        temporary_term = self.parameters[parameter_id].__getattribute__("value")*associated_interactions[0][1] # add any additional terms that are not related to activation/deactivation rates
                    else:
                        temporary_term = terms[parameter_id] # if already in terms then extract the parameter of interest
                    for interaction in associated_interactions: # iterate throat each interaction pair
                        if len(interaction) == 2: # check if first order iteraction between substrate
                            temporary_term *= self.substrates[interaction[0]].__getattribute__("current_value") # apply substrate iteraction to current term
                        else: #  else assume using hill equation
                            if interaction[1] == -1: # if negative iteraction then negative feedback term
                                substrate = self.substrates[interaction[0]].__getattribute__("current_value")
                                Km = self.parameters[interaction[2]].__getattribute__("value")
                                n = self.parameters[interaction[3]].__getattribute__("value")
                                temporary_term = temporary_term * (Km**n/(substrate**n + Km**n))
                            else: # else then positive feedback term
                                substrate = self.substrates[interaction[0]].__getattribute__("current_value")
                                Km = self.parameters[interaction[2]].__getattribute__("value")
                                n = self.parameters[interaction[3]].__getattribute__("value")
                                temporary_term = temporary_term * ((substrate**n + Km**n)/Km**n)
                    terms[parameter_id] = temporary_term # apply calculations of the term of interest to the appropriate location in terms dictionary
                dydt[substrate_id] = sum(list(terms.values()))
            # if substrate is a stimulus
            else:
                max_val = substrate_of_interest.__getattribute__("max_value") # extract max value
                current_val = substrate_of_interest.__getattribute__("current_value") # extract current_value/initial value
                time_ranges = substrate_of_interest.__getattribute__("time_ranges") # extract time range of interaction
                if time_ranges == None: # check if any time ranges specified for application of stimulus
                    between = False
                    after = False
                else: # manually check at each of the time ranges
                    for i, time_range in enumerate(time_ranges): # iterate through each time range
                        # check if current time is in between the current range probed
                        if time >= time_range[0] and time <= time_range[1]:
                            between = True
                            after = False
                            break
                        # check if current time is after the current range probed
                        elif time > time_range[1] and time > time_range[0]:
                            # if this is the last range to probe, then the after rate can be applied
                            if i == len(time_ranges) - 1:
                                after = True
                                between = False
                                break
                            # if this is not the last range to probe and not in next range to probe then after rate can be applied
                            elif time_ranges[i+1][0] > time:
                                after = True
                                between = False
                                break
                            # if this is not the last range to probe and in the next range then allow for check in the next range
                            elif time_ranges[i+1][0] <= time:
                                continue
                        # if not between or after, then before so no rate applied
                        else:
                            between = False
                            after = False
                # check to see if the stimuli has a deprication rate assigned
                if substrate_of_interest.__getattribute__("r") != None and time_ranges != None:
                    # extract deprication rate
                    r = substrate_of_interest.__getattribute__("r").__getattribute__("value")
                    # check to see if the substrate has reached the max value within certain error (1000ths by default)
                    if abs(current_val - max_val) <= 0.001:
                        substrate_of_interest.reached = True # marker for if reached
                    # if between and not reached then allow to reach max value
                    if between and not substrate_of_interest.reached:
                        dydt[substrate_id] = max_val - current_val
                    # if between and reached, then apply deprication rate
                    elif between and substrate_of_interest.reached:
                        dydt[substrate_id] = -1*r*current_val
                    # if after then apply depcrication
                    elif after:
                        dydt[substrate_id] = -1*r*current_val
                    # if before then no rate
                    else:
                        dydt[substrate_id] = 0
                # apply same conditions but without deprication
                else:
                    if between:
                        dydt[substrate_id] = max_val - current_val
                    elif after:
                        dydt[substrate_id] = -current_val
                    else:
                        dydt[substrate_id] = 0
        return list(dydt.values())

    def get_representation_dydt(self):
        dydt = {substrate_id: 0 for substrate_id in list(self.substrates.keys())}
        for substrate_id in list(dydt.keys()):
            substrate_of_interest = self.substrates[substrate_id]
            if substrate_of_interest.other != None:
                other = self.substrates[substrate_of_interest.other]
            else:
                other = None
            if substrate_of_interest.__getattribute__("type") == "non-stimulus":
                k = substrate_of_interest.__getattribute__("k")
                r = substrate_of_interest.__getattribute__("r")
                if substrate_of_interest.other != None:
                    upreg_term = f"{k.__getattribute__('identifier')}[{other.__getattribute__('identifier')}]"
                else:
                    upreg_term = f"{k.__getattribute__('identifier')}"
                downreg_term = f" - {r.__getattribute__('identifier')}[{substrate_of_interest.__getattribute__('identifier')}]"
                terms = {k.__getattribute__("identifier"): upreg_term, r.__getattribute__("identifier"): downreg_term}
                for parameter_id, associated_interactions in self.parsed_interactions[substrate_id].items():
                    if parameter_id not in list(terms.keys()):
                        if associated_interactions[0][1] == -1:
                            temporary_term = f" - {parameter_id}"
                        else:
                            temporary_term = f" + {parameter_id}"
                    else:
                        temporary_term = terms[parameter_id]
                    for interaction in associated_interactions:
                        if len(interaction) == 2:
                            temporary_term += f"[{self.substrates[interaction[0]].__getattribute__('identifier')}]"
                        else:
                            if interaction[1] == -1:
                                substrate = self.substrates[interaction[0]].__getattribute__("identifier")
                                Km = self.parameters[interaction[2]].__getattribute__("identifier")
                                n = self.parameters[interaction[3]].__getattribute__("identifier")
                                temporary_term = temporary_term + f"[{Km}^{n}/({substrate}^{n} + {Km}^{n})]"
                            else:
                                substrate = self.substrates[interaction[0]].__getattribute__("identifier")
                                Km = self.parameters[interaction[2]].__getattribute__("identifier")
                                n = self.parameters[interaction[3]].__getattribute__("identifier")
                                temporary_term = temporary_term + f"[({substrate}^{n} + {Km}^{n})/{Km}^{n}]"
                    terms[parameter_id] = temporary_term
                overall = ""
                key = f"d[{substrate_id}]/dt"
                for term in list(terms.values()):
                    overall += term
                dydt[key] = overall
                del dydt[substrate_id]
            else:
                key = f"d[{substrate_id}]/dt"
                del dydt[substrate_id]
                if substrate_of_interest.r != None:
                    dydt[key] = f"-{substrate_of_interest.r.identifier}[{substrate_id}]"
                else:
                    dydt[key] = 0
        return dydt
    
    def graph(self, time, max_normalize=False, normalize=False, folds_normalize=False, steady_normalize=False, substrates_to_plot=[], path="./figure.png", output_figure=False, parameter_sets=None):
        assert sum([max_normalize, normalize, steady_normalize, folds_normalize]) == 1 or sum([max_normalize, normalize, steady_normalize, folds_normalize]) == 0, "Can only use one normalization strategy at a time (steady state, max for each substrate, steady state folds, folds, or none)"
        if steady_normalize:
            normalize_time = np.linspace(0, time[-1], time[-1]+1)
            stimuli_ranges = []
            stimuli_amts = []
            for s in self.substrates.values():
                stimuli_ranges.append(s.time_ranges)
                stimuli_amts.append(s.max_value)
                s.__setattr__("time_ranges", None)
                s.__setattr__("max_value", 0)
            probe = odeint(self.get_dydt, self.get_initials(), normalize_time)[-1]
            for index, s in enumerate(self.substrates.values()):
                s.__setattr__("time_ranges", stimuli_ranges[index])
                s.__setattr__("max_value", stimuli_amts[index])
        elif normalize:
            normalize_time = np.linspace(0, time[-1], time[-1]+1)
            stimuli_ranges = []
            stimuli_amts = []
            for s in self.substrates.values():
                stimuli_ranges.append(s.time_ranges)
                stimuli_amts.append(s.max_value)
                s.__setattr__("time_ranges", None)
                s.__setattr__("max_value", 0)
            probe = odeint(self.get_dydt, self.get_initials(), normalize_time)[-1]
            for index, s in enumerate(self.substrates.values()):
                s.__setattr__("time_ranges", stimuli_ranges[index])
                s.__setattr__("max_value", stimuli_amts[index])
            folds_y = odeint(self.get_dydt, probe, time, args=(parameter_sets,))
            y = folds_y.copy()
            for t in range(y.shape[0]):
                for index, s in enumerate(self.substrates.values()):
                    if s.__getattribute__("type") == "non-stimulus" and probe[index] > 0:
                        y[t,index] = folds_y[t,index]/probe[index]
                    else:
                        y[t, index] = folds_y[t,index]
        elif folds_normalize:
            probe = self.get_initials()
            folds_y = odeint(self.get_dydt, probe, time, args=(parameter_sets,))
            y = folds_y.copy()
            for t in range(y.shape[0]):
                for index, s in enumerate(self.substrates.values()):
                    if s.__getattribute__("type") == "non-stimulus" and probe[index] > 0:
                        y[t,index] = folds_y[t,index]/probe[index]
                    else:
                        y[t, index] = folds_y[t,index]
        elif max_normalize:
            folds_y = odeint(self.get_dydt, probe, time, args=(parameter_sets,))
            probe = np.max(folds_y, axis=0)
            y = folds_y.copy()
            for t in range(y.shape[0]):
                for index, s in enumerate(self.substrates.values()):
                    if s.__getattribute__("type") == "non-stimulus" and probe[index] > 0:
                        y[t,index] = folds_y[t,index]/probe[index]
                    else:
                        y[t, index] = folds_y[t,index]
        else:
            y = odeint(self.get_dydt, self.get_initials(), time, args=(parameter_sets,))
        for s in self.substrates.values():
            s.__setattr__("reached", False)
        fig = plt.figure()
        if path != None:
            if len(substrates_to_plot) == 0:
                substrates_to_plot = list(self.substrates.keys())
            for i, substrate in tqdm(enumerate(list(self.substrates.values())), desc="Plotting Each Substrate", total=len(self.substrates)):
                if substrate.identifier in substrates_to_plot:
                    plt.plot(time, y[:,i], self.colors[i], label=substrate.__getattribute__("identifier"))
            plt.xlabel("Time (mins)",fontsize=12)
            plt.ylabel("Concentration (AU)",fontsize=12)
            plt.legend(loc="upper right", fontsize=5)
            fig.savefig(path)
        if output_figure:
            return y, fig
        else:
            plt.close(fig)
            return y

    def graph_distributions(self, time, samples, max_normalize=False, normalize=False, substrates_to_plot=[], path="./figure_with_area.png", output_figure=False, initials=None, verbose=True, axis=None, return_min_max=False, parameter_sets=None):
        if initials == None:
            y0s = []
            for _ in tqdm(range(samples),desc="Generating Random Initial",total=samples, disable=~verbose):
                y0 = []
                for i, s in enumerate(self.substrates.values()):
                    if s.__getattribute__("type") == "stimulus":
                        y0.append(0.0)
                    else:
                        y0.append(2**np.random.randn())
                y0s.append(y0)
        else:
            y0s = initials
        ys = []
        for y0 in tqdm(y0s, desc="Generating Simulations", disable=~verbose):
            self.set_initials(y0)
            ys.append(self.graph(time, max_normalize=max_normalize, normalize=normalize, path=None, parameter_sets=parameter_sets))
        min_y = np.mean(ys, axis=0) - np.std(ys, axis=0)*2
        max_y = np.mean(ys, axis=0) + np.std(ys, axis=0)*2
        mean_y = np.mean(ys, axis=0)
        if path != None or axis != None or output_figure:
            if len(substrates_to_plot) == 0:
                substrates_to_plot = list(self.substrates.keys())
            if axis == None:
                temp_fig = plt.figure()
                for i, s in tqdm(enumerate(self.substrates.values()), desc="Plotting Each Substrates", total=len(self.substrates), disable=~verbose):
                    if s.identifier in substrates_to_plot:
                        if s.__getattribute__("type") != "stimulus":
                            plt.fill_between(time, min_y[:,i], max_y[:,i], color=self.colors[i], alpha=0.2)
                        plt.plot(time, mean_y[:,i], color=self.colors[i],label=s.__getattribute__("identifier"),linewidth=1.0)
                plt.xlabel("Time (mins)")
                plt.ylabel("Concentration (AU)")
                plt.legend(loc="upper right", fontsize=5)
                if path != None:
                    temp_fig.savefig(path)
                    plt.close(temp_fig)
            else:
                for i, s in tqdm(enumerate(self.substrates.values()), desc="Plotting Each Substrates", total=len(self.substrates), disable=~verbose):
                    if s.identifier in substrates_to_plot:
                        if s.__getattribute__("type") != "stimulus":
                            axis.fill_between(time, min_y[:,i], max_y[:,i], color=self.colors[i], alpha=0.2)
                        axis.plot(time, mean_y[:,i], color=self.colors[i],label=s.__getattribute__("identifier"),linewidth=1.0)
        if output_figure and return_min_max and axis==None:
            return mean_y, min_y, max_y, temp_fig
        elif output_figure and axis==None:
            return mean_y, temp_fig
        elif return_min_max:
            return mean_y, min_y, max_y
        else:
            return mean_y

    def scale(self, data, time, arguments, scaling_bounds=[1,5], initials=None, number=1, normalize=False, obj_calc=None, mlp=1, plots_path=None):
        '''
        Inputs:
        - Required: data, time frame to fit against, fitting algorithm arguments
        - Optional: number of random conditions, set of initials if certain initial conditions to test, whether to apply fold normalization, objective function if different, mlp processing cores to use
        Outputs:
        - names of parameters, fitted parameters
        Training Data Format: json
        Example:
            [
                {
                    stimuli: ["",],
                    max_values:[_,],
                    time_ranges: [ [[_,_]], ],
                    substrates: 
                    {
                        "": 
                            {
                                time: _,
                            },
                    }
                },
            ]
        '''
        # save original stimuli configurations before running fitting on multiple conditions
        original_stimuli_configs = {}
        for substrate in self.substrates.values():
            if substrate.__getattribute__("type") == "stimulus":
                original_stimuli_configs[substrate.__getattribute__("identifier")] = [substrate.__getattribute__("max_value"), substrate.__getattribute__("time_ranges")]
        bounds, bound_types, names = self.get_bounds_information() # extract information needed for fitting
        pair_names = list(self.find_rate_pairs().keys())
        pairs_unraveled = []
        for pair in pair_names:
            pairs_unraveled.extend(pair.split(","))
        names = [n for n in names if n not in pairs_unraveled]
        names.extend(pair_names)
        bounds = [scaling_bounds for _ in range(len(bounds))]
        bound_types = ["int" for _ in range(len(bound_types))]
        substrate_names = list(self.substrates.keys()) # extract substrate order information for indexing purposes
        # randomly generate starting conditions to generate a more robust fit for long-term system behavior (don't need to apply this if want to fit to a given set of intital conditions)
        if initials == None:
            y0s = [] # instantiate empty list
            for _ in tqdm(range(number),desc="Generating Random Initial",total=number, disable=True): # iterate through the number of randomly generated intial conditions
                y0 = []
                for i, s in enumerate(self.substrates.values()): # iterate through all the substrates
                    if s.__getattribute__("type") == "stimulus": # if the substrate is a stimulus category, use 0 for initial condition
                        y0.append(0.0)
                    else:
                        y0.append(2**np.random.randn()) # else randomly generate one
                y0s.append(y0) # append sample conditions to a list of multiple random conditions
        else:
            y0s = initials
        # calculate loss between model and data to fit against
        def loss(X):
            self.apply_scaling([10**(-1*x) for x in X], names) # setting parameters from current solution space
            cost = 0 # instantiating cost
            count = 0 # instantiating record for calculating average loss
            # iterate through fitting data for each set of conditions data is collected
            for entry in data:
                # apply conditions from current fitting data condition set
                for i, stimulus in enumerate(entry["stimuli"]):
                    self.substrates[stimulus].max_value = entry["max_values"][i]
                    self.substrates[stimulus].time_ranges = entry["time_ranges"][i]
                # generating model predictions with conditions and current paramters solution set
                predictions = self.graph_distributions(time, number, initials=y0s, normalize=normalize, path=None, verbose=False)
                # iterate through each substrate there is data from in that condition
                for substrate_id, substrate_data in entry["substrates"].items():
                    # iterate throat each time point there is data for the substrate in that specific condition
                    for time_point, truth in substrate_data.items():
                        # extract the appropriate prediction
                        prediction = predictions[int(time_point), substrate_names.index(substrate_id)]
                        # calculate cost (either SE or any that one chooses)
                        if obj_calc == None:
                            cost += (prediction - float(truth))**2
                            count += 1
                        else:
                            cost += obj_calc(prediction, truth)
                            count += 1
                self.reset_stimuli()
            return cost/count
        # apply loss function and appropriate parameters to the genetic algorithm implementation
        fitting_model = ga(function = loss,
                           dimension = len(bounds),
                           variable_type = bound_types,
                           variable_boundaries = bounds,
                           algorithm_parameters = arguments
                           )
        
        # run fitting
        fitting_model.run(set_function=ga.set_function_multiprocess(loss, n_jobs=mlp))
        # apply fitted conditions
        self.apply_scaling([10**(-1*s) for s in fitting_model.result.variable], names)
        # plot fit if requested
        if plots_path != None:
            for j, entry in enumerate(data):
                for i, stimulus in enumerate(entry["stimuli"]):
                    self.substrates[stimulus].max_value = entry["max_values"][i]
                    self.substrates[stimulus].time_ranges = entry["time_ranges"][i]
                subs = list(entry["substrates"].keys())
                subs.extend(entry["stimuli"])
                y, fig = self.graph_distributions(time, number, normalize=normalize, output_figure=True, substrates_to_plot=subs)
                plt.figure(fig)
                for substrate_id, time_value_pairs in entry["substrates"].items():
                    index = list(self.substrates.keys()).index(substrate_id)
                    for t, value in time_value_pairs.items():
                        plt.plot(int(t), value, marker=".", color=self.colors[i])
                fig.savefig(os.path.join(plots_path, f"condition_{j}.png"))
                plt.close(fig)
                self.reset_stimuli()
        # restore original user specified configurations
        for stimuli_id, configs in original_stimuli_configs.items():
            self.substrates[stimuli_id].max_value = configs[0]
            self.substrates[stimuli_id].time_ranges = configs[1]
        return names, fitting_model.result.variable

    # this is a function that allows the user to search for solutions spaces of model parameters to align with one dataset
    def fit(self, data, time, arguments, initials=None, number=1, normalize=False, obj_calc=None, mlp=1, plots_path=None):
        '''
        Inputs:
        - Required: data, time frame to fit against, fitting algorithm arguments
        - Optional: number of random conditions, set of initials if certain initial conditions to test, whether to apply fold normalization, objective function if different, mlp processing cores to use
        Outputs:
        - names of parameters, fitted parameters
        Training Data Format: json
        Example:
            [
                {
                    stimuli: ["",],
                    max_values:[_,],
                    time_ranges: [ [[_,_]], ],
                    substrates: 
                    {
                        "": 
                            {
                                time: _,
                            },
                    }
                },
            ]
        '''
        # save original stimuli configurations before running fitting on multiple conditions
        original_stimuli_configs = {}
        for substrate in self.substrates.values():
            if substrate.__getattribute__("type") == "stimulus":
                original_stimuli_configs[substrate.__getattribute__("identifier")] = [substrate.__getattribute__("max_value"), substrate.__getattribute__("time_ranges")]
        bounds, bound_types, names = self.get_bounds_information() # extract information needed for fitting
        substrate_names = list(self.substrates.keys()) # extract substrate order information for indexing purposes
        # randomly generate starting conditions to generate a more robust fit for long-term system behavior (don't need to apply this if want to fit to a given set of intital conditions)
        if initials == None:
            y0s = [] # instantiate empty list
            for _ in tqdm(range(number),desc="Generating Random Initial",total=number, disable=True): # iterate through the number of randomly generated intial conditions
                y0 = []
                for i, s in enumerate(self.substrates.values()): # iterate through all the substrates
                    if s.__getattribute__("type") == "stimulus": # if the substrate is a stimulus category, use 0 for initial condition
                        y0.append(0.0)
                    else:
                        y0.append(2**np.random.randn()) # else randomly generate one
                y0s.append(y0) # append sample conditions to a list of multiple random conditions
        else:
            y0s = initials
        # calculate loss between model and data to fit against
        def loss(X):
            self.set_parameters(X, names) # setting parameters from current solution space
            cost = 0 # instantiating cost
            count = 0 # instantiating record for calculating average loss
            # iterate through fitting data for each set of conditions data is collected
            for entry in data:
                # apply conditions from current fitting data condition set
                for i, stimulus in enumerate(entry["stimuli"]):
                    self.substrates[stimulus].max_value = entry["max_values"][i]
                    self.substrates[stimulus].time_ranges = entry["time_ranges"][i]
                # generating model predictions with conditions and current paramters solution set
                predictions = self.graph_distributions(time, number, initials=y0s, normalize=normalize, path=None, verbose=False)
                # iterate through each substrate there is data from in that condition
                for substrate_id, substrate_data in entry["substrates"].items():
                    # iterate throat each time point there is data for the substrate in that specific condition
                    for time_point, truth in substrate_data.items():
                        # extract the appropriate prediction
                        prediction = predictions[int(time_point), substrate_names.index(substrate_id)]
                        # calculate cost (either SE or any that one chooses)
                        if obj_calc == None:
                            cost += (prediction - float(truth))**2
                            count += 1
                        else:
                            cost += obj_calc(prediction, truth)
                            count += 1
                self.reset_stimuli()
            return cost/count
        # apply loss function and appropriate parameters to the genetic algorithm implementation
        fitting_model = ga(function = loss,
                           dimension = len(bounds),
                           variable_type = bound_types,
                           variable_boundaries = bounds,
                           algorithm_parameters = arguments
                           )
        
        # run fitting
        fitting_model.run(set_function=ga.set_function_multiprocess(loss, n_jobs=mlp))
        # apply fitted conditions
        self.set_parameters(fitting_model.result.variable, names)
        # plot fit if requested
        if plots_path != None:
            for j, entry in enumerate(data):
                for i, stimulus in enumerate(entry["stimuli"]):
                    self.substrates[stimulus].max_value = entry["max_values"][i]
                    self.substrates[stimulus].time_ranges = entry["time_ranges"][i]
                subs = list(entry["substrates"].keys())
                subs.extend(entry["stimuli"])
                y, fig = self.graph_distributions(time, number, normalize=normalize, output_figure=True, substrates_to_plot=subs)
                plt.figure(fig)
                for substrate_id, time_value_pairs in entry["substrates"].items():
                    index = list(self.substrates.keys()).index(substrate_id)
                    for t, value in time_value_pairs.items():
                        plt.plot(int(t), value, marker=".", color=self.colors[i])
                fig.savefig(os.path.join(plots_path, f"condition_{j}.png"))
                plt.close(fig)
                self.reset_stimuli()
        # restore original user specified configurations
        for stimuli_id, configs in original_stimuli_configs.items():
            self.substrates[stimuli_id].max_value = configs[0]
            self.substrates[stimuli_id].time_ranges = configs[1]
        return names, fitting_model.result.variable
