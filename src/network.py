from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from geneticalgorithm2 import geneticalgorithm2 as ga
import numpy as np
from tqdm import tqdm
import pdb
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

        np.random.seed(2824)
        self.colors = list(mcolors.CSS4_COLORS.keys())
        np.random.shuffle(self.colors)
    
    def get_all_parameters(self):
        parameters = {}
        for interaction in list(self.interactions.values()):
            parameters.update(interaction.get_interaction_parameters())
        return parameters
    
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

    def set_currents(self, new_currents):
        for substrate, new_current in zip(list(self.substrates.values()), new_currents):
            substrate.__setattr__("current_value", new_current)

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

    def reset(self):
        for sub_id, substrate in self.substrates.items():
            substrate.__setattr__("current_value", substrate.initial_value)
        for rate_id, rate in self.parameters.items():
            rate.update_rate(self.initial_parameter_values[rate_id])
    
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

    def get_dydt(self, y, time):
        self.set_currents(y)
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
                    upreg_term = 1*k.__getattribute__("value")*other.__getattribute__("current_value")
                else:
                    upreg_term = 1*k.__getattribute__("value")
                downreg_term = -1*r.__getattribute__("value")*substrate_of_interest.__getattribute__("current_value")
                terms = {k.__getattribute__("identifier"): upreg_term, r.__getattribute__("identifier"): downreg_term}
                for parameter_id, associated_interactions in self.parsed_interactions[substrate_id].items():
                    if parameter_id not in list(terms.keys()):
                        temporary_term = self.parameters[parameter_id].__getattribute__("value")*associated_interactions[0][1] # check this
                    else:
                        temporary_term = terms[parameter_id]
                    for interaction in associated_interactions:
                        if len(interaction) == 2:
                            temporary_term *= self.substrates[interaction[0]].__getattribute__("current_value")
                        else:
                            if interaction[1] == -1:
                                substrate = self.substrates[interaction[0]].__getattribute__("current_value")
                                Km = self.parameters[interaction[2]].__getattribute__("value")
                                n = self.parameters[interaction[3]].__getattribute__("value")
                                temporary_term = temporary_term * (Km**n/(substrate**n + Km**n))
                            else:
                                substrate = self.substrates[interaction[0]].__getattribute__("current_value")
                                Km = self.parameters[interaction[2]].__getattribute__("value")
                                n = self.parameters[interaction[3]].__getattribute__("value")
                                temporary_term = temporary_term * ((substrate**n + Km**n)/Km**n)
                    terms[parameter_id] = temporary_term
                dydt[substrate_id] = sum(list(terms.values()))
            else:
                max_val = substrate_of_interest.__getattribute__("max_value")
                current_val = substrate_of_interest.__getattribute__("current_value")
                time_ranges = substrate_of_interest.__getattribute__("time_ranges")
                if time_ranges == None:
                    between = False
                    after = False
                else:
                    for time_range in reversed(time_ranges):
                        if time >= time_range[0] and time <= time_range[1]:
                            between = True
                            after = False
                            break
                        elif time > time_range[1] and time > time_range[0]:
                            after = True
                            between = False
                            break
                        else:
                            between = False
                            after = False
                if substrate_of_interest.__getattribute__("r") != None and time_ranges != None:
                    r = substrate_of_interest.__getattribute__("r").__getattribute__("value")
                    set_it = [time - tr[0] <= 10 for tr in time_ranges]
                    if current_val != max_val and True in set_it and between:
                        dydt[substrate_id] = max_val - current_val
                    elif after or between:
                        dydt[substrate_id] = -1*r*current_val
                    else:
                        dydt[substrate_id] = 0
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
                pass

        return dydt
    
    def graph(self, time, normalize=False, substrates_to_plot=None, path="./figure.png", output_figure=False):
        if normalize:
            normalize_time = np.linspace(0, 500, 501)
            stimuli_ranges = []
            for s in self.substrates.values():
                stimuli_ranges.append(s.time_ranges)
                s.__setattr__("time_ranges", None)
            probe = odeint(self.get_dydt, self.get_initials(), normalize_time)[-1]
            for index, s in enumerate(self.substrates.values()):
                s.__setattr__("time_ranges", stimuli_ranges[index])

            folds_y = odeint(self.get_dydt, probe, time)
            y = folds_y.copy()
            for t in range(y.shape[0]):
                for index, s in enumerate(self.substrates.values()):
                    if s.type == "non-stimulus":
                        y[t,index] = folds_y[t,index]/probe[index]
                    else:
                        y[t, index] = folds_y[t,index]
        else:
            y = odeint(self.get_dydt, self.get_initials(), time)
        fig = plt.figure() 
        if path != None:
            for i, substrate in tqdm(enumerate(list(self.substrates.values())), desc="Plotting Each Substrate", total=len(self.substrates)):
                if substrate.identifier in substrates_to_plot:
                    plt.plot(time, y[:,i], self.colors[i], label=substrate.__getattribute__("identifier"))
            plt.xlabel("Time (mins)",fontsize=12)
            plt.ylabel("Concentration (AU)",fontsize=12)
            plt.legend(loc="upper right", fontsize=5)
            fig.savefig(path)
        plt.close(fig)
        if output_figure:
            return y, fig
        else:
            return y

    def graph_distributions(self, time, samples, normalize=False, substrates_to_plot=[], path="./figure_with_area.png", output_figure=False, initials=None, verbose=True):
        if initials == None:
            y0s = []
            for _ in tqdm(range(samples),desc="Generating Random Initial",total=samples, disable=~verbose):
                y0 = []
                for i, s in enumerate(self.substrates.values()):
                    if s.type == "stimulus":
                        y0.append(0.0)
                    else:
                        y0.append(2**np.random.randn())
                y0s.append(y0)
        else:
            y0s = initials
        ys = []
        for y0 in tqdm(y0s, desc="Generating Simulations", disable=~verbose):
            self.set_initials(y0)
            ys.append(self.graph(time, normalize=normalize, path=None))
        min_y = np.mean(ys, axis=0) - np.std(ys, axis=0)*2
        max_y = np.mean(ys, axis=0) + np.std(ys, axis=0)*2
        mean_y = np.mean(ys, axis=0)
        temp_fig = plt.figure()
        if path != None:
            for i, s in tqdm(enumerate(self.substrates.values()), desc="Plotting Each Substrates", total=len(self.substrates), disable=~verbose):
                if s.identifier in substrates_to_plot:
                    if s.type != "stimulus":
                        plt.fill_between(time, min_y[:,i], max_y[:,i], color=self.colors[i], alpha=0.2)
                    plt.plot(time, mean_y[:,i], color=self.colors[i],label=s.__getattribute__("identifier"),linewidth=1.0)
            plt.xlabel("Time (mins)")
            plt.ylabel("Concentration (AU)")
            plt.legend(loc="upper right", fontsize=5)
            temp_fig.savefig(path)
            print(f"Figure saved: {path}")
        plt.close(temp_fig)
        if output_figure:
            return mean_y, temp_fig
        else:
            return mean_y

    def fit(self, data, time, arguments, number=1, normalize=False, obj_calc=None, mlp=1):
        bounds, bound_types, names = self.get_bounds_information()
        substrate_names = list(self.substrates.keys())
        y0s = []
        for _ in tqdm(range(number),desc="Generating Random Initial",total=number, disable=True):
            y0 = []
            for i, s in enumerate(self.substrates.values()):
                if s.type == "stimulus":
                    y0.append(0.0)
                else:
                    y0.append(2**np.random.randn())
            y0s.append(y0)
        def loss(X):
            self.set_parameters(X, names)
            cost = 0
            predictions = self.graph_distributions(time, number, initials=y0s, normalize=normalize, path=None, verbose=False)
            for substrate_id, substrate_data in data.items():
                for time_point, truth in substrate_data.items():
                    prediction = predictions[int(time_point), substrate_names.index(substrate_id)]
                    if obj_calc == None:
                        cost += (prediction - float(truth))**2
                    else:
                        cost += obj_calc(prediction, truth)
            return cost

        fitting_model = ga(function = loss,
                           dimension = len(bounds),
                           variable_type = bound_types,
                           variable_boundaries = bounds,
                           algorithm_parameters = arguments
                           )
        
        fitting_model.run(set_function=ga.set_function_multiprocess(loss, n_jobs=mlp))
        self.set_parameters(fitting_model.result.variable, names)
        return names, fitting_model.result.variable

