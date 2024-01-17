from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import geneticalgorithm2 as ga

class Network:
    def __init__(self, identifier: str, rates: list, interactions: list, substrates: list):
        self.identifier = identifier
        self.interactions = {interaction.identifier: interaction for interaction in interactions}
        self.substrates = {substrate.identifier: substrate for substrate in substrates}
        self.parameters = self.get_all_parameters()
        self.initial_parameter_values = {parameter_id: parameter.__getattribute__("value") for parameter_id, parameter in self.parameters.items()}
        self.parsed_interactions = self.parse_interactions()
        self.rates = rates
    
    def get_all_parameters(self):
        parameters = {}
        for interaction in list(self.interactions.values()):
            parameters.update(interaction.get_interaction_parameters())
        return parameters
    
    def get_bounds_information(self):
        bounds = []
        bounds_types = []
        for parameter in list(self.parameters.values()):
            if not parameter.fixed:
                bounds.append(parameter.__getattribute__("bounds"))
                bounds_types.append(parameter.__getattribute__("bounds_type"))
        return bounds, bounds_types

    def set_parameters(self, X):
        for i, parameter in enumerate(list(self.parameters.values())):
            if not parameter.fixed:
                parameter.update_rate(X[i])

    def set_currents(self, new_currents):
        for substrate, new_current in zip(list(self.substrates.values()), new_currents):
            substrate.__setattr__("current_value", new_current)

    def get_initials(self):
        y0 = []
        for substrate in list(self.substrates.values()):
            y0.append(substrate.__getattribute__())
        return y0
    
    def reset(self):
        for sub_id, substrate in self.substrates.items():
            substrate.__setattr__("current_value", substrate.initial_value)
        for rate_id, rate in self.parameters.items():
            rate.update_rate(self.initial_parameter_values[rate_id])
    
    def freeze_parameters(self, parameter_ids):
        for parameter_id in parameter_ids:
            self.parameters[parameter_id].__setattr__("fixed", True)

    def parse_interactions(self):
        parsed_interactions = {substrate.identifier: {} for substrate in self.substrates}
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
            if substrate_of_interest.__getattribute__("type") == "non-stimulus":
                k = substrate_of_interest.__getattribute__("k")
                r = substrate_of_interest.__getattribute__("r")
                upreg_term = 1*k
                downreg_term = -1*r*substrate_of_interest.__getattribute__("current_value")
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
                            substrate = self.substrates[interaction[0]].__getattribute__("current_value")
                            Km = self.parameters[interaction[2]].__getattribute__("value")
                            n = self.parameters[interaction[3]].__getattribute__("value")
                            temporary_term = temporary_term * (Km**n/(substrate**n + Km**n))
                    terms[parameter_id] = temporary_term
                dydt[substrate_id] = sum(list(terms.values()))
            else:
                max_val = substrate_of_interest.__getattribute__("max_value")
                current_val = substrate_of_interest.__getattribute__("current_value")
                time_ranges = substrate_of_interest.__getattribute__("time_ranges")
                between = False
                after = False
                if len(time_ranges) == 0:
                    between = True
                    after = False
                else:
                    for time_range in time_ranges:
                        if time >= time_range[0] and time <= time_range[1]:
                            between = True
                            after = False
                            break
                        elif time > time_range[1]:
                            after = True
                            between = False
                        else:
                            after = False
                            between = False
                if between:
                    dydt[substrate_id] = max_val - current_val
                elif after:
                    dydt[substrate_id] = -current_val

        return list(dydt.values())
    
    def graph(self, time, path="./figure.png"):
        colors = mcolors.CSS4_COLORS
        y = odeint(self.get_dydt, self.get_initials(), time)
        fig = plt.figure()
        for i, substrate in list(self.substrates.values()):
            plt.plot(time, y[:,i], colors[i], label=substrate.__getattribute__("identifier"))
        plt.xlabel("Time (mins)",fontsize=12)
        plt.ylabel("Concentration (AU)",fontsize=12)
        plt.legend(loc="upper right", fontsize=5)
        fig.savefig(path)
        plt.close(fig)
        return y

    def fit(self, data, time, arguments, obj_calc=None, mlp=1):
        bounds, bound_types = self.get_bounds_information()
        def objective(X):
            self.set_parameters(X)
            cost = 0
            predictions = odeint(self.get_dydt, self.get_initials(), time)
            for substrate_id, substrate_data in data:
                for time_point, truth in substrate_data:
                    prediction = predictions[int(time_point), substrate_id]
                    if obj_calc == None:
                        cost += (prediction - float(truth))**2
                    else:
                        cost += obj_calc(prediction, truth)
            return cost

        fitting_model = ga(function = objective,
                           dimension = len(bounds),
                           variable_type = bound_types,
                           variable_boundries = bounds,
                           algorithm_parameters = arguments
                           )
        
        fitting_model.run(set_function=ga.set_function_multiprocess(objective, n_jobs=mlp))
        return fitting_model.result.variable