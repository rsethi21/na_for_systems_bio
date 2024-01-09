import pandas as pd

class Network:
    def __init__(self, identifier: str, rates: list, interactions: list, substrates: list):
        self.identifier = identifier
        self.interactions = {interaction.identifier: interaction for interaction in interactions}
        self.substrates = {substrate.identifier: substrate for substrate in substrates}
        self.parameters = self.get_all_parameters()
        self.parsed_interactions = self.parse_interactions()
    
    def get_all_parameters(self):
        parameters = {}
        for interaction in list(self.interactions.values()):
            parameters.update(interaction.get_interaction_parameters())
        return parameters
    
    def set_currents(self, new_currents):
        for substrate, new_current in zip(list(self.substrates.values()), new_currents):
            substrate.__setattr__("current_value", new_current)

    def get_initials(self):
        y0 = []
        for substrate in list(self.substrates.values()):
            y0.append(substrate.__getattribute__())
        return y0

    def freeze_parameters(self, parameter_id_tuple):
        pass

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
                        temporary_term = self.parameters[parameter_id]*associated_interactions[0][1]
                    else:
                        temporary_term = terms[parameter_id]
                    for interaction in associated_interactions:
                        if len(interaction) == 2:
                            temporary_term *= self.substrates[interaction[0]].__getattribute__("current_value")
                        else:
                            substrate = self.substrates[interaction[0]].__getattribute__("current_value")
                            Km = self.parameters[interaction[2]]
                            n = self.parameters[interaction[3]]
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
    
    def graph(self, time):
        pass