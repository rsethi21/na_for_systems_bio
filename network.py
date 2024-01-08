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

    def freeze_parameters(self, parameter_id_tuple):
        pass

    def get_dydt(self, y0, t):
        dydt = {substrate.identifier: {} for substrate in self.substrates}
        for id, interaction in self.interactions.items():
            substrate_id = interaction.__getattribute__("s1").__getattribute__("identifier")
            rate_id = interaction.__getattribute__("rate").__getattribute__("identifier")
            substrate2 = interaction.__getattribute__("s2").__getattribute__("identifier")
            Km = interaction.__getattribute__("Km")
            n = interaction.__getattribute__("n")
            if  Km != None and n != None:
                if rate_id not in dydt[substrate_id].keys():
                    dydt[substrate_id][rate_id] = [(substrate2, Km.__getattribute__("identifier"), n.__getattribute__("identifier"))]
                else:
                    dydt[substrate_id][rate_id].append((substrate2, Km.__getattribute__("identifier"), n.__getattribute__("identifier")))
            else:
                if rate_id not in dydt[substrate_id].keys():
                    dydt[substrate_id][rate_id] = [substrate2.__getattribute__("value")]
                else:
                    dydt[substrate_id][rate_id].append(substrate2.__getattribute__("value"))