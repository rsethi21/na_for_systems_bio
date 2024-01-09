class Interaction:

    def __init__(self, identifier: str, affected: object, effector: object, rate: object, effect: int, dissociation: object = None, hill_coefficient: object = None):
        self.identifier = identifier
        self.s1 = affected
        self.s2 = effector
        self.Km = dissociation
        self.n = hill_coefficient
        self.rate = rate
        if self.effect in [1, -1]:
            self.effect = effect
        else:
            raise ValueError("Effect defines whether the interaction up or downregulates so can only be 1 or -1")

    def update_rate(self, new_rate):
        self.rate.update_rate("value", new_rate)

    def get_interaction_parameters(self):
        parameters = {}
        if self.Km != None:
            parameters[self.n.__getattribute__("identifier")] = self.n
        if self.n != None:
            parameters[self.Km.__getattribute__("identifier")] = self.Km
        parameters[self.rate.__getattribute__("identifier")] = self.rate
        return parameters
        