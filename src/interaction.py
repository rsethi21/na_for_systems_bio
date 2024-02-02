class Interaction:
    """
    Summary
    -------
    The interaction class acts as a data store for parameters and rules that govern how two substrates in the network are affected.

    Attributes
    ----------
    identifier: str
        - a key unique to an interaction between two substrates and can be used to refer to an interaction
        - required
    affected: Substrate
        - substrate object of the substrate whose value will change in this interaction between two substrates
        - required
    effector: Substrate
        - substrate object of the substrate whose value will produce the change in the affected substrate's value
        - required    
    rate: Rate
        - rate object whose value acts as the proportionality constant between the affected and effector substrates
        - required
    effect: int (-1, 1, None)
        - if effect is -1, then the relationship between the substrate effector and affected is inverse meaning that as the effector substrates grows, the affected substrate falls; vice versa for effect of 1; None should only be used if the rate between the two substrates is the intrinsic phophorylation or dephosphoylation rates of the affected substrate and there is no feedback mechanism at play (no dissociation and hill coefficients)
        - optional in special cases (see logic in __init__ function)
    dissociation: Rate
        - rate object that encodes the 50% dissociation rate of two substrates; include if the interaction follows a feedback regulation mechanism; this package utilizes a goodwin oscillator mechanism to establish positive and negative feedback relationships; for more information see: https://pubmed.ncbi.nlm.nih.gov/32212037/
        - optional
    hill_coefficient: Rate
        - rate object that encodes the multiplicity of the two substrates (number of binding sites required); include if the interaction follows a feedback regulation mechanism; this package utilizes a goodwin oscillator mechanism to establish positive and negative feedback relationships; for more information see: https://pubmed.ncbi.nlm.nih.gov/32212037/
        - optional

    Methods
    -------
    __init__():
            initialize the appropriate interaction objects and inform of any exceptions
    update_rate():
            adjust the rate of the interaction; required for finding optimal network parameters
    get_interaction_paramters():
            obtain all associated rates for the interaction of interest

    """
    def __init__(self, identifier: str, affected: object, effector: object, rate: object, effect: int = None, dissociation: object = None, hill_coefficient: object = None):
        self.identifier = identifier
        self.s1 = affected
        self.s2 = effector
        self.Km = dissociation
        self.n = hill_coefficient
        self.rate = rate
        if effect in [1, -1]:
            self.effect = effect
        elif effect == None and dissociation == None and hill_coefficient == None or effect == None and rate.identifier in [affected.k.identifier, affected.r.identifier]:
            self.effect = effect
        else:
            raise ValueError("Effect defines whether the interaction up or downregulates so can only be 1 or -1")

    def update_rate(self, new_rate):
        """
        Summary
        -------
        updates the proportionality constant (encoded as the specified rate object) for the interaction between two substrates; changes rate object value so changes value for all uses of the same object
        
        Parameters
        ----------
        new_rate: value to change rate
        """
        self.rate.update_rate("value", new_rate)

    def get_interaction_parameters(self):
        """
        Summary
        -------
        returns a dictionary of the rates associated with this specific interaction
        
        Parameters
        ----------
        none
        
        Returns
        -------
        parameters: dictionary with keys as the rate identifiers and the values as the associated rate objects
        """
        parameters = {}
        if self.n != None:
            parameters[self.n.__getattribute__("identifier")] = self.n
        if self.Km != None:
            parameters[self.Km.__getattribute__("identifier")] = self.Km
        parameters[self.rate.__getattribute__("identifier")] = self.rate
        return parameters
        