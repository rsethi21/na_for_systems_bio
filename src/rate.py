class Rate:
    """
    Summary
    -------
    The rate class acts as a data store for information that pertains to a proportionality constant governing interaction between substrates

    Attributes
    ----------
    identifier: str
        - a key unique to a rate and can be used to refer to one specific rate
        - required
    value: float
        - substrate object of the substrate whose value will change in this interaction between two substrates
        - required
    fixed: boolean
        - this boolean informs the software of whether this rate will be pliable if the user decides to find optimal parameters using the fit function in network
        - optional; default False
    bounds: list of length 2
        - this is a set of bounds that the rate is allowed to be between and is important in network parameter tuning/fitting if the user decides to use the fit function
        - optional; defaut [1, 10]
    bounds_type: data type for bounds ("int", or "real")
        - this is the data type of the bounds for the rate, these define whether the rates will be integers or float when undergoing the fitting process if the user chooses
        - optional; default "int"

    Methods
    -------
    __init__():
            initialize the appropriate rate objects and inform of any exceptions
    update_rate():
            adjust the value of the rate; required for finding optimal network parameters; will store abs value since rates could be negative or positive depending on the interaction used in
    """
    def __init__(self, identifier: str, value: float, fixed: bool = False, bounds=[1, 10], bounds_type=int):
        self.identifier = identifier
        self.value = value
        self.fixed = fixed
        self.bounds = bounds
        self.bounds_type = bounds_type
    
    def update_rate(self, new_rate):
        """
        Summary
        -------
        updates the proportionality constant (encoded as the specified rate object); changes rate object value so changes value for all uses of the same object
        
        Parameters
        ----------
        new_rate: value to change rate
        """
        if not self.fixed:
            self.__setattr__("value", abs(new_rate))
