class Substrate:
    """
    Summary
    -------
    The substrate class acts as a data store for information that pertains to a substrate protein or stimuli object

    Attributes
    ----------
    identifier: str
        - a key unique to a rate and can be used to refer to one specific rate
        - required
    initial_value: float
        - starting value of the substrate in the simulation
        - required
    current_value: float
        - present value of the substrate in the simulation
        - required
    type: str
        - this encodes the way that the substrate should be handled at integration time; can be either a stimulus or non-stimulus; a stimulus will be applied at certain time ranges to its max amplitude instantaneously as specified by the user
            - if a user specifies the type as stimulus, k and other will be ignored; optionally a user can include r in order to specify a decay rate for a signal however the stimulus may not reach the max value specified by the user
            - if a user specifies the type as non-stimulus, k, r will be required
            - more information about inputs provided in the readme
        - required
    max_value: float
        - max value that a stimuli will approach/reach
        - required
    time_ranges: boolean
        - this boolean informs the software of whether this rate will be pliable if the user decides to find optimal parameters using the fit function in network
        - optional; default False
    k: boolean
        - this boolean informs the software of whether this rate will be pliable if the user decides to find optimal parameters using the fit function in network
        - optional; default False
    r: boolean
        - this boolean informs the software of whether this rate will be pliable if the user decides to find optimal parameters using the fit function in network
        - optional; default False
    other: boolean
        - this boolean informs the software of whether this rate will be pliable if the user decides to find optimal parameters using the fit function in network
        - optional; default False

    Methods
    -------
    __init__():
            initialize the appropriate substrate objects and inform of any exceptions
    """
    def __init__(self, identifier: str, iv: float, type: str, k: float = None, r: float = None, mv: float = None, trs: list = None, other: str = None):
        
        self.identifier = identifier
        self.initial_value = iv
        self.current_value = iv
        self.type = type
        self.max_value = mv
        self.time_ranges = trs
        self.k = k
        self.r = r
        self.other = other
        self.reached = False
