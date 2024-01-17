class Substrate:

    def __init__(self, identifier: str, iv: float, type: str, k: float = None, r: float = None, mv: float = None, trs: list = None):
        self.identifier = identifier
        self.initial_value = iv
        self.current_value = iv
        self.type = type
        # if type.lower() in ["stimulus", "non-stimulus"]:
        #     self.type = type
        # else:
        #     raise ValueError("Type can only be one of the following: 'stimulus' or 'non-stimulus'")
        if self.type == "stimulus":
            self.max_value = mv
        self.time_ranges = trs
        if self.type == "non-stimulus":
            self.k = k
            self.r = r