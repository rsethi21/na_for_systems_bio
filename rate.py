class Rate:
    def __init__(self, identifier: str, value: float, fixed: bool = False, bounds=[1, 10]):
        self.identifier = identifier
        self.value = value
        self.fixed = fixed
        self.bounds = bounds
    
    def update_rate(self, new_rate):
        if not self.fixed:
            self.__setattr__("value", abs(new_rate))