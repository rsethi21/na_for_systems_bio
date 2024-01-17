class Rate:
    def __init__(self, identifier: str, value: float, fixed: bool = False, bounds=[1, 10], bounds_type=int):
        self.identifier = identifier
        self.value = value
        self.fixed = fixed
        self.bounds = bounds
        self.bounds_type = bounds_type
    
    def update_rate(self, new_rate):
        if not self.fixed:
            self.__setattr__("value", abs(new_rate))