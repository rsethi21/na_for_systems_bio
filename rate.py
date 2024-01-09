class Rate:
    def __init__(self, identifier: str, value: float, fixed: bool = False):
        self.identifier = identifier
        self.value = value
        self.fixed = fixed
    
    def update_rate(self, new_rate):
        if not self.fixed:
            self.__setattr__("value", abs(new_rate))