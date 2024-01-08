class Substrate:

    def __init__(self, identifier: str, iv: float, mv: float = None, trs: list = None):
        self.identifier = identifier
        self.initial_value = iv
        self.current_value = iv
        self.max_value = mv
        self.time_ranges = trs