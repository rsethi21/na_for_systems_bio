import pandas as pd
from src.rate import Rate
from src.interaction import Interaction
from src.substrate import Substrate

# adding two values into one with 1 for each rate is total

def parse_rates(rates_csv_path):
    rates_df = pd.read_csv(rates_csv_path)
    rates = []
    for rate_index, rate in rates_df.iterrows():
        parameters = rate.to_dict()
        if pd.isna(parameters["bounds"]):
            del parameters["bounds"]
            del parameters["bounds_type"]
        else:
            if not pd.isna(parameters["bounds_type"]):
                convert_type = parameters["bounds_type"]
                if convert_type == "int":
                    convert_type = int
                elif convert_type == "real":
                    convert_type = float
            else:
                convert_type = float
            parameters["bounds"] = [convert_type(b) for b in parameters["bounds"].split(",")]
        rate_obj = Rate(**parameters)
        rates.append(rate_obj)
    return rates

def parse_substrates(substrates_csv_path, rates):
    substrates_df = pd.read_csv(substrates_csv_path)
    substrates = []
    for substrate_index, substrate in substrates_df.iterrows():
        parameters = substrate.to_dict()
        if not pd.isna(parameters["k"]):
            parameters["k"] = [r for r in rates if r.identifier == parameters["k"]][0]
        else:
            del parameters["k"]
        if not pd.isna(parameters["r"]):
            parameters["r"] = [r for r in rates if r.identifier == parameters["r"]][0]
        else:
            del parameters["r"]
        if not pd.isna(parameters["trs"]):
            range_strings = parameters["trs"].split(";")
            time_ranges = []
            for range_string in range_strings:
                time_ranges.append([int(t) for t in range_string.split(",")])
            parameters["trs"] = time_ranges
        else:
            del parameters["trs"]
        if not pd.isna(parameters["other"]):
            pass
        else:
            del parameters["other"]
        substrate_obj = Substrate(**parameters)
        substrates.append(substrate_obj)
    return substrates

def parse_interactions(interaction_csv_path, rates, substrates):
    interactions_df = pd.read_csv(interaction_csv_path)
    interactions = []
    for interaction_index, interaction in interactions_df.iterrows():
        parameters = interaction.to_dict()
        parameters["affected"] = [s for s in substrates if s.identifier == parameters["affected"]][0]
        parameters["effector"] = [s for s in substrates if s.identifier == parameters["effector"]][0]
        parameters["rate"] = [r for r in rates if r.identifier == parameters["rate"]][0]
        if not pd.isna(parameters["dissociation"]):
            parameters["dissociation"] = [r for r in rates if r.identifier == parameters["dissociation"]][0]
        else:
            del parameters["dissociation"]
        if not pd.isna(parameters["hill_coefficient"]):
            parameters["hill_coefficient"] = [r for r in rates if r.identifier == parameters["hill_coefficient"]][0]
        else:
            del parameters["hill_coefficient"]
        if not pd.isna(parameters["effect"]):
            pass
        else:
            del parameters["effect"]
        interaction_obj = Interaction(**parameters)
        interactions.append(interaction_obj)
    return interactions
