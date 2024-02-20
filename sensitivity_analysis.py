import pandas as pd
import argparse

pairs = {"pAKT": "AKT", "PI3K": "PI3Ks", "PTEN": "pPTEN", "GSK3B": "pGSK3B", "PTEN": "pPTEN", "pAKT": "AKT":, "PI3K": "PI3Ks", "PTEN": "pPTEN", "GSK3B": "pGSK3B", "PTEN": "pPTEN", "pAKT": "AKT", "PI3K": "PI3Ks", "PTEN": "pPTEN", "GSK3B": "pGSK3B", "PTEN": "pPTEN"}
path = "./output/data.csv"
df = pd.read_csv(path, index_col=[0,1])
new_cols = []

for col in df.columns:
    new_cols.append("")
