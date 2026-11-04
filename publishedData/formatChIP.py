#!/bin/python
## Use sci enviroment

import sys
import pandas as pd
from pathlib import Path

# path=Path(sys.argv[1])
# in_file = path / sys.argv[2]
# out_file = path / sys.argv[3]

in_file = Path(sys.argv[1])
out_file = Path(sys.argv[2])

df = pd.read_csv(in_file)
df = df.iloc[:, 1:]
df.loc[:, "Chromosome"] = "chr" + df.loc[:, "Chromosome"].astype(str)
df = df.replace(",", "", regex=True)

df.to_csv(out_file, index=False, sep='\t', header=False)