# This script convert the first .h5 dataframe, passed in the command line, in a LISA type dataframe with name given by the first argument name starting with LISA

import numpy as np
import sys
import pandas as pd
import pickle

df_key = sys.argv[2]
df_nm = sys.argv[1]

BHCat = pd.read_hdf(df_nm, df_key)

BHCat.to_pickle(df_nm[:-3]+'.pkl')
