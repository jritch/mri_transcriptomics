import pandas
import numpy as np

import seaborn.apionly
iris = seaborn.apionly.load_dataset("iris")[0:10]

iris.sort_values(by="petal_width")
