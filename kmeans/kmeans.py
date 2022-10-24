import numpy as np
import seaborn as sns

sns.set_theme()

data = np.random.rand(3,2)
data = sns.load_dataset("tips")
sns.scatterplot(data=data, x="total_bill", y="tip")
print(data)



