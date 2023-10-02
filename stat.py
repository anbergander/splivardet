from sklearn.datasets import load_iris
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scikit_posthocs import posthoc_tukey
from scipy import stats
from statannotations.Annotator import Annotator


temp = pd.read_csv('NM_006579.3_EBP.csv')
temp = temp.dropna()
temp = temp.drop("index", axis=1)
temp.rename(columns={"group": "cohort", "value": "Reads per million"}, inplace=True)
print(data)
print(temp.head())


# First we do a oneway ANOVA as implemented in SciPy


#tukey_df = posthoc_tukey(iris_df, val_col="sepal length (cm)", group_col="species")
#print(tukey_df)
tukey_df = posthoc_tukey(temp, val_col = "Reads per million", group_col = "cohort")
print(tukey_df)

remove = np.tril(np.ones(tukey_df.shape), k=0).astype("bool")
tukey_df[remove] = np.nan

molten_df = tukey_df.melt(ignore_index=False).reset_index().dropna()
print(molten_df)

ax = sns.boxplot(data=temp, x=temp["cohort"], y=temp["Reads per million"], showfliers=False)

pairs = [(i[1]["index"], i[1]["variable"]) for i in molten_df.iterrows()]
p_values = [i[1]["value"] for i in molten_df.iterrows()]

print(pairs)

annotator = Annotator(
    ax, pairs, data=temp, x="cohort", y="Reads per million")
annotator.configure(text_format="star", loc="inside")
annotator.set_pvalues_and_annotate(p_values)
plt.title("Title")

plt.tight_layout()