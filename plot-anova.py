from sklearn.datasets import load_iris
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scikit_posthocs import posthoc_tukey
from scipy import stats
from statannotations.Annotator import Annotator

def plot_posthoc(path):
    temp = pd.read_csv(path)
    temp = temp.dropna()
    temp = temp.drop("index", axis=1)
    temp.rename(columns={"group": "cohort", "value": "Reads per million"}, inplace=True)

    tukey_df = posthoc_tukey(temp, val_col = "Reads per million", group_col = "cohort")
    remove = np.tril(np.ones(tukey_df.shape), k=0).astype("bool")
    tukey_df[remove] = np.nan

    molten_df = tukey_df.melt(ignore_index=False).reset_index().dropna()

    ax = sns.boxplot(data=temp, x=temp["cohort"], y=temp["Reads per million"], showfliers=False)

    pairs = [(i[1]["index"], i[1]["variable"]) for i in molten_df.iterrows()]
    p_values = [i[1]["value"] for i in molten_df.iterrows()]


    annotator = Annotator(
            ax, pairs, data=temp, x="cohort", y="Reads per million")
    annotator.configure(text_format="star", loc="inside")
    annotator.set_pvalues_and_annotate(p_values)
    plt.title("Title")

    plt.tight_layout()