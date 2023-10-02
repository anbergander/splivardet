import pandas as pd
from bioinfokit.analys import stat

temp = pd.read_csv('NM_006579.3_EBP.csv')
temp = temp.dropna()
temp = temp.drop('index', axis=1)
temp = temp.rename(columns={'group': 'cohort'})

#temp.cohort = temp.cohort.astype("str")
#temp.cohort = pd.Categorical(temp.cohort)
#temp['code'] = temp.cohort.cat.codes
#temp['cohort'] = temp['cohort'].astype(str)
#temp = temp.astype({'code':'int'})


print(temp.dtypes)

res = stat()
res.tukey_hsd(df=temp, res_var='value', xfac_var='cohort', anova_model='value ~ C(cohort)')

