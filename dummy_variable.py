from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd


df = pd.DataFrame({"weight": [50, 55, 30, 44],
                   "age": [20, 30, 25, 40],
                   "sex": ["male"] * 2 + ["female"] * 2})
# get dummy variable
# you can set: drop_first = True
df_dummy = pd.get_dummies(df, prefix=["sex"], columns=["sex"])
# you can set: drop_first = True
df_dummy = pd.get_dummies(df, prefix=["sex"], columns=["sex"], drop_first=True)
# another way to get dummy variable
y = df["weight"]
x = df[["age", "sex"]]
oh = OneHotEncoder(sparse=False)  # you can set drop="first"
sex_transformed = oh.fit_transform(x["sex"].to_numpy().reshape(-1, 1))
sex_transdf = pd.DataFrame(sex_transformed, columns=oh.get_feature_names_out())
df_dummy2 = pd.concat([df, sex_transdf], axis=1).drop(["sex"], axis=1)
# fit model
Y = df_dummy["weight"].to_numpy().reshape(-1, 1)
X = df_dummy.loc[:, df_dummy.columns != "weight"]
lr = LinearRegression()
lr.fit(X, Y)
# df.loc[:, ~df.columns.isin(["age", "sex"])], exclude columns of age and sex
print(lr.intercept_, lr.coef_, lr.score(X, Y))

# use statsmodels to fit linear regression
mod = smf.ols(formula="weight ~ age + sex", data=df)
res = mod.fit()
print(res.summary())
print(res.params)
