from sklearn.datasets import load_boston
from sklearn.linear_model import LinearRegression, SGDRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


def lr():
    boston = load_boston()
    x_train, x_test, y_train, y_test = train_test_split(
        boston.data, boston.target, random_state=22)
    transfer = StandardScaler()
    x_train = transfer.fit_transform(x_train)
    x_test = transfer.transform(x_test)
    estimator = LinearRegression()
    estimator.fit(x_train, y_train)
    y_pred = estimator.predict(x_test)
    error = mean_squared_error(y_pred, y_test)
    print(error)


def sgd():
    boston = load_boston()
    x_train, x_test, y_train, y_test = train_test_split(
        boston.data, boston.target, random_state=22)
    transfer = StandardScaler()
    x_train = transfer.fit_transform(x_train)
    x_test = transfer.transform(x_test)
    estimator = SGDRegressor(learning_rate="constant",
                             eta0=0.01, max_iter=10000, random_state=22)
    estimator.fit(x_train, y_train)
    y_pred = estimator.predict(x_test)
    error = mean_squared_error(y_pred, y_test)
    print(error)
    print(estimator.coef_)
    print(estimator.intercept_)


if __name__ == "__main__":
    lr()
    sgd()
