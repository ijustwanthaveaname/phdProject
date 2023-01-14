from json import load
from sklearn.datasets import load_iris
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier




def knn_iris():
    iris = load_iris()
    x_train, x_test, y_train, y_test = train_test_split(
        iris.data, iris.target, test_size=0.2, random_state=0)
    transfer = StandardScaler()
    x_train = transfer.fit_transform(x_train)
    x_test = transfer.transform(x_test)
    estimator = KNeighborsClassifier()
    # cv
    param_dict = {"n_neighbors": list(range(1, 12, 2))}
    estimator = GridSearchCV(estimator, param_grid=param_dict, cv=10)
    estimator.fit(x_train, y_train)
    y_predict = estimator.predict(x_test)
    print(y_test == y_predict)
    score = estimator.score(x_test, y_test)
    print(score)
    print(estimator.best_params_)
    print(estimator.best_score_)
    print(estimator.best_estimator_)
    print(estimator.cv_results_)


def decision_iris():
    iris = load_iris()
    x_train, x_test, y_train, y_test = train_test_split(iris.data, iris.target)
    estimator = DecisionTreeClassifier()
    estimator.fit(x_train, y_train)
    y_predict = estimator.predict(x_test)
    score = estimator.score(x_test, y_test)
    print(y_predict)
    print(score)


def rf():
    iris = load_iris()
    x_train, x_test, y_train, y_test = train_test_split(iris.data, iris.target)
    param_dict = {"n_estimators": [120, 200, 300, 500, 800, 1200],"max_depth": [5, 8 ,15, 25, 30]}
    estimator = RandomForestClassifier()
    estimator = GridSearchCV(estimator, param_grid=param_dict, cv=10)
    estimator.fit()
    
if __name__ == "__main__":
    # knn_iris()
    # decision_iris()
    iris = load_iris()
    x_train, x_test, y_train, y_test = train_test_split(
        iris.data, iris.target, test_size=0.2, random_state=0)
    transfer = StandardScaler()
    x_train = transfer.fit_transform(x_train)
    x_test = transfer.transform(x_test)
    estimator = KNeighborsClassifier()
    # cv
    param_dict = {"n_neighbors": list(range(1, 12, 2))}
    estimator = GridSearchCV(estimator, param_grid=param_dict, cv=10)
    estimator.fit(x_train, y_train)
    y_predict = estimator.predict(x_test)
    print(y_test == y_predict)
    score = estimator.score(x_test, y_test)
    print(score)
    print(estimator.best_params_)
    print(estimator.best_score_)
    print(estimator.best_estimator_)
    print(estimator.cv_results_)
