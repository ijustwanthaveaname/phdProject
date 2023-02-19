from sklearn.ensemble import RandomForestClassifier 
from sklearn.ensemble import VotingClassifier 
from sklearn.linear_model import LogisticRegression 
from sklearn.svm import SVC
from sklearn.datasets import load_iris
from sklearn.metrics import accuracy_score 
from sklearn.model_selection import train_test_split


iris = load_iris()
x_train, x_test, y_train, y_test = train_test_split(iris.data, iris.target, test_size=0.2, random_state=0)
log_clf = LogisticRegression() 
rnd_clf = RandomForestClassifier() 
svm_clf = SVC()
# create emmodel with 3 submodel
voting_clf = VotingClassifier(estimators=[('lr', log_clf), ('rf', rnd_clf), 
  ('svc', svm_clf)],voting='hard') 
voting_clf.fit(x_train, y_train)
from sklearn.metrics import accuracy_score 
    for clf in (log_clf, rnd_clf, svm_clf, voting_clf): 
        clf.fit(x_train, y_train) 
        y_pred = clf.predict(x_test) 
        print(clf.__class__.__name__, accuracy_score(y_test, y_pred))
        

# bagging, 多个同分类器，有放回抽样不同训练集，pasting是无放回抽样
from sklearn.ensemble import BaggingClassifier 
from sklearn.tree import DecisionTreeClassifier
bag_clf = BaggingClassifier(DecisionTreeClassifier(), n_estimators=500,        
  max_samples=100, bootstrap=True, n_jobs=-1) 
bag_clf.fit(x_train, y_train) 
y_pred = bag_clf.predict(x_test)
accuracy_score(y_test, y_pred)

# Out-of-bag 评价
bag_clf = BaggingClassifier(DecisionTreeClassifier(), n_estimators=500,bootstrap=True, n_jobs=-1, oob_score=True)
bag_clf.fit(x_train, y_train) 
bag_clf.oob_score_ 
# evaluation in test set
y_pred = bag_clf.predict(x_test) 
accuracy_score(y_test, y_pred) 
bag_clf.oob_decision_function_

# 随机贴片与随机子空间
# 两个超参数max_features和bootstrap_features控制。
# 他们的工作方式和max_samples和bootstrap一样，但这是对于特征采样而不是实例采样。因此，每一个分类器都会被在随机的输入特征内进行训练。
# 当你在处理高维度输入下此方法尤其有效。对训练实例和特征的采样被叫做随机贴片。
# 保留了所有的训练实例（例如bootstrap=False和max_samples=1.0），但是对特征采样（bootstrap_features=True并且/或者max_features小于 1.0）叫做随机子空间。

# feature importance
from sklearn.datasets import load_iris 
iris = load_iris() 
rnd_clf = RandomForestClassifier(n_estimators=500, n_jobs=-1) 
rnd_clf.fit(iris["data"], iris["target"]) 
for name, score in zip(iris["feature_names"], rnd_clf.feature_importances_): 
    print(name, score) 