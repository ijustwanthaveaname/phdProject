# 数据集获取
from sklearn import datasets
# 数据集划分
from sklearn.model_selection import train_test_split
# 特征抽取
from sklearn.feature_extraction import DictVectorizer
from sklearn.feature_extraction.text import CountVectorizer, TfidfVectorizer
# 预处理，同样用fit_transfrom,输入为numpy数组
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import jieba
# 特征选择
from sklearn.feature_selection import VarianceThreshold
# PCA
from sklearn.decomposition import PCA
# 聚类
from sklearn.cluster import KMeans 
from sklearn.metrics import silhouette_score
# 计算相关系数
from scipy.stats import pearsonr
# fetch_* 大数据集从网上下载
# datasets.fetch_20newsgroups()
# load_* 小数据集本地加载
iris = datasets.load_iris()
print(iris.DESCR)

# 划分数据集
x_train, x_test, y_train, y_test = train_test_split(
    iris.data, iris.target, test_size=0.2, random_state=22)

# 字典特征提取
data = [{'city': c, 'temperature': t}
        for c, t in zip(("北京", "上海", "深圳"), (100, 60, 30))]
transfer = DictVectorizer(sparse=False)  # 是否为稀疏矩阵
data_new = transfer.fit_transform(data)
transfer.get_feature_names()

# 文本特征提取
data = ["Life is short, i like like python",
        "Life is too long, i dislike python"]
transfer = CountVectorizer()
data_new = transfer.fit_transform(data)
print(data_new.toarray())
print(transfer.get_feature_names())


def cut_word(text):
    text = " ".join(list(jieba.cut(text)))
    return text


data = ["你是傻逼我是帅哥", "美女看看我"]
data_new = []
for sent in data:
    data_new.append(cut_word(sent))
    print(data_new)

transfer = CountVectorizer(stop_words=["傻逼"])
data_final = transfer.fit_transform(data_new)
print(data_final.toarray())
print(transfer.get_feature_names())

transfer = TfidfVectorizer(stop_words=["傻逼"])
data_final = transfer.fit_transform(data_new)
print(data_final.toarray())
print(transfer.get_feature_names())
transfer = PCA(n_components=2)
data_new = transfer.fit_transform(iris.data)
km = KMeans(n_clusters=3)
km.fit(data_new)
pred = km.predict(data_new)
silhouette_score(data_new, pred)