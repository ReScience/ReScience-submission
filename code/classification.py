import pandas as pd
import numpy as np
 
from sklearn.model_selection import cross_validate
from sklearn.metrics import accuracy_score, precision_score, f1_score, roc_auc_score,make_scorer
from imblearn.metrics import specificity_score, sensitivity_score


#K-NN
from sklearn.neighbors import KNeighborsClassifier
#SVM (linear and radial)
from sklearn import svm
#Decision tree
from sklearn.tree import DecisionTreeClassifier
#random forests
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
#MLPClassifier 
from sklearn.neural_network import MLPClassifier
#AdaBoostClassifier
from sklearn.ensemble import AdaBoostClassifier
#GaussianNB
from sklearn.naive_bayes import GaussianNB

def methods_classification(n_neighbors=3, 
                           kernel_a='linear',kernel_b='rbf', gamma='auto', 
                           max_depth=5, 
                           n_estimators=10, random_state=42, max_features=1, 
                           max_iter=5000):

    knn = KNeighborsClassifier(n_neighbors=n_neighbors)  # 1
    SVM1 = svm.SVC(kernel=kernel_a)  # 2
    SVM2 = svm.SVC(kernel=kernel_b, gamma=gamma)  # 3
    DT = DecisionTreeClassifier(max_depth=max_depth)  # 4
    RF = RandomForestClassifier(
        n_estimators=n_estimators, random_state=random_state, max_features=max_features)  # 5
    MLP = MLPClassifier(max_iter=max_iter)  # 6
    ADB = AdaBoostClassifier(random_state=random_state)  # 7
    GaussianNB = GaussianNB()  # 8