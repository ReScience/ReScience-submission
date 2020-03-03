import pandas as pd
import numpy as np

from sklearn.model_selection import cross_validate
from sklearn.metrics import accuracy_score, precision_score, f1_score, roc_auc_score, make_scorer
from imblearn.metrics import specificity_score, sensitivity_score
from pathlib import Path



def methods_classification(n_neighbors=3,
                           kernel_a='linear', kernel_b='rbf', gamma='auto',
                           max_depth=5,
                           n_estimators=10, random_state=42, max_features=1,
                           max_iter=5000):
    """

    Parameters
    ----------
    n_neighbors
    kernel_a
    kernel_b
    gamma
    max_depth
    n_estimators
    random_state
    max_features
    max_iter
    """
    # K-NN
    from sklearn.neighbors import KNeighborsClassifier
    # SVM (linear and radial)
    from sklearn import svm
    # Decision tree
    from sklearn.tree import DecisionTreeClassifier
    # random forests
    from sklearn.ensemble import RandomForestClassifier, VotingClassifier
    # MLPClassifier
    from sklearn.neural_network import MLPClassifier
    # AdaBoostClassifier
    from sklearn.ensemble import AdaBoostClassifier
    # GaussianNB
    from sklearn.naive_bayes import GaussianNB

    knn = KNeighborsClassifier(n_neighbors=n_neighbors)  # 1
    SVM1 = svm.SVC(kernel=kernel_a)  # 2
    SVM2 = svm.SVC(kernel=kernel_b, gamma=gamma)  # 3
    DT = DecisionTreeClassifier(max_depth=max_depth)  # 4
    RF = RandomForestClassifier(
        n_estimators=n_estimators, random_state=random_state, max_features=max_features)  # 5
    MLP = MLPClassifier(max_iter=max_iter)  # 6
    ADB = AdaBoostClassifier(random_state=random_state)  # 7
    from sklearn.naive_bayes import GaussianNB
    GaussianNB = GaussianNB()  # 8

    ensemble = VotingClassifier(estimators=[('KNN', knn), ('SVM1', SVM1), ('SVM2', SVM2),
                                            ('dt', DT), ('RF', RF), ('MLP', MLP), 
                                            ('ADB', ADB),
                                            ('GaussianNB', GaussianNB)], voting='hard')

    classifiers = [('KNN', knn), ('SV1', SVM1), ('SVM2', SVM2), ('DT', DT), ('RF', RF),
                   ('MLP', MLP), ('ADB', ADB), ('GaussianNB', GaussianNB),
                   ("Ensemble", ensemble)]

    return classifiers


def makeBF(merge):
    '''
    Function to mark in the latex which were the best results per line.
     After execution it is still necessary to replace "\ t" with "\ t", without regular expression
     "\ n" by "\\\\ \ n" with regular expression. Replace in all files.
    '''
    merge_1 = merge.reset_index()
    merge_1['AVG'] = np.average(merge_1.drop(["m", "Ensemble"], 1), axis=1)
    merge_1 = merge_1.round(3)
    accumulator = pd.DataFrame()
    for m in merge_1['m']:
        tmp_row = merge_1[merge_1['m'] == m]
        tmp_row = tmp_row.drop("m", 1)
        idx_max = tmp_row.values.argmax()

        tmp_row = tmp_row.astype(str)
        tmp_row.iloc[0][idx_max] = "\ textbf{" + tmp_row.iloc[0][idx_max] + "}"
        accumulator = accumulator.append(tmp_row)

    csv = (pd.concat([merge_1['m'], accumulator], axis=1))
    csv = csv[['m', 'KNN', 'SV1', 'SVM2', 'DT', 'RF',
               'MLP', 'ADB', 'GaussianNB', 'AVG', 'Ensemble']]
    return csv


def save_metrics(classifiers, base_fold='../data/boon/featureDataSet',
                 type_loss='l1', base_save='../data/boon/save/'):
    """
    REALLY unoptimized function to save the metrics.

    """
    fold = Path(base_save)

    range_values = [2, 4, 8, 16, 32, 64, 128, 256]
    index = pd.DataFrame(range_values, columns=["m"])
    index.index = index['m']
    merge_acc = pd.DataFrame(index.drop("m", 1))
    merge_pre = pd.DataFrame(index.drop("m", 1))
    merge_spe = pd.DataFrame(index.drop("m", 1))
    merge_sen = pd.DataFrame(index.drop("m", 1))
    merge_f1 = pd.DataFrame(index.drop("m", 1))
    merge_auc = pd.DataFrame(index.drop("m", 1))

    for classifier_name, classifier in classifiers:
        base_acc = pd.DataFrame([], columns=[classifier_name])
        base_pre = pd.DataFrame([], columns=[classifier_name])
        base_spe = pd.DataFrame([], columns=[classifier_name])
        base_sen = pd.DataFrame([], columns=[classifier_name])
        base_f1 = pd.DataFrame([], columns=[classifier_name])
        base_auc = pd.DataFrame([], columns=[classifier_name])

        for m in range_values:

            nameTrain = base_fold+"/train_"+str(m)+"_"+type_loss+".parquet"
            X_train = pd.read_parquet(
                nameTrain, engine='pyarrow').drop(["class"], 1)
            Y_train = pd.read_parquet(nameTrain, engine='pyarrow')['class']

            nameTest = base_fold+"/test_"+str(m)+"_"+type_loss+".parquet"
            X_test = pd.read_parquet(
                nameTest, engine='pyarrow').drop(["class"], 1)
            Y_test = pd.read_parquet(nameTest, engine='pyarrow')['class']

            X = X_train.append(X_test)
            Y = Y_train.append(Y_test)

            scoring = {'accuracy': make_scorer(accuracy_score),
                       'precision': make_scorer(precision_score),
                       'specificity': make_scorer(specificity_score),
                       'sensibility': make_scorer(sensitivity_score),
                       'f1_score': make_scorer(f1_score),
                       'roc_auc_score': make_scorer(roc_auc_score)}

            scores = cross_validate(classifier, X, Y, cv=5, scoring=scoring)

            base_acc = base_acc.append(pd.DataFrame(
                [np.average(scores['test_accuracy'])], columns=[classifier_name]))
            base_pre = base_pre.append(pd.DataFrame(
                [np.average(scores['test_precision'])], columns=[classifier_name]))
            base_spe = base_spe.append(pd.DataFrame(
                [np.average(scores['test_specificity'])], columns=[classifier_name]))
            base_sen = base_sen.append(pd.DataFrame(
                [np.average(scores['test_sensibility'])], columns=[classifier_name]))
            base_f1 = base_f1.append(pd.DataFrame(
                [np.average(scores['test_f1_score'])], columns=[classifier_name]))
            base_auc = base_auc.append(pd.DataFrame(
                [np.average(scores['test_roc_auc_score'])], columns=[classifier_name]))

        base_acc.index = index['m']
        base_pre.index = index['m']
        base_spe.index = index['m']
        base_sen.index = index['m']
        base_f1.index = index['m']
        base_auc.index = index['m']

        merge_acc = pd.concat([merge_acc, base_acc], axis=1)
        merge_pre = pd.concat([merge_pre, base_pre], axis=1)
        merge_spe = pd.concat([merge_spe, base_spe], axis=1)
        merge_sen = pd.concat([merge_sen, base_sen], axis=1)
        merge_f1 = pd.concat([merge_f1, base_f1], axis=1)
        merge_auc = pd.concat([merge_auc, base_auc], axis=1)

    if (~fold.exists()):
        fold.mkdir(parents=True, exist_ok=True)

    makeBF(merge_acc).to_csv(path_or_buf=base_save+"accuracy.csv",
                             sep='&', encoding='utf-8', index=False, header=False)
    makeBF(merge_pre).to_csv(path_or_buf=base_save+"precision.csv",
                             sep='&', encoding='utf-8', index=False, header=False)
    makeBF(merge_spe).to_csv(path_or_buf=base_save+"specificity.csv",
                             sep='&', encoding='utf-8', index=False, header=False)
    makeBF(merge_sen).to_csv(path_or_buf=base_save+"sensibility.csv",
                             sep='&', encoding='utf-8', index=False, header=False)
    makeBF(merge_f1).to_csv(path_or_buf=base_save+"f1_score.csv",
                            sep='&', encoding='utf-8', index=False, header=False)
    makeBF(merge_auc).to_csv(path_or_buf=base_save+"roc_auc_score.csv",
                             sep='&', encoding='utf-8', index=False, header=False)

    return merge_acc, merge_pre, merge_spe, merge_sen, merge_f1, merge_auc