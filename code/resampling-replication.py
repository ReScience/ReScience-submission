#------------------------------------------------------------------------------

import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import pandas as pd

np.random.seed(1) # for reproducable results

#------------------------------------------------------------------------------

# Define a function to calculate prediction errors from LDA model
def peLDA(trainPA, trainEV, testPA, testEV):
    # Create and fit model
    LDA.fit(trainEV, trainPA)
    # Predict against test data
    predicts = LDA.predict(testEV)
    # Calculate and return the prediction errors
    return(np.abs(testPA - predicts))

# Define a function to create bootstrap samples from the dataframe
def bootSample(df):
    # Create bootstrap sample
    boot = np.random.choice(np.arange(df.shape[0]), size=df.shape[0])
    # Identify samples not in bootstrap sample
    xboot = np.where(np.isin(df["ID"], boot, invert=True))[0]
    # Create train data by selecting bootstrap sample
    trainPA = df.iloc[boot,1]
    trainEV = df.iloc[boot,2:]
    # Create test data by selecting non-bootstrap sample
    testPA = df.iloc[xboot,1]
    testEV = df.iloc[xboot,2:]
    return(boot, xboot, trainPA, trainEV, testPA, testEV)

# Define a function to compute f-fold cross-validated prediction errors
def kfcv(df, k):
    # Randomly reorder the dataframe
    rdf = df.sample(frac=1)
    # Create empty list to hold mean for each fold   
    kfcvfold = []
    # Create folds for each sample
    sampleFolds = np.repeat(range(k), int(rdf.shape[0] / k))
    # Calculate the prediction error for each sample
    for i in range(rdf.shape[0]):
        # Identify the fold the sample belongs to
        fold = sampleFolds[i]
        # Identify indices of samples not in the sample's fold
        xfold = np.where(np.isin(sampleFolds, fold, invert=True))[0]
        # Create training and testing data
        trainPA = rdf.iloc[xfold,1]
        trainEV = rdf.iloc[xfold,2:]
        testPA = rdf.iloc[i,1]
        testEV = rdf.iloc[i,2:].values.reshape(1,-1)
        # Calculate prediction errors, and append mean to the fold results
        errors = peLDA(trainPA, trainEV, testPA, testEV)
        kfcvfold.append(np.mean(errors))
    # Return the mean prediction error across all folds
    return(np.mean(kfcvfold))    

# Define a function to compute hold-out validated prediction errors
def hv(df, p):
    # Calculate index that makes split point
    split = int(df.shape[0] * p)
    # Define a function to return samples
    def hsample(df, split):
        # Randomly reorder the dataframe
        rdf = df.sample(frac=1)
        # Create training and testing data
        testPA = rdf.iloc[:split,1]
        testEV = rdf.iloc[:split,2:]
        trainPA = rdf.iloc[split:,1]
        trainEV = rdf.iloc[split:,2:]
        return(testPA, testEV, trainPA, trainEV)
    # Create the samples
    testPA, testEV, trainPA, trainEV = hsample(df, split)
    # Replace the sample if training data is all presences or absences
    while np.sum(trainPA) == s or np.sum(trainPA) == 0:
       testPA, testEV, trainPA, trainEV = hsample(df, split)
    # Calculate prediction errors
    errors = peLDA(trainPA, trainEV, testPA, testEV)
    # Return the mean prediction error across all samples
    return(np.mean(errors))
    
#------------------------------------------------------------------------------

# Create empty lists to hold results from each computational experiment
resubstitution = []
jackknife = []
cv10fold = []
cv3fold = []
bootstrap = []
hov = []
hocvp368 = []
hocvp200 = []

# For 1000 computational experiments
for ce in range(1000):
    
    print(ce)
    
    # Specify sample size
    s = 30
    # Create an empty dataframe for data
    colnames = ["ID","pa","ev1","ev2","ev3","ev4",
                "ev5","ev6","ev7","ev8","ev9","ev10"]
    df = pd.DataFrame(columns = colnames)
    # Genreate sample ID and random presence and absence
    df.iloc[:,0] = range(s)
    df.iloc[:,1] = np.random.binomial(1, 0.5, size=s)
    # Generate normally distributed random explanatory variables
    for i in range(10):
        df.iloc[:,i+2] = np.random.normal(size=s)
            
    #--------------------------------------------------------------------------
    # RESUBSTITUTION
    #--------------------------------------------------------------------------
    
    # Create training and testing data
    trainPA = df.iloc[:,1]
    trainEV = df.iloc[:,2:]
    testPA = trainPA
    testEV = trainEV
    # Calculate prediction errors
    errors = peLDA(trainPA, trainEV, testPA, testEV)
    # Calculate resubsititution prediction error, and append to results
    resubstitution.append(np.mean(errors))
    
    #--------------------------------------------------------------------------
    # BOOTSTRAP
    #--------------------------------------------------------------------------
    
    # Create some empty arrays to hold results from each bootstrap
    xBootSamples = 0
    xBootErrors = 0
    # For each bootstrap sample
    for b in range(200):
        # Create the bootstrap samples
        boot, xboot, trainPA, trainEV, testPA, testEV = bootSample(df)
        # Replace the sample if training data is all presences or absences
        while np.sum(trainPA) == s or np.sum(trainPA) == 0:
            boot, xboot, trainPA, trainEV, testPA, testEV = bootSample(df)
        # Calculate prediction errors, and append the bootstrap sample results
        errors = peLDA(trainPA, trainEV, testPA, testEV)
        xBootSamples = xBootSamples + len(xboot)
        xBootErrors = xBootErrors + np.sum(errors)
    # Calculate the mean across all samples, and append to results
    bootstrap.append(xBootErrors / xBootSamples)
    
    #--------------------------------------------------------------------------
    # K-FOLD CROSS-VALIDATION
    #--------------------------------------------------------------------------
    
    # Do 10-fold cross-validation
    cv10fold.append(kfcv(df, k = 10))
    # Do 3-fold cross-validation
    cv3fold.append(kfcv(df, k = 3))
    # Do jackknife cross-validation
    jackknife.append(kfcv(df, k = s))
    
    #--------------------------------------------------------------------------
    # HOLD-OUT CROSS-VALIDATION
    #--------------------------------------------------------------------------
    
    # Do hold-out cross-validation for single split
    hov.append(hv(df, p=0.368))
    # Do hold-out cross-validations for multiple splits
    hcv_H = []
    for H in range(200):
        hcv_H.append(hv(df, p=0.368))
    hocvp368.append(np.mean(hcv_H))
    hcv_H = []
    for H in range(200):
        hcv_H.append(hv(df, p=0.200))
    hocvp200.append(np.mean(hcv_H))
    
#------------------------------------------------------------------------------

# Plot smooted prediction errors
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)

from scipy import stats
x = np.linspace(0,1,100)
bw = 1
resampMethods = [r"$Err^{R}$", r"$Err^{H=1}_{p=0.368}$", r"$Err^{H=200}_{p=0.368}$", 
                 r"$Err^{H=200}_{p=0.200}$", r"$Err^{k=3}$", r"$Err^{k=10}$", 
                 r"$Err^{k=n}$", r"$Err^{B=200}$"]
predErrors = [resubstitution, hov, hocvp368, hocvp200, cv3fold, cv10fold, 
              jackknife, bootstrap]
orange = (0.90,0.62,0)
blue = (0.34,0.71,0.91)
teal = (0,0.62,0.45)
linecolour = [blue, teal, teal, teal, orange, orange, orange, "black"]
linestyle = ['-', '-', '--', ':', '-', '--', ':', '-', '-']
plt.figure(figsize=(7, 5), dpi=150)
plt.xlabel('Prediction error estimate')
plt.ylabel('Kernel density')
for rm in range(8):
    kernel = stats.gaussian_kde(predErrors[rm], bw_method=bw)
    KDE = kernel(x)
    plt.plot(x, KDE, c=linecolour[rm], linestyle=linestyle[rm], linewidth=1)
plt.legend(resampMethods, loc='upper right', frameon=False)    
plt.axvline(x=0.5, ymin=0, ymax = 200, linewidth=1, color='lightgrey')
plt.savefig('resampling-results.png', bbox_inches='tight')

#------------------------------------------------------------------------------

# Export summary data as text file
colnames = ["method", "mean-pe", "std-pe"]
resampMethods = ["R", "H=1,p=0.368", "H=200,p=0.368", "H=200,p=0.200", "K=3", 
                 "K=10", "K=n", "B=200"]
resultsdf = pd.DataFrame(columns = colnames)
resultsdf.iloc[:,0] = resampMethods
for rm in range(8):
    resultsdf.iloc[rm,1] = '{:.3f}'.format(np.round(np.mean(predErrors[rm]), 3))
    resultsdf.iloc[rm,2] = '{:.3f}'.format(np.round(np.std(predErrors[rm]), 3))
resultsdf.to_csv("resampling-results.txt", sep='\t', index=False)

#------------------------------------------------------------------------------
