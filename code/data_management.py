from pathlib import Path
from wget import download
from zipfile import ZipFile
from pandas import read_csv
from numpy import zeros, ones, concatenate, array
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split

def zip_with_unique(base, list_suffix):
    """ Auxiliary function to generate a paired 
    list considering a unique element.
    
    An adaptation of the convolution function (`zip`) 
    to map a single value in a sequence tuple.  
    This mapping is a surjective-only.

    An adaptation of the "scalar product" 
    to make a zip with vectors of different sizes. 
    The first vector must have 1 item, 
    while the second must-have n items.
    
    Parameters
    ----------
    base: array-like, shape (1, )
        One (1) base prefix that will be paired with the suffixes.
        
    list_suffix : array-like, shape (n_suffix,)
        A suffix list that will be paired with one prefix.
    
    Returns
    -------
    list: array-like, shape (n_suffix)
        Base added with the suffix.
    """
    
    return list(base + suffix for suffix in list_suffix)


def download_bonn(path_data='data/boon/') -> [str]:
    """
    Adapted from mne-tools.github.io/mne-features/auto_examples/plot_seizure_example.html
    Code changes were:
        * Adding more folders;
        * Control for folder creation;
    :rtype: [str]
    """
    fold = Path(path_data)
    child_fold = ['setA', 'setB', 'setC', 'setD', 'setE']
    base_url = 'http://epileptologie-bonn.de/cms/upload/workgroup/lehnertz/'
    urls_suffix = ['Z.zip', 'O.zip', 'N.zip', 'F.zip', 'S.zip']

    path_child_fold = zip_with_unique(path_data, child_fold)

    if fold.exists():
        print("Folder already exists")
        check_child_folders = [Path(child).exists()
                               for child in path_child_fold]

        if all(check_child_folders):
            print("Subfolders already exist")

            return path_child_fold
    else:
        print("Creating folder")
        # Create parent directory
        fold.mkdir(parents=True, exist_ok=True)
        # This way, the child directory will also be created.
        for child in path_child_fold:
            Path(child).mkdir(parents=True, exist_ok=True)

        urls = zip_with_unique(base_url, urls_suffix)

        print("Downloading and unzipping the files")

        for url, path in list(zip(urls, path_child_fold)):
            file_directory = download(url, path)

            with ZipFile(file_directory, "r") as zip_ref:
                zip_ref.extractall(path)

    return path_child_fold


def read_boon(path_child_fold) -> array:
    """Function for reading the boon database, and return X and y.
    Also adapted from:
    https://mne-tools.github.io/mne-features/auto_examples/plot_seizure_example.html
    Parameters
    ----------

    path_child_fold : TO-DO

    Returns
    -------
    X : array-like, shape (n_samples, n_features)
        Data vectors, where n_samples is the number of samples
        and n_features is the number of features.
    y : array-like, shape (n_samples,)
        Target values.

    """

    data_segments = list()
    labels = list()

    for path in path_child_fold:

        f_names = [s for s in Path(path).iterdir() if str(s).lower().endswith('.txt')]

        for f_name in f_names:
            _data = read_csv(f_name, sep='\n', header=None)

            data_segments.append(_data.values.T[None, ...])

        if ('setE' in path) or ('setC' in path) or ('setD' in path):

            labels.append(ones((len(f_names),)))
        else:
            labels.append(zeros((len(f_names),)))

    X = concatenate(data_segments).squeeze()
    y = concatenate(labels, axis=0)

    return X, y


def preprocessing_split(X, y, test_size=.20, random_state=42):
    """Function to perform the train and test split 
    and normalize the data set with Min-Max.
    
    Parameters
    ----------
        
    X : array-like, shape (n_samples, n_features)
        Training vectors, where n_samples is the number of samples
        and n_features is the number of features.
        
    y : array-like, shape (n_samples,)
        Target values.
    
    test_size : float
        value between 0 and 1 to indicate the 
        percentage that will be used in the test.
    
    random_state : int
        seed to be able to replicate split
        
    Returns
    -------
    
    TO-DO: Explanation that will be 
    the separation between training and testing.


    """

    X_train, X_test, Y_train, Y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state)

    # MinMax Scaler

    minMax = MinMaxScaler()
    minMax = minMax.fit(X_train)

    X_train = minMax.transform(X_train)
    X_test = minMax.transform(X_test)

    X_train = X_train[:, :4096]
    X_test = X_test[:, :4096]

    return X_train, X_test, Y_train, Y_test
