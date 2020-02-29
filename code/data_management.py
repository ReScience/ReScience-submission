from pathlib import Path
from wget import download
from zipfile import ZipFile
from pandas import read_csv
from numpy import zeros, ones, concatenate, array


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
        print(check_child_folders)

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
    
    Parameters
    ----------

    path_child_fold : TO-DO
    
    Returns
    -------
    X : array-like, shape (n_samples, n_features)
        Training vectors, where n_samples is the number of samples
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

    X = concatenate(data_segments)
    y = concatenate(labels, axis=0)

    return X, y

