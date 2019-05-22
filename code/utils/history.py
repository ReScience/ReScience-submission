import pickle as pkl
from os.path import isfile

class History:
    def __init__(self, file):
        self._file = file
        if not isfile(file):
            self._hist = {}

        else:
            with open(file, 'rb') as f:
                self._hist = pkl.load(f)

    def push(self, key, timestamp, value):
        if not key in self._hist:
            self._hist[key] = {}
        self._hist[key][timestamp] = value


    def dump(self):
        with open(self._file, 'wb') as f:
            pkl.dump(self._hist, f)

    def cut(self, timestamp):
        for key, h in self._hist.items():
            t_to_remove = [t for t in h.keys() if t >= timestamp]
            for t in t_to_remove:
                h.pop(t)
            