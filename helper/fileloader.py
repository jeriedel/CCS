import pandas as pd
import numpy as np
import h5py
import re
import os

class FileLoader:
    floatRgx = re.compile(r'[+-]*[0-9]+\.*[0-9]*[e,E]*[+-]*[0-9]*')
    chrRgx   = re.compile(r"[^eE\-\+\d\.\s\t]")

    @classmethod
    def getFloats(cls, iterable):
        data = []
        for line in iterable:
            if re.search(FileLoader.chrRgx, line):
                pass
            else:
                data.append(re.findall(FileLoader.floatRgx, line))
        return data
            
    @classmethod
    def cleanContent(cls, data):
        cleaned = []
        for line in data:
            if len(line) < 2:
                pass
            else:
                clean = [float(x) for x in line]
                cleaned.append(clean)
        return cleaned
        
    @classmethod
    def openFile(cls, p2f):
        p2f = os.path.expanduser(p2f)
        with open(p2f, 'r') as f:
            content = f.readlines()
        return content
    
    @classmethod
    def generalInterface(cls, p2f):
        content = FileLoader.openFile(p2f)
        data    = FileLoader.getFloats(content)
        clean   = FileLoader.cleanContent(data)
        return clean
        
    @staticmethod
    def getFileAsArray(p2f, cols=None):
        data = FileLoader.generalInterface(p2f)
        data = np.asarray(data)
        if cols is not None:
            try:
                data = data[:, cols]
            except Exception:
                return np.array()
        return data
    
    @staticmethod
    def getFileAsDataFrame(p2f, cols=None):
        data = FileLoader.generalInterface(p2f)
        data = pd.DataFrame(data=data)
        
        if cols is not None:
            try:
                data = data.iloc[:, cols]
                print(data)
            except Exception:
                return pd.DataFrame()
            
        return data

class ImobFileLoader(FileLoader):
    @staticmethod
    def getFileAsDataFrame(p2f):
        pass

