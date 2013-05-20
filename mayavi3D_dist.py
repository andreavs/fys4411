import sys, os, numpy as np

from math import sqrt
from mayavi import mlab
import matplotlib.pyplot as plt
#from pylab import *


def parseCML():
    if len(sys.argv[1]) < 2:
        raise Exception("Path to 3D file must be given as first cml arg.")
    else:
        path = sys.argv[1]
        
        # if not os.path.exists(path):
        #     raise Exception("%s does not exist on your file system." % path)
        
    return path        
        
def loadArmaCube(path):
    """reading a armadillo binary cube representing a 3D histogram"""
    
    with open(path, 'rb') as binFile:
        binFile.readline()
        n, m, l = [int(n_i) for n_i in binFile.readline().strip().split()]
        
        data = np.fromfile(binFile, dtype=np.float64)
        
    print ("Loaded %d data points" % data.shape),
    
    data.resize(n, m, l)            
    
    print "reshaped to a %s array. " % str(data.shape)
    
    return data

def loadArmaMat(path):
    with open(path, 'rb') as binFile:
        binFile.readline()
        n, m = [int(n_i) for n_i in binFile.readline().strip().split()]
        
        data = np.fromfile(binFile, dtype=np.float64)
        
    print ("Loaded %d data points" % data.shape),
    
    data.resize(n, m)            
    
    print "reshaped to a %s array. " % str(data.shape)
    
    return data

def loadArmaVec(path):
    with open(path, 'rb') as binFile:
        binFile.readline()
        n, m = [int(n_i) for n_i in binFile.readline().strip().split()]
        
        data = np.fromfile(binFile, dtype=np.float64)
        
    print ("Loaded %d data points" % data.shape),
    
    data.resize(n, m)            
    
    print "reshaped to a %s array. " % str(data.shape)
    
    return data

def earthSpherify(data):
    """Creates a spherical representation of the data with a slice to the center"""
    
    n, m, l = data.shape    
    
    f = lambda i, n: (i*2.0/(n-1) - 1)**2
    f = np.vectorize(f)

    D_i = f(xrange(0, n), n)
    D_j = f(xrange(0, m), m)
    D_k = f(xrange(0, l), l)
    nBins = 100
    rhist = np.linspace(0,1,nBins)
    hist = np.zeros(nBins)
    #Create the sphere
    for i , d_i in enumerate(D_i):
        for j, d_j in enumerate(D_j):
            for k, d_k in enumerate(D_k):
                r = sqrt(d_i + d_j + d_k)
                if r > 2:
                    data[i, j, k] = 0
                else:
                    pos = int(r*nBins)
                    hist[pos]+=data[i,j,k];
    #Create the slice
    plt.plot(rhist,hist)
    plt.show()
    data[n/2:, m/2:, l/2:] = 0;
    
    return data;

def plot1d(data):
    plt.plot(data);
    plt.show()

def plot2d(data):
    plt.imshow(data)
    plt.show()
    


def main():
    
    path = parseCML()
    path = "VMC-build-desktop-Qt_4_8_1_in_PATH__System__Release/results/oneBody/" + path + "res"

    if sys.argv[2] == "1":
        path = path + "1d" + sys.argv[3] + ".mat"
        data = loadArmaVec(path)
        plot1d(data)

    elif sys.argv[2] == "2":
        path = path + "2d" + sys.argv[3] + ".mat"
        data = loadArmaVec(path)
        plot2d(data)
    
    elif sys.argv[2] == "3":
        path = path + "3d.mat"
        data = earthSpherify(loadArmaCube(path))
        print np.max(data)
        mlab.pipeline.volume(mlab.pipeline.scalar_field(data), vmin=0, vmax=1)
        mlab.show()
        #data = loadArmaCube(path)


    
if __name__ == "__main__":
    main()
