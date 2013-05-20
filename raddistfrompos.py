import sys, os, numpy as np

from math import sqrt
from mayavi import mlab
import matplotlib.pyplot as plt
#from pylab import *

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

def makeRadialHist1(data):
    n, m, l = data.shape    
    #print data
    #f = lambda i, n: (i*2.0/(n-1) - 1)**2
    #f = np.vectorize(f)

    #D_i = f(xrange(0, n), n)
    #D_j = f(xrange(0, m), m)
    #D_k = f(xrange(0, l), l)
    nBins = 100
    radii = np.zeros((l,n));
    print l,m,n
    hist = np.zeros(nBins)
    for i in xrange(l):
    	#if(i%100 == 0):
    		#print i
    	for j in xrange(n):
    		for k in xrange(m):
    			radii[i,j] += data[j,k,i]**2;
    radii = np.sqrt(radii)
    
    a = np.histogram(radii,1000)
    print a
    b = a[1][-1]
    a = a[0]

    t = np.linspace(0,b,1000)
    plt.hist(radii,1000)
    plt.figure()
    plt.plot(t, a)
    plt.xlabel('radial distance')
    plt.ylabel('non-normalized distribution')
    plt.show()

def make1dHist(data):
    t = np.linspace(-5,5,200)
    plt.plot(t, data)
    plt.xlabel('x-position')
    plt.ylabel('non-normalized density')
    plt.show()

def make2dHist(data):
    
    plt.imshow(data, extent=[-5,5,-5,5])
    plt.xlabel('x-axis')
    plt.ylabel('y-axis')
    plt.colorbar()
    plt.show()

def make3dHist(data):
    # n, m, l = data.shape
    # results = np.zeros((n*l,3))
    # for i in range(n):
    #     for j in range(m):
    #         results[i*l:(i+1)*l,j] = data[i,j,:]

    # H, edges = np.histogramdd(results, bins = (100,100,100))
    # print H.shape
    data = earthSpherify(data)
    #mlab.contour3d(data)
    data = data/np.max(data)*3.74
    mlab.pipeline.volume(mlab.pipeline.scalar_field(data), vmin=0, vmax=1.0)
    mlab.show()

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
                if r > 1:
                    data[i, j, k] = 0
                else:
                    pos = int(r*nBins)
                    #hist[pos]+=data[i,j,k];
    #Create the slice
    # plt.plot(rhist,hist)
    # plt.show()
    data[n/2:, m/2:, :] = 0
    
    return data;

def makeRadialHist2(data):
    n, m = data.shape    
    # nBins = 100
    # radii = np.zeros(n*m);
    # print l,m,n
    # hist = np.zeros(nBins)	
    # print sum(sum(np.isnan(data)))
    a = np.histogram(data,100,normed=True)
    b = a[1][-1]
    c = a[1][0]
    a = a[0]
    t = np.linspace(c,b,100)
    plt.plot(t, a)
    plt.xlabel('x-position')
    plt.ylabel('non-normalized density')
    #plt.show()

def loadArmaVec(path):
    with open(path, 'rb') as binFile:
        binFile.readline()
        n, m = [int(n_i) for n_i in binFile.readline().strip().split()]
        
        data = np.fromfile(binFile, dtype=np.float64)
        
    print ("Loaded %d data points" % data.shape),
    
    data.resize(n, m)            
    
    print "reshaped to a %s array. " % str(data.shape)
    
    return data

if __name__ == "__main__":
    path1 = "VMC-build-desktop-Qt_4_8_1_in_PATH__System__Release/results/" + sys.argv[1] + "radii.mat"
    path2 = "VMC-build-desktop-Qt_4_8_1_in_PATH__System__Release/results/" + sys.argv[1] + "positions.mat"
    path3 = "VMC-build-desktop-Qt_4_8_1_in_PATH__System__Release/results/" + sys.argv[1] + "x.mat"
    path4 = "VMC-build-desktop-Qt_4_8_1_in_PATH__System__Release/results/" + sys.argv[1] + "xy.mat"
    path5 = "VMC-build-desktop-Qt_4_8_1_in_PATH__System__Release/results/" + sys.argv[1] + "xyz.mat"
    path6 = "VMC-build-desktop-Qt_4_8_1_in_PATH__System__Release/results/" + sys.argv[1] + "radii1.mat"
    path7 = "VMC-build-desktop-Qt_4_8_1_in_PATH__System__Release/results/" + sys.argv[1] + "radii2.mat"
    path8 = "VMC-build-desktop-Qt_4_8_1_in_PATH__System__Release/results/" + sys.argv[1] + "radii3.mat"
    path9 = "VMC-build-desktop-Qt_4_8_1_in_PATH__System__Release/results/" + sys.argv[1] + "radii4.mat"
    #make2dHist(loadArmaCube(path2), 0)
    #makeRadialHist1(loadArmaCube(path2))
    makeRadialHist2(loadArmaMat(path6))
    plt.hold('on')
    makeRadialHist2(loadArmaMat(path7))
    makeRadialHist2(loadArmaMat(path8))
    makeRadialHist2(loadArmaMat(path9))
    # data = loadArmaVec(path3)
    # make1dHist(data);
    # data = loadArmaMat(path4)
    # make2dHist(data);
    # data = loadArmaCube(path5)
    # make3dHist(data);
    plt.xlabel('Distance [a.u]')
    plt.ylabel('Probability')
    plt.legend(['with correlations', 'without correlations, alpha=10', 'without correlations, alpha=7.81', 'with correlations, alpha=10, beta=0.104'])
    plt.show()