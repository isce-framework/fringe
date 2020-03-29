#!/usr/bin/env python3

import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt


Baseline = np.array([   0.        ,  202.23261006 , 247.67336194 , 157.46570522 ,  -3.03759106,
  146.51305542,  432.34588205 , 129.9950411  , 304.6576168  , 342.39451334,
   77.48752496,   35.375345   ,  49.16112925 ,-291.19738937 , 525.84622291,
  131.16981696,  -73.74287393 ,-267.98532137 , 126.57065401 ,  51.46784891,
 -137.74090081,  204.51096416 , 519.05312568 ,  75.56740286 ,-253.08744983,
  321.84328088,  372.77095751 , 322.70861019 , 239.31234022 , -81.07387604,
  132.89024919,  -71.70800096 , 210.09229394 , -47.12478361 , 312.45312082,
  152.03018981,  366.71227489 , 310.07070069 ,-216.6017471  , -85.26060589])


def grid_search(cpxPhase, K, DEM_error_range):

    gamma = np.zeros_like(DEM_error_range, dtype=np.float64)
    for ii, dz in enumerate(DEM_error_range):
        pred_cpxphase = np.exp(1J*K*dz)
        dif = cpxPhase * np.conjugate(pred_cpxphase)
        tempCoh = np.abs(np.sum(dif)/len(dif))
        gamma[ii] = tempCoh
        #print(ii, tempCoh)

    maxIdx = np.argmax(gamma)
    #print(maxIdx, gamma[maxIdx])
    #print("##############################")
    return DEM_error_range[maxIdx]

def grid_search_2D(cpxPhase, Kb, Kt, DEM_error_range, rate_range):

    gamma = np.zeros((len(DEM_error_range), len(rate_range)), dtype=np.float64)
    for ii, dz in enumerate(DEM_error_range):
        for jj,rate in enumerate(rate_range):
            pred_cpxphase = np.exp(1J*(Kb*dz+Kt*rate))
            dif = cpxPhase * np.conjugate(pred_cpxphase)
            tempCoh = np.abs(np.sum(dif)/len(dif))
            gamma[ii,jj] = tempCoh
    
    maxIdx = np.unravel_index(np.argmax(gamma, axis=None), gamma.shape)
    dz = DEM_error_range[maxIdx[0]]
    rate = rate_range[maxIdx[1]]
    #plt.imshow(gamma)
    #plt.show()
    #print("##############################")
    return dz, rate



def estimate_dem_error(cpxPhase, wavelength, bperp, rng, theta, dt):

    PI = np.pi
    K = 4.0 * PI * bperp / (wavelength * rng * np.sin(theta))
    Kt = 4.0 * PI * dt/wavelength

    DEM_error_range = np.arange(-200,200,2)
    rate_range = np.arange(-0.3,0.3,0.01)
    dz, rate = grid_search_2D(cpxPhase, K, Kt, DEM_error_range, rate_range)

    #dz = grid_search(cpxPhase, K, DEM_error_range)
    
    DEM_error_range = np.arange(dz-10,dz+10,1) 
    rate_range = np.arange(rate-0.01,rate+0.01,0.002)
    dz, rate = grid_search_2D(cpxPhase, K, Kt, DEM_error_range, rate_range)


    DEM_error_range = np.arange(dz-2,dz+2,0.1)
    rate_range = np.arange(rate-0.005,rate+0.005,0.001)
    dz, rate = grid_search_2D(cpxPhase, K, Kt, DEM_error_range, rate_range)

    #DEM_error_range = np.arange(dz-1,dz+1,.1)
    #dz = grid_search(cpxPhase, K, DEM_error_range)

    return dz, rate


class stack(object):
    def __init__(self, ySize=100, xSize=100, numDates=40):
        self.numDates= numDates
        self.ySize = ySize
        self.xSize = xSize
        self.wavelength = 0.056 
        self.r = 850000 # avergae range
        self.theta = 34.0 * np.pi/180.0

    def configure(self, dt = None, dB= None):
        self.stackSize = (self.ySize, self.xSize, self.numDates)
        if dt is None:
            dt = np.arange(self.numDates)*12.0/365.0

        if dB is None:
            dB = np.random.randn(self.numDates)*200

        self.dB = dB - dB[0]
        self.dt = dt
        self.signal = np.zeros(self.stackSize)
        self.dz = np.zeros((self.ySize, self.xSize))
        self.rate = np.zeros((self.ySize, self.xSize))
        self.coherence = np.ones((self.ySize, self.xSize))
        self.coherenceThreshold = 1.0
    def add_linear_signal(self, rate = None): #, dz=None):
    
        if rate is None:
            xCenter = int(self.xSize/2)
            yCenter = int(self.ySize/2)
            x,y = np.meshgrid(np.arange(self.xSize),np.arange(self.ySize))
            r = np.sqrt((x-xCenter)**2 + (y-yCenter)**2)
            r = (np.max(r) - r)

            radius = np.min((xCenter, yCenter))
            ind = r>radius
            

            rate = r.copy()
            
            min_rate = np.min(rate[ind])
            max_rate = np.max(rate[ind])
            rate = (rate-min_rate)/(max_rate-min_rate)

            rate[r<radius] = 0
            rate = rate * 0.1

        self.rate = rate

        for i in range(self.numDates):
            self.signal[:,:,i] = self.signal[:,:,i] + (4*np.pi/self.wavelength)*self.rate*self.dt[i]

    def add_dem_error(self, max_demError=50.0, dz=None):

        if dz is None:
            x,y = np.meshgrid(np.arange(self.xSize),np.arange(self.ySize))
            dz = 100*(np.sin(2*2*np.pi*x/self.xSize) + np.cos(2*2*np.pi*x/self.xSize))#+ np.sin(2*np.pi*x/2) + np.sin(2*np.pi*x/6) + np.sin(2*np.pi*x/10)

        dz = (dz - np.min(dz))/(np.max(dz)-np.min(dz))

        self.dz = max_demError*dz

        demErrorCoef = 4.0*np.pi/(self.wavelength*self.r*np.sin(self.theta))
        #self.signal = np.zeros(self.stackSize)
        for i in range(self.numDates):
            self.signal[:,:,i] = self.signal[:,:,i] + demErrorCoef*self.dB[i]*self.dz

    def add_nonlinear_deformation(self):
        for i in range(self.numDates):
            acceleration = self.rate/4.0
            self.signal[:,:,i] = self.signal[:,:,i] + (4*np.pi/self.wavelength)*acceleration*(self.dt[i]**2)

    def add_atmosphere(self, max_atmosphere):
        x,y = np.meshgrid(np.arange(self.xSize),np.arange(self.ySize))
        randomCoef_X = np.random.randn(self.numDates) 
        randomCoef_Y = np.random.randn(self.numDates)
        for i in range(1,self.numDates):
            ramp = randomCoef_X[i]*x + randomCoef_Y[i]*y
            ramp = max_atmosphere*(ramp-np.min(ramp))/(np.max(ramp)-np.min(ramp))
            self.signal[:,:,i] = self.signal[:,:,i] + ramp

    def simulate_coherence(self):
        coh = np.random.rand(self.ySize, self.xSize)
        coh = (coh-np.min(coh))/(np.max(coh) - np.min(coh))

        self.coherence = coh

    def wrapped_phase(self, coherenceThreshold = 0.3):
        self.coherenceThreshold = coherenceThreshold
        self.phase = np.zeros(self.stackSize)
        for i in range(self.numDates):
            ph = st.signal[:,:,i]
            ph[st.coherence<coherenceThreshold] = np.nan
            self.phase[:,:,i] = np.angle(np.exp(1j*ph))


    def DelaunayTriangulation(self):
        ind = self.coherence>self.coherenceThreshold
        x,y = np.meshgrid(np.arange(self.xSize),np.arange(self.ySize))
        # coherent points
        points = np.vstack((x[ind],y[ind])).T
        self.tri = Delaunay(points)
        
        #self.tri = tri

        plt.imshow(self.phase[:,:,-1])
        plt.triplot(points[:,0], points[:,1], self.tri.simplices.copy())
        plt.plot(points[:,0], points[:,1], 'o')
        plt.show()

        self.tri_coords = points[self.tri.simplices]

        tt = self.tri_coords[-1]
        arc01 = self.phase[tt[0][1], tt[0][0],:] - self.phase[tt[1][1], tt[1][0],:]
        arc01_c = np.exp(1J*self.phase[tt[0][1], tt[0][0],:]) * np.exp(-1J*self.phase[tt[1][1], tt[1][0],:])

        arc_unw = self.signal[tt[0][1], tt[0][0],:] - self.signal[tt[1][1], tt[1][0],:]

        A=np.ones((self.numDates,2))
        A[:,0] = (4*np.pi/self.wavelength)*self.dt
        A[:,1] = 4.0*np.pi*self.dB/(self.wavelength*self.r*np.sin(self.theta))
        X = np.dot(np.linalg.pinv(A), arc01)
        X_c = np.dot(np.linalg.pinv(A), np.angle(arc01_c))

        #predicted_L = np.dot(A,X_c)
        #dz_gridSearch = coarse_grid_search(arc01_c, self.wavelength, self.dB, self.r, self.theta) 
        dz_gridSearch, rate_gridSearch = estimate_dem_error(arc01_c, self.wavelength, self.dB, self.r, self.theta, self.dt)
        print("%%%%%%%%%%%%%%%%%%%%%%%")
        print("simulated values (difference on the arc)")
        print("rate: " , self.rate[tt[0][1], tt[0][0]] - self.rate[tt[1][1], tt[1][0]])
        print("dZ: " , self.dz[tt[0][1], tt[0][0]] - self.dz[tt[1][1], tt[1][0]])
        print("%%%%%%%%%%%%%%%%%%%%%%%")
        '''print("estimated values")
        print("observation = wrapped phase diff")
        print("rate", X[0])
        print("dz", X[1])
        print("    ***************      ")
        print("observation = angle( (exp(1J*Phi) * exp(-1J*Phi) )")
        print("rate", X_c[0])
        print("dz", X_c[1])
        print("    ***************      ")
        '''
        print("grid search")
        print("rate: ", rate_gridSearch)
        print("dz: ", dz_gridSearch )
        print()
        print("%%%%%%%%%%%%%%%%%%%%%%%")
        #arc12 = data[tt[1][1], tt[1][0]] - data[tt[2][1], tt[2][0]]
        #arc20 = data[tt[2][1], tt[2][0]] - data[tt[0][1], tt[0][0]]
        #print(tt)
        #print(arc01)
        #plt.plot(arc01)
        #plt.plot(np.angle(arc01_c))
        #plt.plot(arc_unw, '--')
        #plt.plot((arc_unw-predicted_L), '-^')
        #plt.show()

    def estimate_arcs(self):
        #coordinates of all triangle vertices
        #tri_coords = points[tri.simplices]

        # A design matrix for Least Square estimation
        #A = np.ones((self.numDates,2))
        #A[:,0] = (4*np.pi/self.wavelength)*self.dt
        #A[:,1] = 4.0*np.pi*self.dB/(self.wavelength*self.r*np.sin(self.theta))
        #A1 = np.linalg.pinv(A)

        #startVertex = []
        #endVertex = []
        arc = []
        rateDiff = []
        dzDiff = [] 
        st = [0,1,2]
        en = [1,2,0]

        for ii in range(self.tri.nsimplex):
            print(ii, " of ", self.tri.nsimplex)
            tt = self.tri_coords[ii]
            for jj in range(3):
                ss = st[jj]
                ee = en[jj]
                #startVertex.append((tt[ss][1], tt[ss][0]))
                #endVertex.append((tt[ee][1], tt[ee][0]))
                arcIdx = (tt[ss][1], tt[ss][0], tt[ee][1], tt[ee][0])
                if arcIdx not in arc:
                    arc.append(arcIdx)
                    arc_phase = np.exp(1J*self.phase[tt[ss][1], tt[ss][0],:]) * np.exp(-1J*self.phase[tt[ee][1], tt[ee][0],:])
                    dz_gridSearch, rate_gridSearch = estimate_dem_error(arc_phase, self.wavelength, self.dB, self.r, self.theta, self.dt)
                    rateDiff.append(rate_gridSearch)
                    dzDiff.append(dz_gridSearch)
            #startVertex.append((tt[0][1], tt[0][0]))
            #endVertex.append((tt[1][1], tt[1][0]))
            #arc01 = np.exp(1J*self.phase[tt[0][1], tt[0][0],:]) * np.exp(-1J*self.phase[tt[1][1], tt[1][0],:])
            
            #X = np.dot(A1, np.angle(arc01))
            #rateDiff.append(X[0])
            #dzDiff.append(X[1])

        return arc, rateDiff, dzDiff

    def integrate_arcs_sparse(self, arc, rateDiff, dzDiff):
        rows = np.zeros((1,2*len(arc)))
        cols = np.zeros((1,2*len(arc)))
        data = np.zeros((1,2*len(arc)))
        #A = np.zeros((len(arc), self.xSize*self.ySize))
        for ii,aa in enumerate(arc):
            Yst = aa[0]
            Xst = aa[1]
            Yend = aa[2]
            Xend = aa[3]
            st = Yst*self.xSize + Xst
            end = Yend*self.xSize + Xend
            print("*******")
            print(aa)
            print(st, end)
            A[ii, st] = -1.0
            A[ii, end] = 1.0
            cols[0,2*ii] = ii
            rows[0,2*ii] = st

            data[0,2*ii] = -1

            cols[0,2*ii+1] = ii
            rows[0,2*ii+1] = end
            data[0,2*ii+1] = 1

        A = coo_matrix((data, (row, col)), shape=((len(arc), self.xSize*self.ySize))) 
        

    def integrate_arcs(self, arc, rateDiff, dzDiff):
        A = np.zeros((len(arc), self.xSize*self.ySize))
        for ii,aa in enumerate(arc):
            Yst = aa[0]
            Xst = aa[1]
            Yend = aa[2]
            Xend = aa[3]
            st = Yst*self.xSize + Xst
            end = Yend*self.xSize + Xend
            print("*******")
            print(aa)
            print(st, end)
            A[ii, st] = -1.0
            A[ii, end] = 1.0



        Yref = 0
        Xref = 0
        refPixel = Yref*self.xSize + Xref
        A = np.hstack((A[:,0:refPixel], A[:,(refPixel+1):]))

        rankA = np.linalg.matrix_rank(A)
        print("size of design mat: ", A.shape)
        print("rank of A: ", rankA)
        rate = np.dot(np.linalg.pinv(A), np.array(rateDiff))
        dz = np.dot(np.linalg.pinv(A), np.array(dzDiff))

        dz = np.hstack((0,dz))
        rate = np.hstack((0,rate))

        rate = np.reshape(rate, (self.ySize, self.xSize))
        dz = np.reshape(dz, (self.ySize, self.xSize))

        return rate , dz

    def integrate_arcs_coh(self, arc, rateDiff, dzDiff):
        A = np.zeros((len(arc), self.xSize*self.ySize))
        pixelsIdx = [False] * self.xSize*self.ySize 
        for ii,aa in enumerate(arc):
            Yst = aa[0]
            Xst = aa[1]
            Yend = aa[2]
            Xend = aa[3]
            st = Yst*self.xSize + Xst
            end = Yend*self.xSize + Xend
            print("*******")
            print(aa)
            print(st, end)
            A[ii, st] = -1.0
            A[ii, end] = 1.0
            pixelsIdx[st] = True
            pixelsIdx[end] = True
                    
        A = A[:, pixelsIdx]
        
        Yref = 0
        Xref = 0
        refPixel = Yref*self.xSize + Xref
        
        #refPixel = end
        #Yref = int(end/self.xSize)
        #Xref = end-Yref*self.xSize
        A = np.hstack((A[:,0:refPixel], A[:,(refPixel+1):]))

        rankA = np.linalg.matrix_rank(A) 
        print("size of design mat: ", A.shape)
        print("rank of A: ", rankA)
        rate = np.dot(np.linalg.pinv(A), np.array(rateDiff))
        dz = np.dot(np.linalg.pinv(A), np.array(dzDiff))
        
        dz = np.hstack((0,dz))
        rate = np.hstack((0,rate))

        dz_all = np.zeros((self.xSize*self.ySize))
        rate_all = np.zeros((self.xSize*self.ySize))

        dz_all[pixelsIdx] = dz
        rate_all[pixelsIdx] = rate

        rate_all = np.reshape(rate_all, (self.ySize, self.xSize))
        dz_all = np.reshape(dz_all, (self.ySize, self.xSize))

        rate_all[self.coherence<self.coherenceThreshold] = np.nan
        dz_all[self.coherence<self.coherenceThreshold] = np.nan

        refPixel = end
        Yref = int(end/self.xSize)
        Xref = end-Yref*self.xSize

        rate_all = rate_all - rate_all[Yref, Xref]
        dz_all = dz_all - dz_all[Yref, Xref]

        return rate_all , dz_all, Xref, Yref
#####################################

st = stack()
st.xSize = 256
st.ySize = 256
st.configure(dB=Baseline)

st.add_linear_signal()
st.add_dem_error(max_demError=50.0)
st.add_atmosphere(max_atmosphere=10) # maximum atmosphere over the scene in radians
st.add_nonlinear_deformation()
st.simulate_coherence()
st.wrapped_phase(coherenceThreshold = 0.3)
st.DelaunayTriangulation()
arc, rateDiff, dzDiff = st.estimate_arcs()
#rate , dz = st.integrate_arcs(arc, rateDiff, dzDiff)
rate , dz, Xref, Yref = st.integrate_arcs_coh(arc, rateDiff, dzDiff)

fig = plt.figure(figsize=(8,6))
ax1 = fig.add_subplot(2,3,1)
cax1=ax1.imshow(st.dz-st.dz[Yref, Xref])
ax1.set_title("simulated DEM error")
cbar = fig.colorbar(cax1)

ax2 = fig.add_subplot(2,3,2)
cax2=ax2.imshow(-1*dz)
ax2.set_title("estimated DEM error")
cbar = fig.colorbar(cax2)

ax3 = fig.add_subplot(2,3,3)
cax3=ax3.imshow(st.dz -st.dz[Yref,Xref] + dz, vmin = -1, vmax = 1)
ax3.set_title("difference")
cbar = fig.colorbar(cax3)

ax4 = fig.add_subplot(2,3,4)
cax4=ax4.imshow(st.rate-st.rate[Yref,Xref])
cbar = fig.colorbar(cax4)
ax4.set_title("simulated rate")

ax5 = fig.add_subplot(2,3,5)
cax5=ax5.imshow(-1*rate)
cbar = fig.colorbar(cax5)
ax5.set_title("estimated rate")

ax6 = fig.add_subplot(2,3,6)
cax6=ax6.imshow(st.rate -st.rate[Yref,Xref] + rate, vmin = -0.01, vmax = 0.01)
cbar = fig.colorbar(cax6)
ax6.set_title("difference")


plt.savefig("simulation.png")

plt.show()

#print(st.stackSize)
#print(st.dt)
#print(st.dB)

#import matplotlib.pyplot as plt
'''plt.imshow(st.displacement[:,:,-1])
plt.imshow(np.angle(np.exp(1j*st.phase[:,:,-1])))
plt.show()
plt.plot(st.phase[125,125,:])
plt.show()

plt.imshow(st.coherence, cmap='gray')
plt.show()
'''

#ph = st.phase[:,:,-1]
#ph[st.coherence<0.3] = np.nan
#ph = np.angle(np.exp(1j*ph))
#plt.imshow(st.phase[:,:,-1])
#plt.show()

#plt.plot(st.signal[int(st.ySize/2), int(st.xSize/2),:])
#plt.show()

# Simulate displacement 

# Simulate DEM error

# Simulate wrapped phase

# Simulate average coherence (temporal coherence)

# Delaunay triangulation of coherent pixels

# Compue edge phases

# estimate v and dz for edges


