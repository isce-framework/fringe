import DesignMatrix
import numpy as np

# instantiate the design matrix class



pairs=['1_2','1_3','1_4', '2_3', '2_4', '2_5', 
          '3_4', '3_5', '3_6','4_5', '4_6', '4_7', 
          '5_6', '5_7', '5_8', '6_7', '6_8', '6_9',
          '7_8', '7_9', '8_9']

# a linear signal
signal = np.arange(9)
signal2 = np.arange(9) + 1
signal = signal[1:] 
signal2 = signal2[1:]

# get a design matrix object
dm = DesignMatrix.DesignMatrix()
dm.pairs = pairs
dm.configure()

# get the SBAS design matrix
dm.timeseries()
A = dm.getG()

# simulate interferometric phases
pairs_phase = np.dot(A,signal)
pairs_phase2 = np.dot(A,signal2)
#print(pairs_phase)

# simulate unwrapping error for one pair
unw_err = np.zeros_like(pairs_phase)
unw_err[5] = 2*np.pi

# add unw error to the interferometric phases
pairs_phase_err = pairs_phase + unw_err
pairs_phase2_err = pairs_phase2 + unw_err

# get the phase closure design matrix
dm.closure()
G = dm.getG()

#compute closure
closure = np.dot(G, pairs_phase_err)
closure2 = np.dot(G, pairs_phase2_err)

#################
# unw err estimation based on closure
C = -2.0*np.pi*G 

X1 = np.dot(np.linalg.pinv(C), closure)
U1 = np.round(X1)

X2 = np.dot(np.linalg.pinv(C), closure2)
U2 = np.round(X2)

pairs_phase_corrected = pairs_phase_err + 2.0*np.pi*U1
pairs_phase2_corrected = pairs_phase2_err + 2.0*np.pi*U2
###############
#residual of the correction and simulated pairs phases before adding unw err
residual_before = pairs_phase_err - pairs_phase
residual_after = pairs_phase_corrected - pairs_phase
residual_after2 = pairs_phase2_corrected - pairs_phase2
closure_after = np.dot(G, pairs_phase_corrected)

print("----------------------")
print("residual after correction ")
print(residual_after)
print(residual_after2)

#######################
numPixels = 2
CC = C.copy()
for ii in range(1,numPixels):
    CC = np.vstack((CC,C))


L = np.hstack((closure, closure2))

X_region = np.dot(np.linalg.pinv(CC), L)
U_region = np.round(X_region)

print("----------------------")
print("estimated ambiguity for the region:")
print("U_region", U_region)
pairs_phase_corrected = pairs_phase_err + 2.0*np.pi*U_region
pairs_phase2_corrected = pairs_phase2_err + 2.0*np.pi*U_region

residual_after = pairs_phase_corrected - pairs_phase
residual_after2 = pairs_phase2_corrected - pairs_phase2

print("----------------------")
print("residual after correction region")
print(residual_after)
print(residual_after2)


