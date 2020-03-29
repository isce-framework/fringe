import DesignMatrix
import numpy as np

# instantiate the design matrix class



pairs=['1_2','1_3','1_4', '2_3', '2_4', '2_5', 
          '3_4', '3_5', '3_6','4_5', '4_6', '4_7', 
          '5_6', '5_7', '5_8', '6_7', '6_8', '6_9',
          '7_8', '7_9', '8_9']

# a linear signal
signal = np.arange(9)
signal = signal[1:] 

# get a design matrix object
dm = DesignMatrix.DesignMatrix()
dm.pairs = pairs
dm.configure()

# get the SBAS design matrix
dm.timeseries()
A = dm.getG()

# simulate interferometric phases
pairs_phase = np.dot(A,signal)
print(pairs_phase)

# simulate unwrapping error for one pair
unw_err = np.zeros_like(pairs_phase)
unw_err[5] = 2*np.pi

# add unw error to the interferometric phases
pairs_phase_err = pairs_phase + unw_err

# get the phase closure design matrix
dm.closure()
G = dm.getG()

#compute closure
closure = np.dot(G, pairs_phase_err)

#################
# unw err estimation based on closure
C = -2.0*np.pi*G 

X = np.dot(np.linalg.pinv(C), closure)
U = np.round(X)

pairs_phase_corrected = pairs_phase_err + 2.0*np.pi*U

###############
#residual of the correction and simulated pairs phases before adding unw err
residual_before = pairs_phase_err - pairs_phase
residual_after = pairs_phase_corrected - pairs_phase
closure_after = np.dot(G, pairs_phase_corrected)

print("----------------------")
print("residual before correction ")
print(residual_before)
print("residual after correction ")
print(residual_after)

print("----------------------")
print("closure before correction ")
print(closure)
print("closure after correction ")
print(closure_after)


np.testing.assert_almost_equal(residual_after, np.zeros(residual_after.shape))

