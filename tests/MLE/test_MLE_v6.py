#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
import argparse
from mle import estimate_MLE_4 as mle4

def createParser():
    '''
    Create command line parser.
    '''

    parser = argparse.ArgumentParser( description='A simulation for time-series analysis')
    parser.add_argument('-t', '--time_length', type=int, dest='timeLength', required=True,
            help='Duration of time-series in Days.')
    parser.add_argument('-i', '--acquisition_interval', type=int, dest='dt', default=12,
            help='Acquisition interval of time-series in Days.')
    parser.add_argument('-r', '--signal_linear_rate', type=float, dest='linearRate', default=1.0,
            help='Linear rate of the simulated signal over time. [radian/year]')
    parser.add_argument('-s', '--signal_random_sigma', type=float, dest='randomStd', default=1.0,
            help='Sigma od the random phase behaviour over time. [radian]')
    parser.add_argument('-k', '--signal_seasonal_mode', type=float, dest='seasonalMode', default=1,
            help='Seasonal mode of the phase behaviour over time. annual = 1, semiannual=2')
    parser.add_argument('-T', '--decorrelation_tau0', type=float, dest='Tau0', default=72,
            help='Tau0 in the exponential decay decorrelation')
    parser.add_argument('-g', '--decorrelation_gamma0', type=float, dest='gamma0', default=0.95,
            help='gamma0 in the exponential decay decorrelation')
    parser.add_argument('-d', '--decorrelation_model', type=str, dest='decorModel', default="e",
            help='decorrelation model: "e": exponential, "es": exponential + seasonal, "er": exponential + random')
    parser.add_argument('-n', '--number_of_neighbor_samples', type=int, dest='neighborSamples', default=200,
            help='Number of samples in one neighborhood. i.e., Pixels which are statistically homogeneous')
    parser.add_argument('-N', '--number_of_nearestNeighbor_epochs', type=int, dest='nearestNeighborStBAS', default=4,
            help='Number of nearest neighbor acquisitions for simulating small temporal baselines')

    return parser

def cmdLineParse(iargs = None):
    '''
    Command line parser.
    '''

    parser = createParser()
    inps = parser.parse_args(args=iargs)

    return inps



def simulate_coherence_matrix_random(t, gamma0, Tau0, ph):
    length = t.shape[0]
    C = np.ones((length,length), dtype=np.complex64)
    for ii in range(length):
        r = np.random.randn(length)
        r = np.abs(r)/np.max(np.abs(r))
        for jj in range (ii+1,length):
            C[ii,jj] = r[jj]*gamma0*np.exp((t[ii]-t[jj])/Tau0)*np.exp(1j*(ph[ii] - ph[jj]))
            C[jj,ii] = np.conj(C[ii,jj])

    return C

def simulate_coherence_matrix_seasonal(t, gamma0, Tau0, ph):
    length = t.shape[0]
    C = np.ones((length,length), dtype=np.complex64)
    s = np.cos(2.0*np.pi*(t)/120.0)
    S = (s-np.min(s))/(np.max(s)-np.min(s))

    for ii in range(length):
        #s = np.cos(2.0*np.pi*(t-t[ii])/120.0)
        #S = (s-np.min(s))/(np.max(s)-np.min(s))
        S1 = S*S[ii]
        ind = S1<0.1
        S1[ind] = 0.1
        for jj in range (ii+1,length):
            #C[ii,jj] = gamma0*np.exp((t[ii]-t[jj])/Tau0)*np.exp(1j*(ph[ii] - ph[jj]))
            #C[ii,jj] = S[ii]*S[jj]*np.exp(1j*(ph[ii] - ph[jj]))
            C[ii,jj] = S1[jj]*gamma0*np.exp((t[ii]-t[jj])/Tau0)*np.exp(1j*(ph[ii] - ph[jj]))
            C[jj,ii] = np.conj(C[ii,jj])

    return C

def simulate_coherence_matrix_exponential(t, gamma0, Tau0, ph):
    length = t.shape[0]
    C = np.ones((length,length), dtype=np.complex64)
    for ii in range(length):
        for jj in range (ii+1,length):
            C[ii,jj] = gamma0*np.exp((t[ii]-t[jj])/Tau0)*np.exp(1j*(ph[ii] - ph[jj]))
            C[jj,ii] = np.conj(C[ii,jj])

    return C

def simulate_noise(corr_matrix):
    N = corr_matrix.shape[0]

    nsar = corr_matrix.shape[0]
    w, v = np.linalg.eigh(corr_matrix)
    msk = (w < 1e-3)
    w[msk] = 0.
    #corr_matrix =  np.dot(v, np.dot(np.diag(w), np.matrix.getH(v)))

    #C = np.linalg.cholesky(corr_matrix)
    C = np.dot(v, np.dot(np.diag(np.sqrt(w)), np.matrix.getH(v)))
    Z = (np.random.randn(N) +1j*np.random.randn(N)) / np.sqrt(2)
    noise = np.dot(C,Z)

    return noise



def covarinace(C1, C2):

    A1 = np.sum(np.abs(C1)**2)
    A2 = np.sum(np.abs(C2)**2)

    cov = np.sum(C1*np.conjugate(C2))/(np.sqrt(A1)*np.sqrt(A2))
    return cov


def simulate_neighborhood_stack(corr_matrix, neighborSamples=200):
    numberOfSlc = corr_matrix.shape[0]
    # A 2D matrix for a neighborhood over time. Each column is the neighborhood caomplex data for each acquisition date

    neighbor_stack = np.zeros((neighborSamples , numberOfSlc), dtype=np.complex64)

    for ii in range(neighborSamples):
        cpxSLC = simulate_noise(corr_matrix)
        neighbor_stack[ii,:] = cpxSLC

    return neighbor_stack


def compute_covariance_matrix(neighbor_stack):
    numberOfSamples , numberOfSlc = neighbor_stack.shape
    cov_mat = np.zeros((numberOfSlc , numberOfSlc), dtype=np.complex64)
    for ti in range(numberOfSlc):
        for  tj in range(ti+1, numberOfSlc):
            cov = covarinace(neighbor_stack[:,ti], neighbor_stack[:,tj])
            cov_mat[ti,tj] = cov
            cov_mat[tj,ti] = np.conjugate(cov)
        cov_mat[ti,ti] = 1.0

    return cov_mat

def wrap_phase(ph):
    ph = ph - np.round(ph/(2*np.pi)) * 2*np.pi
    return ph



def is_semi_definite(M):
    ###Check if matrix M is positive semi-definite
    w1, v1 = np.linalg.eigh(M)
    msk = (w1 < 0)

    if np.sum(msk)>0:
        print("Matrix is NOT positive semi-definite")
        print("Range of eigen values {0} to {1}".format(np.min(w1),np.max(w1)))
    else:
        print("Matrix is positive semi-definite")
        print("Range of eigen values {0} to {1}".format(np.min(w1),np.max(w1)))
    return None


def plot_result(t, signal_phase, ph, ph_SB, cov_matrix, cov_matrix_SB, corr_mat, Gmat, Gmat_SB):

    fig = plt.figure()

    ax = fig.add_subplot(2,3,1)
    ax.imshow(np.abs(corr_mat), vmin =0, vmax = 1)
    ax.set_title('simulated correlation matrix')

    ax = fig.add_subplot(2,3,2)
    ax.imshow(np.abs(cov_matrix), vmin =0, vmax = 1)
    ax.set_title('Estimated full correlation matrix')

    ax = fig.add_subplot(2,3,3)
    ax.imshow(np.abs(cov_matrix_SB), vmin =0, vmax = 1)
    ax.set_title('Estimated nearest neighbor corr matrix')


    ax = fig.add_subplot(2,3,4)
    ax.imshow(np.angle(corr_mat), vmin=-np.pi, vmax=np.pi)
    ax.set_title('simulated correlation phase')

    ax = fig.add_subplot(2,3,5)
    ax.imshow(np.angle(cov_matrix), vmin=-np.pi, vmax=np.pi)
    ax.set_title('Estimated correlation phase from full cov mat')

    ax = fig.add_subplot(2,3,6)
    ax.imshow(np.angle(cov_matrix_SB), vmin=-np.pi, vmax=np.pi)
    ax.set_title('estimated correlation phase from SB cov matrix')


    #################################

    w, v = np.linalg.eigh(Gmat)
    ph0 = np.angle(v[:,0])
    ph0 = wrap_phase(ph0 - ph0[0])

    w_SB, v_SB = np.linalg.eigh(Gmat_SB)
    ph0_SB = np.angle(v_SB[:,0])
    ph0_SB = wrap_phase(ph0_SB - ph0_SB[0])


    #print('simulated phase:')
    #print(signal_phase)

    #print('estimated phase (MLE) using full cov mat:')
    #print(ph)

    #print('estimated phase (MLE) using SB cov mat:')
    #print(ph_SB)

    #print('first guess phase:')
    #print(ph0)

    #print('phase from single Master Network (PS):')
    #print(np.angle(cov_matrix[:,0]))


    fig = plt.figure()
    ax = fig.add_subplot(1,2,1)
    ax.plot(t, signal_phase)
    ax.plot(t, ph)
    #ax.plot(t, np.angle(cov_matrix[:,0]))
    ax.plot(t, ph0)
    ax.plot(t, ph_SB)
    ax.legend(['simulated', 'estimated MLE', 'Eig', 'SBt'])
    ax.set_xlabel('time [days]')
    ax.set_ylabel('phase [rad]')
    ax.set_ylim([-4,4])
    #plt.plot(t, np.angle(cov_matrix[:,0]))
    #plt.plot(t, np.angle(cov_matrix[:,10]))
    ax = fig.add_subplot(1,2,2)
    ax.plot(t, signal_phase - signal_phase)
    ax.plot(t, ph - signal_phase)
    #ax.plot(t, np.angle(cov_matrix[:,0]) - signal_phase)
    ax.plot(t, ph0 - signal_phase)
    ax.plot(t, ph_SB - signal_phase)
    ax.set_ylim([-6.5,6.5])
    ax.legend(['simulated', 'estimated MLE', 'Eig', 'SBt'])
    #ax.plot(t,np.abs(cov_matrix[:,0]),'o')
    #plt.imshow(np.abs(cov_matrix))
    plt.show()

def simulate_covariance_matrix(t, signal_phase, gamma0, Tau0, method='exponential'):
    #simulating a full covariance matrix
    if method == 'e':  #'exponential':
        corr_mat = simulate_coherence_matrix_exponential(t, gamma0, Tau0, signal_phase)

    elif method == 'es': #'exponential_seasonal':
        corr_mat = simulate_coherence_matrix_seasonal(t, gamma0, Tau0, signal_phase)

    elif method == 'er': #'exponential_random':
        corr_mat = simulate_coherence_matrix_random(t, gamma0, Tau0, signal_phase)

    return corr_mat

def simulate_phase_timeSeries(time_series_length=365, acquisition_interval=12, signal_rate=1, std_random=0, k=1):
    # time_series_length: length of time-series in days
    # acquisition_interval: time-differense between subsequent acquisitions (days)
    # signal_rate: linear rate of the signal (rad/year)
    # k: seasonal parameter, 1 for annaual and  2 for semi-annual
    t = np.arange(0, time_series_length, acquisition_interval)
    signal_phase =  signal_rate*(t-t[0])/365.0
    if k>0:
       seasonal = np.sin(2*np.pi*k*t/365.) + np.cos(2*np.pi*k*t/365.0)
       signal_phase = signal_phase + seasonal

    # adding random temporal signal (which simulates atmosphere + DEM error + ...)
    signal_phase = signal_phase + std_random*np.random.randn(len(t))
    signal_phase = signal_phase - signal_phase[0]
    # wrap the phase to -pi to p
    signal_phase = wrap_phase(signal_phase)

    return signal_phase, t

def make_StBAS_covariance_mat(cov_matrix, numberOfNearestNeighbor=4):
    # Creating a StBAS covariance matrix from full covarinace matrix
    # StBAS covariance matrix: a covarinace matrix that has value only
    #                           around the main diagonal with a buffer of
    #                           nearest neighbor
    numSar = cov_matrix.shape[0]
    cov_matrix_SB = cov_matrix.copy()
    for i in range(numSar):
        for j in range(i + numberOfNearestNeighbor, numSar):
            cov_matrix_SB[i,j] = 0.0+ 1j*0.0
            cov_matrix_SB[j,i] = 0.0+ 1j*0.0

    return cov_matrix_SB


def main(iargs=None):
    inps = cmdLineParse(iargs)
    # setting up the simulation parameters

    signal_phase, t = simulate_phase_timeSeries(time_series_length=inps.timeLength, acquisition_interval=inps.dt, signal_rate=inps.linearRate, std_random=inps.randomStd, k=inps.seasonalMode)

    corr_mat = simulate_covariance_matrix(t, signal_phase, gamma0 = inps.gamma0, Tau0 = inps.Tau0, method=inps.decorModel)

    print("***** Simulating the neighbothood stack ******")
    neighbor_stack = simulate_neighborhood_stack(corr_mat, neighborSamples = inps.neighborSamples)

    ##Estimated cov_matrix
    cov_matrix = compute_covariance_matrix(neighbor_stack)

    ##################################
    # change covariance matrix to Small temporal Baseline network
    # (Keeping only 4 nearest neighbor pairs)
    cov_matrix_SB = make_StBAS_covariance_mat(cov_matrix, numberOfNearestNeighbor=inps.nearestNeighborStBAS)

    ##################################
    #MLE estimation of wrapped phase time-series using Full cov matrix
    print("******************************************")
    print("MLE estimation using full cov matrix")
    ph, Gmat = mle4(cov_matrix)

    print("******************************************")
    print("MLE estimation using SB cov matrix")
    ph_SB, Gmat_SB = mle4(cov_matrix_SB)

    plot_result(t, signal_phase, ph, ph_SB, cov_matrix, cov_matrix_SB, corr_mat, Gmat, Gmat_SB)

    print("******************************************")

    #Change  the covarainace matrix to have the known coherence
    cov_matrix = np.abs(corr_mat)*np.exp(1.0j*np.angle(cov_matrix))

    cov_matrix_SB = make_StBAS_covariance_mat(cov_matrix, numberOfNearestNeighbor=4)

    print("MLE estimation using full cov matrix (simulated known coherence is used)")
    ph, Gmat = mle4(cov_matrix)

    plot_result(t, signal_phase, ph, ph_SB, cov_matrix, cov_matrix_SB, corr_mat, Gmat, Gmat_SB)



if __name__ == '__main__':

    main()
