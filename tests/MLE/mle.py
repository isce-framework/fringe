import numpy as np
import matplotlib.pyplot as plt

def wrap_phase(ph):
    ph = ph - np.round(ph/(2*np.pi)) * 2*np.pi
    return ph




def estimate_MLE_5(cor_matrix):

    ###Correlation matrix
    ###This matrix needs to be positive semi-definite
    nsar = cor_matrix.shape[0]
    w, v = np.linalg.eigh(cor_matrix)
    msk = (w < 1e-6)
    w[msk] = 0.
    cor_matrix =  np.dot(v, np.dot(np.diag(w), np.matrix.getH(v)))

    print('Eigen value properties of COV')
    print('Range: ', np.min(w), np.max(w))
    print('Number used: {0} out of {1}'.format(nsar - np.sum(msk), nsar))


    ###Absolute value of the correlation matrix
    ###This matrix needs to be positive semi-definite
    coh = np.abs(cor_matrix)
    coh = 0.5 * (coh + coh.T)
    w1, v1 = np.linalg.eigh(coh)

    msk = (w1 > 1e-6)
    wdiag1 = np.zeros(w1.size)
    wdiag1[msk] = 1 / w1[msk]
    coh_inv = np.dot(v1, np.dot(np.diag(wdiag1), np.matrix.getH(v1)))

    print('Eigen value properties of COH')
    print('Range: ', np.min(w1), np.max(w1))
    print('Range inv: ', np.min(wdiag1), np.max(wdiag1))
    print('Number used: {0} out of {1}'.format(np.sum(msk), nsar))



    ### Generating the Hadamard product
    ### This is the matrix used for the MLE iterations
    Gmat = coh_inv * cor_matrix
    Gmat = 0.5 * (Gmat + np.matrix.getH(Gmat))
    w, v = np.linalg.eigh(Gmat)
    #w, v = np.linalg.eigh(cor_matrix)

    ###Init with zeros
    #ph0 = np.zeros(nsar)

    ###Init with random
    #ph0 = np.random.randn(nsar)

    ###Init with cor_matrix
    #ph0 = -np.angle(cor_matrix[:,0])

    ###Init with eigen

    ph0 = np.angle(v[:,0])

    print('Eigen value properties of Hadamard')
    print('Hadamard: ' ,np.min(w), np.max(w))


    N = coh.shape[0]
    ind = range(N)
    ph = np.zeros(ph0.shape)
    ph0 = wrap_phase(ph0 - ph0[0])
    iteration = 0
    plt.figure()
    niter = 25
    while iteration < niter:
        print('iteration: ', iteration)
        for p in range(N):
            psum = 0.0+1.j*0.0

            for n in range(N):
                if (n != p):
                    psum += (coh_inv[p,n]*cor_matrix[p,n]*np.exp(1j*ph0[n]))

            ph[p] = np.angle(psum)*(np.abs(psum)>1e-5)

        ph = np.angle(np.exp(1j*ph)*np.exp(-1j*ph[0]))
        #ph = ph - ph[0]
        residual = (np.angle(np.exp(1j*ph)*np.exp(-1j*ph0)))
        #print(residual)

        if (iteration==0):
            plt.subplot(niter+1,1,1)
            plt.plot(ph0)

        ph0 = ph.copy()
        iteration+=1

        plt.subplot(niter+1, 1, iteration+1)
        plt.plot(ph)
        if np.max(np.abs(residual)) < 1.0*np.pi/180. :
            iteration = niter

    return ph, Gmat



def estimate_MLE_4(cor_matrix):

    ###Correlation matrix
    ###This matrix needs to be positive semi-definite
    nsar = cor_matrix.shape[0]
    w, v = np.linalg.eigh(cor_matrix)
    msk = (w < 1e-3)
    w[msk] = 0.
    cor_matrix =  np.dot(v, np.dot(np.diag(w), np.matrix.getH(v)))

    print('Eigen value properties of COV')
    print('Range: ', np.min(w), np.max(w))
    print('Number used: {0} out of {1}'.format(nsar - np.sum(msk), nsar))


    ###Absolute value of the correlation matrix
    ###This matrix needs to be positive semi-definite
    coh = np.abs(cor_matrix)
    coh = 0.5 * (coh + coh.T)
    w1, v1 = np.linalg.eigh(coh)

    msk = (w1 > 1e-3)
    wdiag1 = np.zeros(w1.size)
    wdiag1[msk] = 1 / w1[msk]
    coh_inv = np.dot(v1, np.dot(np.diag(wdiag1), np.matrix.getH(v1)))

    print('Eigen value properties of COH')
    print('Range: ', np.min(w1), np.max(w1))
    print('Range inv: ', np.min(wdiag1), np.max(wdiag1))
    print('Number used: {0} out of {1}'.format(np.sum(msk), nsar))



    ### Generating the Hadamard product
    ### This is the matrix used for the MLE iterations
    Gmat = coh_inv * cor_matrix
    Gmat = 0.5 * (Gmat + np.matrix.getH(Gmat))
    w, v = np.linalg.eigh(Gmat)


    ###Init with zeros
    #ph0 = np.zeros(nsar)

    ###Init with random
    #ph0 = np.random.randn(nsar)

    ###Init with cor_matrix
    #ph0 = -np.angle(cor_matrix[:,0])

    ###Init with eigen
    ph0 = np.angle(v[:,0])


    print('Eigen value properties of Hadamard')
    print('Hadamard: ' ,np.min(w), np.max(w))


    N = coh.shape[0]
    ind = range(N)
    ph = np.zeros(ph0.shape)
    ph0 = wrap_phase(ph0 - ph0[0])
    iteration = 0
    plt.figure()
    niter = 25
    while iteration < niter:
        print('iteration: ', iteration)
        for p in range(N):
            psum = 0.0+1.j*0.0

            for n in range(N):
                if (n != p):
                    psum += (coh_inv[p,n]*cor_matrix[p,n]*np.exp(1j*ph0[n]))

            ph[p] = np.angle(psum)*(np.abs(psum)>1e-5)

        ph = np.angle(np.exp(1j*ph)*np.exp(-1j*ph[0]))
        #ph = ph - ph[0]
        residual = (np.angle(np.exp(1j*ph)*np.exp(-1j*ph0)))
        #print(residual)

        if (iteration==0):
            plt.subplot(niter+1,1,1)
            plt.plot(ph0)

        ph0 = ph.copy()
        iteration+=1

        plt.subplot(niter+1, 1, iteration+1)
        plt.plot(ph)
        if np.max(np.abs(residual)) < 10.0*np.pi/180. :
            iteration = niter

    return ph, Gmat


def estimate_MLE_3(cor_matrix):

    ###Correlation matrix
    ###This matrix needs to be positive semi-definite
    nsar = cor_matrix.shape[0]
    w, v = np.linalg.eigh(cor_matrix)
    msk = (w < 1e-3)
    w[msk] = 0.
    cor_matrix =  np.dot(v, np.dot(np.diag(w), np.matrix.getH(v)))

    print('Eigen value properties of COV')
    print('Range: ', np.min(w), np.max(w))
    print('Number used: {0} out of {1}'.format(nsar - np.sum(msk), nsar))


    ###Absolute value of the correlation matrix
    ###This matrix needs to be positive semi-definite
    coh = np.abs(cor_matrix)
    coh = 0.5 * (coh + coh.T)
    w1, v1 = np.linalg.eigh(coh)

    msk = (w1 > 1e-3)
    wdiag1 = np.zeros(w1.size)
    wdiag1[msk] = 1 / w1[msk]
    coh_inv = np.dot(v1, np.dot(np.diag(wdiag1), np.matrix.getH(v1)))

    print('Eigen value properties of COH')
    print('Range: ', np.min(w1), np.max(w1))
    print('Range inv: ', np.min(wdiag1), np.max(wdiag1))
    print('Number used: {0} out of {1}'.format(np.sum(msk), nsar))



    ### Generating the Hadamard product
    ### This is the matrix used for the MLE iterations
    Gmat = coh_inv * cor_matrix
    Gmat = 0.5 * (Gmat + np.matrix.getH(Gmat))
    w, v = np.linalg.eigh(Gmat)


    ###Init with zeros
    #ph0 = np.zeros(nsar)

    ###Init with random
    #ph0 = np.random.randn(nsar)

    ###Init with cor_matrix
    #ph0 = -np.angle(cor_matrix[:,0])

    ###Init with eigen
    ph0 = np.angle(v[:,0])


    print('Eigen value properties of Hadamard')
    print('Hadamard: ' ,np.min(w), np.max(w))


    N = coh.shape[0]
    ind = range(N)
    ph = np.zeros(ph0.shape)
    ph0 = wrap_phase(ph0 - ph0[0])
    iteration = 0
    plt.figure()
    niter = 25
    while iteration < niter:
        print('iteration: ', iteration)
        for p in range(N):
            psum = 0.0+1.j*0.0

            for n in range(N):
                if (n != p):
                    psum += (coh_inv[p,n]*cor_matrix[p,n]*np.exp(1j*ph0[n]))

            ph[p] = np.angle(psum)

        ph = np.angle(np.exp(1j*ph)*np.exp(-1j*ph[0]))
        #ph = ph - ph[0]
        residual = (np.angle(np.exp(1j*ph)*np.exp(-1j*ph0)))
        #print(residual)


        if (iteration==0):
            plt.subplot(niter+1,1,1)
            plt.plot(ph0)

        ph0 = ph.copy()
        iteration+=1

        plt.subplot(niter+1, 1, iteration+1)
        plt.plot(ph)

    return ph, Gmat







