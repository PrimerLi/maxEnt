import default
import kernel

def f(alpha, G, A, omega_index, omega_n, omega, C_real, C_imag):
    import numpy as np
    import numpy.linalg 

    result = -alpha*(1 + np.log(A[omega_index]/default.D(omega[omega_index])))

    C_real_inv = numpy.linalg.inv(C_real)
    C_imag_inv = numpy.linalg.inv(C_imag)
    
    Niom = len(omega_n)

    vector_right = np.zeros(Niom)
    for nw in range(Niom):
        KA = 0.0
        domega = omega[1] - omega[0]
        for i in range(len(omega)):
            KA = KA + kernel.K_real(omega_n[nw], omega[i])*A[i]*domega
        vector_right[nw] = G[nw].real - KA
    vector_left = np.zeros(Niom)
    for nw in range(Niom):
        vector_left[nw] = kernel.K_real(omega_n[nw], omega[omega_index])
    result = result + vector_left.dot(C_real_inv).dot(vector_right)

    vector_right = np.zeros(Niom)
    for nw in range(Niom):
        KA = 0.0
        domega = omega[1] - omega[0]
        for i in range(len(omega)):
            KA = KA + kernel.K_imag(omega_n[nw], omega[i])*A[i]*domega
        vector_right[nw] = G[nw].imag - KA
    vector_left = np.zeros(Niom)
    for nw in range(Niom):
        vector_left[nw] = kernel.K_imag(omega_n[nw], omega[omega_index])
    result = result + vector_left.dot(C_imag_inv).dot(vector_right)
    return result
