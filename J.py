import default
import kernel

def J(alpha, A, omega_index, omega_index_prime, omega_n, omega, C_real, C_imag):
    import numpy as np
    import numpy.linalg

    C_real_inv = numpy.linalg.inv(C_real)
    C_imag_inv = numpy.linalg.inv(C_imag)

    result = 0.0
    if (omega_index == omega_index_prime):
        result = -alpha/A[omega_index]

    Niom = len(omega_n)
    vector_right = np.zeros(Niom)
    vector_left = np.zeros(Niom)
    for nw in range(Niom):
        vector_right[nw] = -kernel.K_real(omega_n[nw], omega[omega_index_prime])
        vector_left[nw] = kernel.K_real(omega_n[nw], omega[omega_index])
    result = result + vector_left.dot(C_real_inv).dot(vector_right)

    vector_right = np.zeros(Niom)
    vector_left = np.zeros(Niom)
    for nw in range(Niom):
        vector_right[nw] = -kernel.K_imag(omega_n[nw], omega[omega_index_prime])
        vector_left[nw] = kernel.K_imag(omega_n[nw], omega[omega_index])
    result = result + vector_left.dot(C_imag_inv).dot(vector_right)
    return result
