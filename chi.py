import kernel

def chi(G, A, omega_n, omega, C_real_inv, C_imag_inv):
    import numpy as np
    import numpy.linalg

    Niom = len(G)
    Nomega = len(omega)
    G_real = np.zeros(Niom)
    G_imag = np.zeros(Niom)
    for i in range(Niom):
        G_real[i] = G[i].real
        G_imag[i] = G[i].imag

    result = 0.0
    vector_left = np.zeros(Niom)
    vector_right = np.zeros(Niom)
    for nw in range(Niom):
        KA = 0.0
        domega = omega[1] - omega[0]
        for i in range(Nomega):
            KA = KA + kernel.K_real(omega_n[nw], omega[i])*A[i]*domega
        vector_left[nw] = G_real[nw] - KA
    vector_right[:] = vector_left[:]
    result = result + vector_left.dot(C_real_inv).dot(vector_right)

    vector_left = np.zeros(Niom)
    for nw in range(Niom):
        KA = 0.0
        domega = omega[1] - omega[0]
        for i in range(Nomega):
            KA = KA + kernel.K_imag(omega_n[nw], omega[i])*A[i]*domega
        vector_left[nw] = G_imag[nw] - KA
    result = result + vector_left.dot(C_imag_inv).dot(vector_left)
    return result
