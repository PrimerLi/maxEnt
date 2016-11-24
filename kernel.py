import numpy as np

def K_real(omega_n, omega):
    return -omega/(omega_n**2 + omega**2)*1/(2*np.pi)

def K_imag(omega_n, omega):
    return -omega_n/(omega_n**2 + omega**2)*1.0/(2*np.pi)
