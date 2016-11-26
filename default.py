import printFile

def gauss(x, mu, sigma):
    import numpy as np
    return 1.0/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))

def trapezoid(f, lower, upper):
    s = 0.0
    N = 1000
    dx = (upper - lower)/float(N)
    for i in range(N):
        xa = lower + i*dx
        xb = lower + (i+1)*dx
        s = s + 0.5*(f(xa) + f(xb))*dx
    return s

def D(omega):
    return gauss(omega, 0, 1)

def A(omega): 
    return 0.5*gauss(omega, -1.2, 0.5) + 0.5*gauss(omega, 1.2, 0.4)

def defaultG(omega_n):
    import numpy as np
    
    lower = -20
    upper = -lower
    N = 600
    domega = (upper - lower)/float(N)

    def integrand(omega_n, omega):
        return 0.5/np.pi*(A(omega))/(-omega + omega_n*1j)

    s = 0.0 + 0.0j
    for i in range(N):
        omega_a = lower + i*domega
        omega_b = lower + (i+1)*domega
        s = s + 0.5*domega*(integrand(omega_n, omega_a) + integrand(omega_n, omega_b))
    return s

def main():
    import os
    import sys
    import numpy as np

    beta = 10.0
    Niom = 100
    omega_n = []
    for nw in range(Niom):
        omega_n.append((2*nw+1)*np.pi/beta)

    G_real = []
    G_imag = []
    for nw in range(Niom):
        temp = defaultG(omega_n[nw])
        G_real.append(temp.real)
        G_imag.append(temp.imag)

    randomNumbers = np.random.normal(0, 0.0001, Niom)
    for nw in range(len(G_real)):
        G_real[nw] = G_real[nw] + randomNumbers[nw]
    randomNumbers = np.random.normal(0, 0.00015, Niom)
    for nw in range(len(G_imag)):
        G_imag[nw] = G_imag[nw] + randomNumbers[nw]
    printFile.printFile(omega_n, G_real, "G_real.txt")
    printFile.printFile(omega_n, G_imag, "G_imag.txt")
    ofile = open("G.txt", "w")
    for i in range(len(omega_n)):
        ofile.write(str(omega_n[i]) + "    " + str(G_real[i]) + "    " + str(G_imag[i]) + "\n")
    ofile.close()
    ofile = open("G_error.txt", "w")
    for i in range(len(omega_n)):
        ofile.write(str(omega_n[i]) + "    " + str(0.0001) + "   " + str(0.00015) + "\n")
    ofile.close()

    omega = []
    spectral = []
    omega_lower = -5
    omega_upper = 5
    Nomega = 60
    delta = (omega_upper - omega_lower)/float(Nomega)
    for i in range(Nomega+1):
        omega.append(omega_lower + i*delta)
    for i in range(len(omega)):
        spectral.append(A(omega[i]))
    printFile.printFile(omega, spectral, "spectral.txt")
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
