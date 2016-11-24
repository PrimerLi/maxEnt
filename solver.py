import default
import kernel
import f
import J
import chi

def readFiles():
    import os
    import sys
    import numpy as np
    import sys

    omega_n = []
    G = []
    G_real = []
    G_imag = []

    try:
        ifile = open("G_real.txt", "r")
    except:
        sys.exit("G_real.txt does not exist. ")
    for index, string in enumerate(ifile):
        a = string.split()
        omega_n.append(float(a[0]))
        G_real.append(float(a[1]))
    ifile.close()

    try:
        ifile = open("G_imag.txt", "r")
    except:
        sys.exit("G_imag.txt does not exist. ")
    for index, string in enumerate(ifile):
        a = string.split()
        G_imag.append(float(a[1]))
    ifile.close()

    for i in range(len(G_real)):
        G.append(G_real[i] + G_imag[i]*1j)

    return omega_n, G

def diff(a, b):
    import numpy as np
    s = 0.0
    for i in range(len(a)):
        s = s + (a[i] - b[i])**2
    return np.sqrt(s)

def newton(alpha, G, omega_n, omega, A_initial, C_real, C_imag):
    import numpy as np
    import numpy.linalg

    Nomega = len(omega)
    Niom = len(omega_n)
    Jacobian = np.zeros((Nomega, Nomega))
    function = np.zeros(Nomega)
    A_updated = np.zeros(Nomega)
    iterationMax = 30
    counter = 0

    eps = 0.001
    while(True):
        counter = counter + 1 
        if (counter > iterationMax):
            break
        for i in range(Nomega):
            for j in range(Nomega):
                Jacobian[i, j] = J.J(alpha, A_initial, i, j, omega_n, omega, C_real, C_imag)
        for i in range(Nomega):
            function[i] = f.f(alpha, G, A_initial, i, omega_n, omega, C_real, C_imag)
        A_updated[:] = A_initial[:] - numpy.linalg.inv(Jacobian).dot(function)[:]
        error = diff(A_initial, A_updated)
        if (error < eps):
            break
        print "counter = ", counter, ", diff = ", error
        A_initial[:] = A_updated[:]

    return A_initial

def printFile(x, y, fileName):
    ofile = open(fileName, "w")
    for i in range(len(x)):
        ofile.write(str(x[i]) + "    " + str(y[i]) + "\n")
    ofile.close()

def main():
    import os
    import sys
    import numpy as np
    import numpy.linalg

    omega_n, G = readFiles()
    Niom = len(omega_n)

    N = 40
    omega_lower = -5
    omega_upper = -omega_lower
    domega = (omega_upper - omega_lower)/float(N)
    omega = np.zeros(N+1)
    for i in range(N+1):
        omega[i] = omega_lower + i*domega
    Nomega = len(omega)
    A_initial = np.zeros(Nomega)
    for i in range(len(A_initial)):
        A_initial[i] = default.D(omega[i])
    
    printFile(omega, A_initial, "A_initial.txt")

    C_real = np.zeros((Niom, Niom))
    C_imag = np.zeros((Niom, Niom))
    for i in range(Niom):
        C_real[i, i] = 0.001**2
        C_imag[i, i] = 0.0015**2

    if (True):
        alpha = []
        chi_values = []
        for i in range(8):
            alpha.append(1.0e1*(i+2))
        for i in range(len(alpha)):
            print "alpha = ", alpha[i]
            A_updated = newton(alpha[i], G, omega_n, omega, A_initial, C_real, C_imag)
            chi_values.append(chi.chi(G, A_updated, omega_n, omega, C_real, C_imag))
            printFile(omega, A_updated, "A_updated_alpha_" + str(alpha[i]) + ".txt")
        printFile(alpha, chi_values, "chi_alpha.txt")
    if (False):
        alpha = 1.5e4
        A_updated = newton(alpha, G, omega_n, omega, A_initial, C_real, C_imag)
        printFile(omega, A_updated, "A_updated_alpha_" + str(alpha) + ".txt")
        print chi.chi(G, A_updated, omega_n, omega, C_real, C_imag)
    return 0

main()
