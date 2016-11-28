import default
import kernel
import f
import J
import chi

def readFiles(Greal, Gimag):
    import os
    import sys
    import numpy as np
    import sys

    omega_n = []
    G = []
    G_real = []
    G_imag = []

    try:
        ifile = open(Greal, "r")
    except:
        sys.exit(Greal + " does not exist. ")
    for index, string in enumerate(ifile):
        a = string.split()
        omega_n.append(float(a[0]))
        G_real.append(float(a[1]))
    ifile.close()

    try:
        ifile = open(Gimag, "r")
    except:
        sys.exit(Gimag + " does not exist. ")
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

def newton(alpha, G, omega_n, omega, A_initial, C_real_inv, C_imag_inv):
    import numpy as np
    import numpy.linalg

    Nomega = len(omega)
    Niom = len(omega_n)
    Jacobian = np.zeros((Nomega, Nomega))
    function = np.zeros(Nomega)
    A_updated = np.zeros(Nomega)
    iterationMax = 100
    counter = 0

    eps = 0.0001
    while(True):
        counter = counter + 1 
        if (counter > iterationMax):
            break
        for i in range(Nomega):
            for j in range(Nomega):
                Jacobian[i, j] = J.J(alpha, A_initial, i, j, omega_n, omega, C_real_inv, C_imag_inv)
        for i in range(Nomega):
            function[i] = f.f(alpha, G, A_initial, i, omega_n, omega, C_real_inv, C_imag_inv)
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

    Greal = "G_cc_real.txt"
    Gimag = "G_cc_imag.txt"
    omega_n, G = readFiles(Greal, Gimag)
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
    if (not os.path.exists("A_initial.txt")):
        for i in range(len(A_initial)):
            A_initial[i] = default.D(omega[i])
        printFile(omega, A_initial, "A_initial.txt")
    else:
        omega = []
        A_initial = []
        ifile = open("A_initial.txt", "r")
        for index, string in enumerate(ifile):
            a = string.split()
            omega.append(float(a[0]))
            A_initial.append(float(a[1]))
        ifile.close()
        Nomega = len(omega)

    C_real = np.zeros((Niom, Niom))
    C_imag = np.zeros((Niom, Niom))

    ifile = open("CM_cc_real.txt", "r")
    for (index, string) in enumerate(ifile):
        a = string.split()
        rowIndex = int(a[0])-1
        colIndex = int(a[1])-1
        if (True):
            C_real[rowIndex, colIndex] = float(a[2])
    ifile.close()
    ifile = open("CM_cc_imag.txt", "r")
    for (index, string) in enumerate(ifile):
        a = string.split()
        rowIndex = int(a[0])-1
        colIndex = int(a[1])-1
        if (True):
            C_imag[rowIndex, colIndex] = float(a[2])
    ifile.close()
    for i in range(Niom):
        C_real[i, i] = 0.001**2
        C_imag[i, i] = 0.001**2
    C_real_inv = numpy.linalg.inv(C_real)
    C_imag_inv = numpy.linalg.inv(C_imag)

    if (True):
        alpha = []
        chi_values = []
        for i in range(10):
            alpha.append(100*np.exp(-i*0.5))
        for i in range(len(alpha)):
            print "alpha = ", alpha[i]
            A_updated = newton(alpha[i], G, omega_n, omega, A_initial, C_real_inv, C_imag_inv)
            chi_values.append(chi.chi(G, A_updated, omega_n, omega, C_real_inv, C_imag_inv))
            printFile(omega, A_updated, "A_updated_alpha_" + str(alpha[i]) + ".txt")
            os.system("cp A_updated_alpha_" + str(alpha[i]) + ".txt A_initial.txt")
        printFile(alpha, chi_values, "chi_alpha.txt")
    else:
        alpha = 0.01
        print "alpha = ", alpha
        A_updated = newton(alpha, G, omega_n, omega, A_initial, C_real_inv, C_imag_inv)
        printFile(omega, A_updated, "A_updated_alpha_" + str(alpha) + ".txt")
        os.system("cp A_updated_alpha_" + str(alpha) + ".txt" + "A_initial.txt")
        print chi.chi(G, A_updated, omega_n, omega, C_real_inv, C_imag_inv)
    return 0

main()
