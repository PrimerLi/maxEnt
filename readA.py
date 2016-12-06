import chi
import entropy
import numpy as np

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


def readA(fileName):
    omega = []
    A = []
    ifile = open(fileName, "r")
    for i, string in enumerate(ifile):
        a = string.split()
        omega.append(float(a[0]))
        A.append(float(a[1]))
    ifile.close()
    return omega, A

def main():
    import os
    import sys
    import numpy.linalg
    
    Greal = "G_cc_real.txt"
    Gimag = "G_cc_imag.txt"
    omega_n, G = readFiles(Greal, Gimag)
    Niom = len(omega_n)
    
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
        C_real[i, i] = 0.004**2
        C_imag[i, i] = 0.004**2
    C_real_inv = numpy.linalg.inv(C_real)
    C_imag_inv = numpy.linalg.inv(C_imag)


    alpha = []
    ifile = open("alpha.txt", "r")
    for i, string in enumerate(ifile):
        alpha.append(float(string))
    ifile.close()
    
    spectrals = []
    probability = []
    chi_values = []
    entropy_values = []
    for i in range(len(alpha)):
        print alpha[i]
        fileName = "A_updated_alpha_" + str(alpha[i]) + ".txt"
        omega, A = readA(fileName)
        spectrals.append(np.asarray(A))
        chi_values.append(chi.chi(G, A, omega_n, omega, C_real_inv, C_imag_inv))
        entropy_values.append(entropy.entropy(omega, A))
        probability.append(np.exp(alpha[i]*entropy_values[i] - 0.5*chi_values[i]))

    spectral_mean = np.zeros(len(spectrals[0]))
    for i in range(1, len(spectrals)):
        spectral_mean[:] = spectral_mean[:] + ((-alpha[i] + alpha[i-1])*probability[i]*spectrals[i])[:]
    s = 0.0
    for i in range(1, len(alpha)):
        s = s + (-alpha[i] + alpha[i-1])*probability[i]
    print "s = ", s
    ofile = open("Chi_log_alpha.txt", "w")
    for i in range(len(alpha)):
        ofile.write(str(np.log(alpha[i])) + "   " + str(chi_values[i]) + "\n")
    ofile.close()
    ofile = open("entropy_log_alpha.txt", "w")
    for i in range(len(alpha)):
        ofile.write(str(np.log(alpha[i])) + "    " + str(entropy_values[i]) + "\n")
    ofile.close()
    ofile = open("P_log_alpha.txt", "w")
    for i in range(len(alpha)):
        ofile.write(str(np.log(alpha[i])) + "   " + str(probability[i]) + "\n")
    ofile.close()
    ofile = open("bryan.txt", "w")
    for i in range(len(omega)):
        ofile.write(str(omega[i]) + "    " + str(spectral_mean[i]/s) + "\n")
    ofile.close()
    return 0

main()
