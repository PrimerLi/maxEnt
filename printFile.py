def printFile(x, y, fileName):
    ofile = open(fileName, "w")
    for i in range(len(x)):
        ofile.write(str(x[i]) + "    " + str(y[i]) + "\n")
    ofile.close()
