import os

def main():
    secperstep = []
    cores = []
    for name in os.listdir(os.getcwd()):
        if name[:4] == 'time' and name[len(name)-4:] == '.txt':
            print(name)
            f = open( name, 'r' )
            line = readl(f)
            #
            # Change the below core = line!!! if using mpi+omp to multiply by the number of cores per node!
            #
            core = int(line[0])
            iarry = []
            timearry = []
            line = readl(f) # Comment line
            line = readl(f) # Comment line wrap-around
            line = readl(f)
            while(line != []):
                iarry.append(int(line[0]))
                timearry.append(float(line[1]))
                line = readl(f)
                line = readl(f)
            f.close()
            istart = iarry.pop(0)
            timearry.pop(0)
            iarry = [ x - istart for x in iarry ]
            for i in range(0,len(timearry)):
                timearry[i] = timearry[i]/iarry[i]
            avg = sum(timearry)/len(timearry)
            #print(core, avg)
            cores.append(core)
            secperstep.append(avg)

    cores.sort()
    secperstep.sort()
    secperstep.reverse()
    steppersec = [ 1/x for x in secperstep ]
    steppersec.sort()
    conversion = cores[0]/steppersec[0]
    #print conversion
    #print cores
    #print secperstep
    #print steppersec
    speedup = [ conversion*x for x in steppersec ]
    f = open("output_speedup.txt", 'w')
    f.write("cores sec/step step/sec speedup\n")
    for i in range(0,len(cores)):
        f.write(str(cores[i]) +" "+ str(secperstep[i]) +" "+ str(steppersec[i]) +" "+ str(speedup[i]) + '\n')


def readl(f):
    line = f.readline()
    line = line.split(" ")
    while( '' in line ):
        line.remove('')
    return line

if __name__ == "__main__":
    main();
