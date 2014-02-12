

import sys


""" This script will take a chi2.txt file as an input, which is an output
    of the HRMC simulation. However, first run double_to_single_lines.py
    on th chi2.txt file. This script will then parse that file for the
    step, energy, and chi2 columns.
    It will then output the normalized versions of each of these. """

def main():
    chi2file = sys.argv[1]
    
    print("Reading in content...")
    with open(chi2file,'r') as f:
        content = f.readlines()

    print("Reducing content...")
    content = [line.split()[0:3] for line in content]
    
    print("Parsing content...")
    print("Parsing step...")
    header = content.pop(0)
    step = [int(line[0]) for line in content]
    print("Parsing chi2...")
    chi2 = [float(line[1]) for line in content]
    print("Parsing energy...")
    energy = [float(line[2]) for line in content]

    minE = min(energy)
    energy = [x - minE for x in energy]

    print("Calculating scales...")
    minchi2 = min(chi2)
    maxchi2 = max(chi2)
    minE = min(energy)
    maxE = max(energy)

    print("Rescaling...")
    chi2 = [x/maxchi2 for x in chi2]
    energy = [x/maxE for x in energy]

    #for i in range(0,len(chi2)):
    #    print('{0}\t{1}'.format(chi2[i],energy[i]))
    #    print('{0}'.format(chi2[i]/energy[i]))

    c_over_e  = [chi2[i]/energy[i] for i in range(1,len(chi2))]
    c_minus_e = [chi2[i]-energy[i] for i in range(1,len(chi2))]
    e_minus_c = [energy[i]-chi2[i] for i in range(1,len(chi2))]

    print(' '.join(header)+', chi2/energy, chi2-energy, energy-chi2')
    for i,s in enumerate(step):
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(s,chi2[i],energy[i],c_over_e[i],c_minus_e[i],e_minus_c[i]))



if __name__ == "__main__":
    main()
