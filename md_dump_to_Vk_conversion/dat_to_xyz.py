import os
import time
import sys
from decimal import *

def main():
    try:
        xyz_to_dat(sys.argv[1])
    except:
        print("\nMini help printout:\nThe first (and only) argument must be the .xyz file you want to convert to .dat.")
        print(".dat to .xyz is currently not supported due to the different syntax in .dat files.")
        print("Note also, that no blank lines in the .xyz file are allowed. The format MUST be:\n\nA comment line\nxsize ysize zsize\nList of atoms on following lines\n-1\n")

def dat_to_xyz(name):
    df = open( name, 'r' )
    xyz = open(name[:-3] + 'xyz', 'w')
    line = df.readline() #Comment line
    xyz.write(line)
    line = df.readline() #Blank line
    if(line.strip() != ''): sys.exit("Expected a blank line")
    natoms = int(df.readline().strip().split(" ")[0]) #natoms
    line = df.readline() #3 atom types
    xsize = sum( [ abs(float(x)) for x in df.readline().split(" ")[0:2] ] )
    ysize = sum( [ abs(float(x)) for x in df.readline().split(" ")[0:2] ] )
    zsize = sum( [ abs(float(x)) for x in df.readline().split(" ")[0:2] ] )
    print xsize, ysize, zsize
    xyz.write(str(xsize) + ' ' + str(ysize) + ' ' + str(zsize) + '\n')
    line = df.readline() #Blank line
    if(line.strip() != ''): sys.exit("Expected a blank line")
    line = df.readline() #Blank line
    if(line.strip() != ''): sys.exit("Expected a blank line")
    line = df.readline() #Atoms
    line = df.readline() #Blank line
    if(line.strip() != ''): sys.exit("Expected a blank line")
    line = df.readline().strip() #First atom line

    while(len(line)!=0):
        # The types given in the dump file need to be converted to the actual elements' atomic numbers. In my case, they are all the same, so I will not code the universal way to convert them. You can find the masses for the types in the .dat files and then look up which mass is which element. For my case, 1 = 91.22 = Zr = 40, 2 = 63.55 = Cu = 29, and 3 = 26.98 = Al = 13.
        split_line = line.split(" ")[1:]
        if(split_line[0] == "1"):
            split_line[0] = "40"
        elif(split_line[0] == "2"):
            split_line[0] = "29"
        elif(split_line[0] == "3"):
            split_line[0] = "13"
        xyz.write( " ".join(split_line) + '\n' )
        line = df.readline().strip()
    xyz.write('-1')

def xyz_to_dat(name):
    xyzf = open( name, 'r' )
    line = xyzf.readline() #Comment line.
    line = xyzf.readline() #Model size line.
    line = xyzf.readline() #First atom line.
    xyzf = open( name, 'r' )
    dat = open(name[:-3] + 'dat', 'w')
    line = xyzf.readline() #Get comment line
    dat.write(line) #Write comment line
    dat.write('\n') #Write blank line
    natoms = 0
    line = xyzf.readline() #Get model size
    sline = line.split(' ')
    while( '' in sline): sline.pop('')
    xsize = float(sline[0])
    ysize = float(sline[1])
    zsize = float(sline[2])
    print xsize, ysize, zsize

    # Atoms start now.
    atoms = []
    elements = []
    line = xyzf.readline()
    while(line!= '-1\n'):
        natoms += 1
        sline = line.split(" ")
        while( '' in sline): sline.pop('')
        for i in range(0,len(sline)):
            #print len(sline[i].strip()) - sline[i].strip().find('.')
            getcontext().prec = 6
            sline[i] = str(Decimal(sline[i].strip()))
        if sline[0] not in elements:
            elements.append(sline[0])
        sline[0] = str(elements.index(sline[0])+1)
        sline = [str(natoms)] + sline
        #print sline
        atoms.append(" ".join(sline) + '\n')
        line = xyzf.readline() #Get another atom.

    dat.write('    '+str(natoms)+' atoms\n')
    dat.write('    '+str(len(elements))+' atom types\n')
    dat.write(str(-xsize/2)+' '+str(xsize/2)+' xlo xhi\n')
    dat.write(str(-ysize/2)+' '+str(ysize/2)+' ylo yhi\n')
    dat.write(str(-zsize/2)+' '+str(zsize/2)+' zlo zhi\n')
    dat.write('\n')
    dat.write('Masses\n')
    for elem in elements:
        dat.write(str(elements.index(elem)+1) + ' ' + get_mass(int(elem)) + '\n')
    dat.write('\n')
    dat.write('Atoms\n')
    for aline in atoms:
        dat.write(aline)
    dat.write('\n')


def get_mass(x):
        if(x == 1): return str(1.01)
        if(x == 2): return str(4.00)
        if(x == 3): return str(6.94)
        if(x == 4): return str(9.01)
        if(x == 5): return str(10.81)
        if(x == 6): return str(12.01)
        if(x == 7): return str(14.01)
        if(x == 8): return str(16.00)
        if(x == 9): return str(19.00)
        if(x == 10): return str(20.18)
        if(x == 11): return str(22.99)
        if(x == 12): return str(24.31)
        if(x == 13): return str(26.98)
        if(x == 14): return str(28.09)
        if(x == 15): return str(30.97)
        if(x == 16): return str(32.07)
        if(x == 17): return str(35.45)
        if(x == 19): return str(39.10)
        if(x == 18): return str(39.95)
        if(x == 20): return str(40.08)
        if(x == 21): return str(44.96)
        if(x == 22): return str(47.87)
        if(x == 23): return str(50.94)
        if(x == 24): return str(52.00)
        if(x == 25): return str(54.94)
        if(x == 26): return str(55.85)
        if(x == 28): return str(58.69)
        if(x == 27): return str(58.93)
        if(x == 29): return str(63.55)
        if(x == 30): return str(65.39)
        if(x == 31): return str(69.72)
        if(x == 32): return str(72.64)
        if(x == 33): return str(74.92)
        if(x == 34): return str(78.96)
        if(x == 35): return str(79.90)
        if(x == 36): return str(83.80)
        if(x == 37): return str(85.47)
        if(x == 38): return str(87.62)
        if(x == 39): return str(88.91)
        if(x == 40): return str(91.22)
        if(x == 41): return str(92.91)
        if(x == 42): return str(95.94)
        if(x == 43): return str(98.00)
        if(x == 44): return str(101.07)
        if(x == 45): return str(102.91)
        if(x == 46): return str(106.42)
        if(x == 47): return str(107.87)
        if(x == 48): return str(112.41)
        if(x == 49): return str(114.82)
        if(x == 50): return str(118.71)
        if(x == 51): return str(121.76)
        if(x == 53): return str(126.90)
        if(x == 52): return str(127.60)
        if(x == 54): return str(131.29)
        if(x == 55): return str(132.91)
        if(x == 56): return str(137.33)
        if(x == 57): return str(138.91)
        if(x == 58): return str(140.12)
        if(x == 59): return str(140.91)
        if(x == 60): return str(144.24)
        if(x == 61): return str(145.00)
        if(x == 62): return str(150.36)
        if(x == 63): return str(151.96)
        if(x == 64): return str(157.25)
        if(x == 65): return str(158.93)
        if(x == 66): return str(162.50)
        if(x == 67): return str(164.93)
        if(x == 68): return str(167.26)
        if(x == 69): return str(168.93)
        if(x == 70): return str(173.04)
        if(x == 71): return str(174.97)
        if(x == 72): return str(178.49)
        if(x == 73): return str(180.95)
        if(x == 74): return str(183.84)
        if(x == 75): return str(186.21)
        if(x == 76): return str(190.23)
        if(x == 77): return str(192.22)
        if(x == 78): return str(195.08)
        if(x == 79): return str(196.97)
        if(x == 80): return str(200.59)
        if(x == 81): return str(204.38)
        if(x == 82): return str(207.20)
        if(x == 83): return str(208.98)
        if(x == 84): return str(209.00)
        if(x == 85): return str(210.00)
        if(x == 86): return str(222.00)
        if(x == 87): return str(223.00)
        if(x == 88): return str(226.00)
        if(x == 89): return str(227.00)
        if(x == 91): return str(231.04)
        if(x == 90): return str(232.04)
        if(x == 93): return str(237.00)
        if(x == 92): return str(238.03)
        if(x == 95): return str(243.00)
        if(x == 94): return str(244.00)
        if(x == 96): return str(247.00)
        if(x == 97): return str(247.00)
        if(x == 98): return str(251.00)
        if(x == 99): return str(252.00)
        if(x == 100):return str(257.00)
        if(x == 101):return str(258.00)
        if(x == 102):return str(259.00)
        if(x == 103):return str(262.00)



if __name__ == "__main__":
    main();
