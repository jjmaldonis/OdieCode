import os
import time

def main():
    for name in os.listdir(os.getcwd()):
        if name[len(name)-4:] == '.dat':
            parse_dat(name)



def parse_dat(name):
    df = open( name, 'r' )
    xyz = open(name[:-3] + 'xyz', 'w')
    line = df.readline() #Comment line
    xyz.write(line)
    line = df.readline() #Blank line
    natoms = int(df.readline().strip().split(" ")[0]) #natoms
    line = df.readline() #3 atom types
    xsize = sum( [ abs(float(x)) for x in df.readline().split(" ")[0:2] ] )
    ysize = sum( [ abs(float(x)) for x in df.readline().split(" ")[0:2] ] )
    zsize = sum( [ abs(float(x)) for x in df.readline().split(" ")[0:2] ] )
    print xsize, ysize, zsize
    xyz.write(str(xsize) + ' ' + str(ysize) + ' ' + str(zsize) + '\n')
    line = df.readline() #Blank line
    line = df.readline() #Blank line
    line = df.readline() #Atoms
    line = df.readline() #Blank line
    line = df.readline().strip() #First atom line

    while(len(line)!= 0):
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


if __name__ == "__main__":
    main();
