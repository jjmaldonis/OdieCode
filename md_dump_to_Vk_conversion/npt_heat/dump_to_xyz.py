import os
import time

def main():
    dirs = [ name for name in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), name)) ]
    for name in dirs:
        parse_dump(name)



def parse_dump(dirname):
    df = open( dirname + '/' + dirname + '_npt_heat.dump', 'r' )
    xyz = None
    line = df.readline()
    while(len(line)!= 0):
        # If we get to a new timestep create a new file and close the old one.
        if( line == "ITEM: TIMESTEP\n" ):
            if(xyz != None):
                xyz.write("-1")
                xyz.close()
            timestep = int(df.readline())
            line = df.readline() # ITEM: NUMBER OF ATOMS
            natoms = int(df.readline())
            line = df.readline() # ITEM: BOX BOUNDS pp pp pp
            xsize = sum( [ abs(float(x)) for x in df.readline().split(" ") ] )
            ysize = sum( [ abs(float(x)) for x in df.readline().split(" ") ] )
            zsize = sum( [ abs(float(x)) for x in df.readline().split(" ") ] )
            line = df.readline() # ITEM: ATOMS id type xs ys zs
            xyz = open( dirname + '/xyz_files/' + dirname + '_npt_heat_'+str(timestep)+'.xyz', 'w')
            print dirname + '/xyz_files/' + dirname + '_npt_heat_'+str(timestep)+'.xyz'
            xyz.write( str(natoms) + " atoms !Comment line\n")
            xyz.write( str(xsize) + " " + str(ysize) + " " + str(zsize) + " !Box size\n" )
            linenum = 9
        else:
            # We have an atom. Save it.
            # The types given in the dump file need to be converted to the actual elements' atomic numbers. In my case, they are all the same, so I will not code the universal way to convert them. You can find the masses for the types in the .dat files and then look up which mass is which element. For my case, 1 = 91.22 = Zr = 40, 2 = 63.55 = Cu = 29, and 3 = 26.98 = Al = 13.
            split_line = line.split(" ")[1:]
            #Reset atomic numbers
            if(split_line[0] == "1"):
                split_line[0] = "40"
            elif(split_line[0] == "2"):
                split_line[0] = "29"
            elif(split_line[0] == "3"):
                split_line[0] = "13"
            # Translate atomic positions
            split_line[1] = str( (float(split_line[1])-0.5)*xsize )
            split_line[2] = str( (float(split_line[2])-0.5)*ysize )
            split_line[3] = str( (float(split_line[3])-0.5)*zsize )
            # PBC on atomic positions
            if(float(split_line[1]) > xsize/2.0):
                split_line[1] = str(float(split_line[1]) - xsize)
            if(float(split_line[2]) > ysize/2.0):
                split_line[2] = str(float(split_line[2]) - ysize)
            if(float(split_line[3]) > zsize/2.0):
                split_line[3] = str(float(split_line[3]) - zsize)
            if(float(split_line[1]) < -xsize/2.0):
                split_line[1] = str(float(split_line[1]) + xsize)
            if(float(split_line[2]) < -ysize/2.0):
                split_line[2] = str(float(split_line[2]) + ysize)
            if(float(split_line[3]) < -zsize/2.0):
                split_line[3] = str(float(split_line[3]) + zsize)

            xyz.write( " ".join( split_line) + '\n')
        line = df.readline()
    xyz.write("-1")
    xyz.close()


if __name__ == "__main__":
    main();
