import random
import math

#THIS IS THE CONFIGURATION USED TO CREATE PEI'S MODEL AFTER THE MD SIM PAUL RAN. THE CUTOFF IS TOO SMALL, IT SHOULD BE MORE LIKE 2.2940A BUT THAT DOESNT WORK BECAUSE RANDOM PLACEMENT OF THE ATOMS WONT LET THEM ALL IN. I THEN RAN AN MD SIM ON THE MODEL THIS FILE CREATED (WITH THE CUTOFF BELOW OF 2.0) THAT LOWERS THE MODEL CONFIGURATION TO ITS CLOSEST LOCAL MINIMUM (ENERGY WISE).

# This function generates a random model with natoms atoms. Below you can specify natoms, the composition, the box size, the cutoff distance, and the output filename.
# This function is completely random. No monte carlo or other type of algorithm is used. That means this may not be the best way to make a model if the cutoff distance, model size, and number of atoms doesn't leave much extra space available. A MC approach would be better.
def main():
    natoms = 1523
    al_natoms = int(round(0.08 * natoms))
    cu_natoms = int(round(0.38 * natoms))
    zr_natoms = natoms - al_natoms - cu_natoms #Fill up whats left here
    #zr_natoms = 0.54
    print(str(al_natoms))
    print(str(cu_natoms))
    print(str(zr_natoms))
    print(str(natoms))

    cutoff = 2.272

    lx = 20*math.sqrt(2)
    ly = 20*math.sqrt(2)
    lz = 20*math.sqrt(2)
    atoms = []
    f = open('Zr54Cu38Al8_1523atoms_post_md.xyz','w')
    f.write('Zr54Cu38Al8 1523atoms post md. atom density = ' + str(1000/math.pow(24.581,3)) + '\n')
    f.write( str(lx) + ' ' + str(ly) + ' ' + str(lz) + "\n")

    i = 0
    while(i<al_natoms):
        good = True
        x = random.uniform(-lx/2,lx/2)
        y = random.uniform(-ly/2,ly/2)
        z = random.uniform(-lz/2,lz/2)
        for atom in atoms:
            xx = atom[1]-x
            yy = atom[2]-y
            zz = atom[3]-z
            xx = xx-lx*int(round(xx/lx))
            yy = yy-ly*int(round(yy/ly))
            zz = zz-lz*int(round(zz/lz))
            r2 = xx**2 + yy**2 + zz**2
            r = math.sqrt(r2)
            if( r < cutoff ):
                good = False
        if(good):
            atoms.append( (13, x,y,z) )
            i+=1
        else:
            print("Stuck at "+ str(i))
            #print(atoms)

    while(i<al_natoms+cu_natoms):
        good = True
        x = random.uniform(-lx/2,lx/2)
        y = random.uniform(-ly/2,ly/2)
        z = random.uniform(-lz/2,lz/2)
        for atom in atoms:
            xx = atom[1]-x
            yy = atom[2]-y
            zz = atom[3]-z
            xx = xx-lx*int(round(xx/lx))
            yy = yy-ly*int(round(yy/ly))
            zz = zz-lz*int(round(zz/lz))
            r2 = xx**2 + yy**2 + zz**2
            r = math.sqrt(r2)
            if( r < cutoff ):
                good = False
        if(good):
            atoms.append( (29, x,y,z) )
            i+=1
        else:
            print("Stuck at "+ str(i))

    while(i<al_natoms+cu_natoms+zr_natoms):
        good = True
        x = random.uniform(-lx/2,lx/2)
        y = random.uniform(-ly/2,ly/2)
        z = random.uniform(-lz/2,lz/2)
        for atom in atoms:
            xx = atom[1]-x
            yy = atom[2]-y
            zz = atom[3]-z
            xx = xx-lx*int(round(xx/lx))
            yy = yy-ly*int(round(yy/ly))
            zz = zz-lz*int(round(zz/lz))
            r2 = xx**2 + yy**2 + zz**2
            r = math.sqrt(r2)
            if( r < cutoff ):
                good = False
        if(good):
            atoms.append( (40, x,y,z) )
            i+=1
        else:
            print("Stuck at "+ str(i))

    for atom in atoms:
        f.write( str(atom[0]) + " " + str(atom[1]) + " " + str(atom[2]) + " " + str(atom[3]) + "\n" )
    f.write("-1")

    f.close()


    

if __name__ == "__main__":
    main();
