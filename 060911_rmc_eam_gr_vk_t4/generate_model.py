import random
import math

def main():
    natoms = 1425*8
    al_natoms = int(round(0.08 * natoms))
    cu_natoms = int(round(0.38 * natoms))
    zr_natoms = natoms - al_natoms - cu_natoms #Fill up whats left here
    #zr_natoms = 0.54
    print(str(al_natoms))
    print(str(cu_natoms))
    print(str(zr_natoms))
    print(str(natoms))

    cutoff = 1.76

    lx = 56.56900
    ly = 56.56900
    lz = 56.56900
    atoms = []
    f = open('double_model.xyz','w')
    f.write('double model\n')
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
