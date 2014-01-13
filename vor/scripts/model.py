
from atom import Atom
from hutch import Hutch
import sys
import math

class Model(object):
    """ xyz model file class """
    
    def __init__(self, *args, **kwargs):
        """ sets:
                self.comment
                self.lx
                self.ly
                self.lz
                self.atoms
                self.natoms """

        super(Model,self).__init__()
        if(len(args) == 1):
            modelfile = args[0]
            if modelfile[-4:] == '.xyz':
                self.read_xyz(modelfile)
            else:
                raise Exception("Unknown input model type!")
        elif(len(args) == 5):
            self.comment = args[0]
            self.lx = args[1]
            self.ly = args[2]
            self.lz = args[3]
            self.atoms = args[4]
            self.natoms = len(args[4])
        else:
            raise Exception("Unknown input parameters to Model()!")
        self.hutch = Hutch(self)


    def read_xyz(self,modelfile):
        with open(modelfile) as f:
            content = f.readlines()

        self.comment = content.pop(0) # Comment line
        content.pop(-1) # '-1' line
        self.lx,self.ly,self.lz = tuple([float(x) for x in content.pop(0).strip().split()])

        self.natoms = len(content)
        content = [x.strip().split() for x in content]
        for i in range(0,len(content)):
            for j in range(0,len(content[i])):
                try:
                    content[i][j] = int(content[i][j])
                except:
                    try:
                        content[i][j] = float(content[i][j])
                    except:
                        pass
        self.atoms = []
        for i,atom in enumerate(content):
            self.atoms.append(Atom(i,atom[0],atom[1],atom[2],atom[3]))

    def write_our_xyz(self,outfile):
        of = open(outfile,'w')
        of.write(self.comment)
        of.write("{0} {1} {2}\n".format(self.lx, self.lx, self.lz))
        for atom in self.atoms:
            of.write(atom.ourxyz()+'\n')
        of.write('-1')
        of.close()
        
    def add_atom(self,atom):
        self.atoms.append(atom)
        self.hutch.add_atom(atom)
        self.natoms += 1

    def generate_neighbors(self,cutoff):
        self.neighs = {}
        for atom in self.atoms:
            self.neighs[atom] = self.get_atoms_in_cutoff(atom,cutoff)
            #for i in range(0,len(self.neighs[atom])):
            #    # Neighbor syntax: (neighbor atom, dist to main atom)
            #    self.neighs[atom][i] = (self.neighs[atom][i],self.dist(atom,self.neighs[atom][i]))

    def get_atoms(self):
        return self.atoms
    def get_natoms(self):
        return self.natoms
    def get_box_size(self):
        return (self.lx,self.ly,self.lz)

    def get_atoms_in_cutoff(self,atom,cutoff):
        x_start = atom.coord[0] - cutoff
        y_start = atom.coord[1] - cutoff
        z_start = atom.coord[2] - cutoff
        x_end   = atom.coord[0] + cutoff
        y_end   = atom.coord[1] + cutoff
        z_end   = atom.coord[2] + cutoff
        if(x_start < -self.lx/2.0): x_start = x_start + self.lx #PBC
        if(y_start < -self.ly/2.0): y_start = y_start + self.ly #PBC
        if(z_start < -self.lz/2.0): z_start = z_start + self.lz #PBC
        if(x_end > self.lx/2.0): x_end = x_end - self.lx #PBC
        if(y_end > self.ly/2.0): y_end = y_end - self.ly #PBC
        if(z_end > self.lz/2.0): z_end = z_end - self.lz #PBC
        atom_s = Atom(0,0,x_start,y_start,z_start)
        hutch_s = self.hutch.hutch_position(atom_s)
        i_start = hutch_s[0]
        j_start = hutch_s[1]
        k_start = hutch_s[2]
        atom_e = Atom(0,0,x_end,y_end,z_end)
        hutch_e = self.hutch.hutch_position(atom_e)
        i_end = hutch_e[0]
        j_end = hutch_e[1]
        k_end = hutch_e[2]
        list = []

        #print(atom)
        #print(self.hutch.hutch_position(atom))
        #print(self.hutch.get_atoms_in_hutch(self.hutch.hutch_position(atom)))
        #print(atom)

        for i in range(0,self.hutch.nhutchs):
            if(i_start <= i_end):
                if(i < i_start or i > i_end): continue
            else:
                if(i < i_start and i > i_end): continue
            for j in range(0,self.hutch.nhutchs):
                if(j_start <= j_end):
                    if(j < j_start or j > j_end): continue
                else:
                    if(j < j_start and j > j_end): continue
                for k in range(0,self.hutch.nhutchs):
                    if(k_start <= k_end):
                        if(k < k_start or k > k_end): continue
                    else:
                        if(k < k_start and k > k_end): continue

                    list = list + self.hutch.get_atoms_in_hutch((i,j,k))
                    #print(self.hutch.get_atoms_in_hutch((i,j,k)))
                    #print((i,j,k))
        #print(atom)
        list.remove(atom)
        list = [atomi for atomi in list if self.dist(atom,atomi) <= cutoff]
        #for atomi in list:
        #    if(self.dist(atom,atomi) > cutoff):
        #        list.remove(atomi)
        #    else:
        #        print("dist={0}".format(self.dist(atom,atomi)))
        return list

    def dist(self,atom1,atom2):
        x = (atom1.coord[0] - atom2.coord[0])
        y = (atom1.coord[1] - atom2.coord[1])
        z = (atom1.coord[2] - atom2.coord[2])
        if(x > self.lx/2): x = self.lx - x
        if(y > self.lx/2): y = self.ly - y
        if(z > self.lx/2): z = self.lz - z
        if(x < -self.lx/2): x = self.lx + x
        if(y < -self.lx/2): y = self.ly + y
        if(z < -self.lx/2): z = self.lz + z
        x2 = x**2
        y2 = y**2
        z2 = z**2
        return math.sqrt(x2+y2+z2)
    def get_all_dists(self):
        dists = []
        for atomi in self.atoms:
            for atomj in self.atoms[self.atoms.index(atomi)+1:]:
                dists.append([atomi,atomj,self.dist(atomi,atomj)])
        return dists
    def nearest_neigh(self,atom):
        """ returns an atoms nearest neighbor """
        hutch = self.hutch.hutch_position(atom)
        atoms = self.hutch.get_atoms_in_hutch(hutch)[:]
        if atom in atoms: atoms.remove(atom)

        # This generation isn't perfect but it will work
        rots = [(1,0,0),(0,1,0),(0,0,1)]
        i = 0
        while len(atoms) == 0:
            hutch = ((hutch[0]+rots[i][0])%self.hutch.nhutchs,(hutch[1]+rots[i][1])%self.hutch.nhutchs,(hutch[2]+rots[i][2])%self.hutch.nhutchs)
            i = (i+1) % 3
            atoms = self.hutch.get_atoms_in_hutch(hutch)
            if atom in atoms: atoms.remove(atom)
        start = atoms[0]
        #print(atom)
        #print(start)
        #print(self.dist(atom,start))

        atoms = self.get_atoms_in_cutoff(atom,self.dist(atom,start))
        d = float("inf")
        for atomi in atoms:
            dt = self.dist(atom,atomi)
            if dt < d:
                d = dt
                a = atomi
        return a

    def save_vp_dict(self,vp_dict):
        self.vp_dict = vp_dict


def main():
    m = Model(sys.argv[1])
    dists = []
    for atom in m.atoms:
        dists.append(m.dist(atom,m.nearest_neigh(atom)))
    print(sum(dists)/len(dists))
    print(max(dists))
    print(min(dists))

if __name__ == "__main__":
    main()

