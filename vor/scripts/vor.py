import sys
import os
import random
import string
import subprocess

class Vor(object):
    """ Basically this is a python script that I will use to call and store into memory the outputs from the voronoi.f90 code. It won't be pretty but it should work. """

    def __init__(self):
        """ constructor """
        super(Vor,self).__init__()


    def gen_paramfile(self,modelfile,cutoff):
        with open(modelfile) as f:
            content = f.readlines()
        content = [line.strip().split() for line in content]
        content.pop(0) # remove comment header
        content.pop(-1) # remove last '-1' line from model file
        # convert to ints and floats
        for j in range(0,len(content)):
            for i in range(0,len(line)):
                try:
                    content[j][i] = int(content[j][i])
                except:
                    try:
                        content[j][i] = float(content[j][i])
                    except:
                        pass
        box = content[0][0]
        content.pop(0) # remove world size line
        natoms = len(content)
        atom_types = []
        natom_types = {}
        for atom in content:
            if atom[0] not in atom_types:
                atom_types.append(atom[0])
                natom_types[atom[0]] = 1
            else:
                natom_types[atom[0]] += 1
        natom_types = [natom_types[i] for i in atom_types]
        for i in range(0,len(natom_types)):
            natom_types[i] = str(natom_types[i])
            atom_types[i] = str(atom_types[i])
        out = []
        out.append("# atom types")
        out.append(str(len(atom_types)) + ' ' + ' '.join(atom_types))
        out.append("# steps, # total atoms, # atoms of type1, # atoms of type2")
        out.append('1 '+str(natoms)+' '+" ".join(natom_types))
        out.append("# box size, # cut-off of neighbor")
        out.append(str(box)+' '+str(cutoff))
        return "\n".join(out)


    def run(self, modelfile, cutoff, outbase=None):
        """ Runs vor from command line. """

        if(type(cutoff) == type('hi')): 
            print(cutoff)
            cutoff = float(cutoff)

        if outbase is None:
            self.randstr = ''.join(random.choice(string.ascii_lowercase + string.ascii_uppercase + string.digits) for x in range(8))
        else:
            self.randstr = outbase

        opf = open(self.randstr+'.vparm','w')
        opf.write(self.gen_paramfile(modelfile,cutoff))
        opf.close()

        p = subprocess.Popen(['/export/home/jjmaldonis/OdieCode/vor/vor',self.randstr+'.vparm',modelfile,self.randstr], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.poutput = p.stdout.read()
        self.perr = p.stderr.read()
        self.preturncode = p.wait()
        if(self.preturncode != 0):
            self.del_files(self.randstr)
            raise Exception("Voronoi.f90 failed! "+self.perr)
        print("Voronoi.f90 exit status: "+str(self.preturncode))

    def save(self,outbase):
        """ Saves to memory the stat and index .out files with base 'outbase' """
        with open(outbase+'_index.out') as f:
            self.index = f.readlines()
        with open(outbase+'_stat.out') as f:
            self.stat = f.readlines()
        self.stat[0] = self.stat[0].strip() + self.stat[1]
        self.stat.pop(1)
        self.statheader = self.stat.pop(0) # Remove header too

    def del_files(self,outbase):
        """ Deletes the outbase files from the cwd. """
        os.remove(outbase+'_index.out')
        os.remove(outbase+'_stat.out')
        os.remove(self.randstr+'.vparm')

    def get_stat(self,):
        return self.stat
    def get_index(self,):
        return self.index

    def runall(self,modelfile,cutoff):
        self.run(modelfile,cutoff)
        self.save(self.randstr)
        self.del_files(self.randstr)


def main():
    vorrun = Vor()
    # sys.argv = [modelfile, cutoff, outbase (optional)]
    print(sys.argv)
    if(len(sys.argv) == 3):
        vorrun.run(sys.argv[1],sys.argv[2])
    elif(len(sys.argv) == 4):
        vorrun.run(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        sys.exit("Wrong number of inputs.")
    vorrun.save(vorrun.randstr)
    #for line in vorrun.index: print(line.strip())
    #print(vorrun.statheader.strip())
    #for line in vorrun.stat: print(line.strip())
    print(vorrun.poutput)
    print("Return code: "+str(vorrun.preturncode))
    print("Outbase: "+vorrun.randstr)
    #vorrun.del_files(vorrun.randstr)


if __name__ == "__main__":
    main()



