import os
import time

def main():
    xyzs = os.listdir('xyz_files/')
    for xyz in xyzs:
        os.system('qsub submit.sh ' + 'xyz_files/'+xyz + ' vor_results/'+xyz[:-4])
        #print(xyz)


if __name__ == "__main__":
    main();

