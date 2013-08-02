import os
import time

def main():
    dirs = [ name for name in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), name)) ]
    for dir in dirs:
        print dir
        for xyz in os.listdir(dir + '/' + 'xyz_files'):
            run_femsim(xyz, dir)



def run_femsim(xyz, dir):
    timestep = xyz[xyz.find('1'):-4]
    if timestep[0:12] == '10_npt_heat_':
        timestep = timestep[12:]
    #print 'Running femsim on ' + dir+'/xyz_files/'+dir+'_npt_heat_'+timestep+'.xyz '
    if(timestep == '16500000'):
        os.system('qsub submit.sh "./rmc_for_md_sims ' + dir+'/xyz_files/'+dir+'_npt_heat_'+timestep+'.xyz ' + dir+'/param_file.in ' + dir+'/vk/vk_'+timestep+'.out"')
    
if __name__ == "__main__":
    main();

