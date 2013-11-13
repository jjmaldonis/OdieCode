import os
import time

def main():
    dirs = [ name for name in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), name)) ]
    for dir in dirs:
        print dir
        lines = ['k_'+dir]
        for vk in os.listdir(dir + '/' + 'vk'):
            timestep = vk[vk.find('1'):-4]
            if timestep[0:12] == '10_npt_heat_':
                timestep = timestep[12:]
            vkf = open( dir + '/vk/vk_' + timestep + '.out', 'r' )
            line = vkf.readline()
            ln = 1
            while(len(line)!= 0):
                line = line.split(" ")
                while( '' in line ):
                    line.remove('')
                # line[1] is V(k) for this line, line[0] is k.
                # I am assuming all the ks are the same for each timestep.
                try:
                    lines[ln] = lines[ln] + ' ' + line[1][:-1]
                except:
                    lines.append(line[0] + ' ' + line[1][:-1])
                ln = ln + 1
                line = vkf.readline()
            vkf.close()
            lines[0] = lines[0] + ' ' + dir + '_' + timestep
        #print lines
        vkcf = open( dir + '/vk_compiled.out', 'w')
        for line in lines:
            vkcf.write(line + '\n')
        vkcf.close()


if __name__ == "__main__":
    main();

