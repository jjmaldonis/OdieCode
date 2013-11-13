import os
import sys


def main():
    lines = []
    f = open(sys.argv[1],'r')
    line = f.readline() # comment line
    comment = line

    while( len(line) != 0 ):
        line = f.readline()
        line = line[:-1] + f.readline()
        lines.append(line)

    f.close()
    f = open(sys.argv[1]+'.new','w')
    f.write(comment)
    for line in lines:
        f.write(line)



if __name__ == '__main__':
    main()
