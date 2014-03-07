import sys

if __name__ == '__main__':
    with open(sys.argv[1],'r') as f:
        lines = f.readlines()
    lines = [line.split() for line in lines]
    lines = [[line[0],line[1],line[2]] for line in lines]
    lines = [' '.join(line) for line in lines]

    for line in lines:
        print(line)

