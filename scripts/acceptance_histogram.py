
# This program will generate random numbers to test against the probabilistic acceptance of del_chi in the HRMC simulation. I will output all numbers that are accepted based on the random generation, and a histogram of those numbers will show the probability acceptance rate.


import math
import random
def main():
    temp = 200000
    beta = 1.0 / ( (8.6171e-5) * temp )
    f = open('acceptances.txt','w')
    i = 0
    while(i < 5000000):
        del_chi = random.random() * 1000
        rand = random.random()
        if( math.log(1 - rand ) < -1 * del_chi * beta ):
            # Accept the move.
            f.write(str(rand)+','+str(del_chi)+'\n')
            print i
            i = i + 1

if __name__ == "__main__":
    main();
