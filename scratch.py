import time
import pickle
import matplotlib.image as mpimg

start = time.time()
#with open('Mus_musculus.GRCm38.dna.chromosome.1.fa') as genome:
#    for i, line in enumerate(genome):
#        if '>' in line:
#            print("line {}: {}".format(i, line.rstrip()))
#print("Elapsed time: {}".format(time.time() - start))


def main():
    with open("a.pkl", mode="w+b") as p_out:
        thing = mpimg.imread('A_small.png')
        pickle.dump(thing, p_out)
if __name__ == "__main__":
    main()