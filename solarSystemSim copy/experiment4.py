from orbitalMotion import body
from orbitalMotion import simulation

def main():
    #if runtime is increased, more alignments occur.
    sim = simulation()
    sim.run("Beeman")
    sim.printAlignments()


if __name__ == "__main__":
    main()