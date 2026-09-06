from orbitalMotion import body
from orbitalMotion import simulation

def main():
    sim = simulation()
    sim.run("Beeman")
    sim.compareOrbitalPeriods()


if __name__ == "__main__":
    main()