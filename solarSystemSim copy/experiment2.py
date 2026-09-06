from orbitalMotion import body
from orbitalMotion import simulation

def main():
    sim = simulation(runTime = 7500000000)
    sim.run("Beeman")
    sim.plotEnergyGraph("Beeman")
    sim.run("EC")
    sim.plotEnergyGraph("Beeman, EC")
    sim.run("DE")
    sim.plotEnergyGraph("Beeman, EC, DE")


if __name__ == "__main__":
    main()