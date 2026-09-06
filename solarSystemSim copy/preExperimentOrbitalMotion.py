
import os
import json
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation


class body():
    def __init__(self, name, mass, orbitalRad, colour, G, sunMass):
        self.name = name
        self.mass = mass
        self.colour = colour
        self.orbitalRad = orbitalRad
        self.G = G
        if orbitalRad == 0:
            self.position = np.zeros(2)
            self.velocity = np.zeros(2)
        else:
            self.position = np.array([orbitalRad,0])
            self.velocity = np.array([0,(self.G*sunMass/orbitalRad)**(1/2)])
        self.force = np.zeros(2)
        self.allAccelerations = [np.zeros(2)]
        self.allPositions = [self.position.copy()]
        self.orbitFound = False
        self.orbitPeriod = None

    

    
    def newForce(self, bodies):
        self.force = np.zeros(2) 
        for body in bodies:
            if body != self:
                relativePos = body.position - self.position
                relativeDistance = ((body.position[0] - self.position[0])**2 + (body.position[1] - self.position[1])**2)**(1/2)
                forceMagnitude= (self.G*body.mass*self.mass)/(relativeDistance)**2
                #possible divivision by zero if collision occurs
                forceDirection = relativePos/relativeDistance
                partForce = forceMagnitude * forceDirection
                self.force = self.force + partForce


    def calcAcceleration(self):
        self.allAccelerations.append(self.force/self.mass)
    
    def newPosition(self, timeStep):
        #make sure to trigger calc acceleration once, at the very start, manually
        self.position = self.position + self.velocity * timeStep + (1/6)*(4*self.allAccelerations[-1] - self.allAccelerations[-2])*timeStep**2
        
    def newVelocity(self, timeStep):
        self.velocity = self.velocity + (1/6)*(2*self.allAccelerations[-1]+5*self.allAccelerations[-2] - self.allAccelerations[-3])*timeStep

    def calcKineticEnergy(self):
        self.kineticEnergy = 1/2 * self.mass * ((self.velocity[0])**2 + (self.velocity[1])**2)
    

class simulation():
    def __init__(self, jsonFile = 'planetInfo.json'):
        with open(jsonFile) as f:
            parameters = json.load(f)
        self.timeStep = parameters["timestep"]
        self.runTime = parameters["runtime"]
        self.G = parameters["gravitational_constant"]
        sunMass = parameters["mass of sun"]
        self.currentTime = 0
        self.bodies = []
        for data in parameters["bodies"]:
            newBody = body(data["name"], data["mass"], data["orbital_radius"], data["colour"], self.G, sunMass)
            self.bodies.append(newBody)


    def step(self):
        self.currentPositions = []
        for body in self.bodies:
            body.newForce(self.bodies)
            body.calcAcceleration()
            body.newVelocity(self.timeStep)
            body.newPosition(self.timeStep)
            body.allPositions.append(body.position.copy())
        self.currentTime = self.currentTime + self.timeStep
    
    def separateCoords(self):
        for body in self.bodies:
            body.xCoords = np.array([position[0] for position in body.allPositions])
            body.yCoords = np.array([position[1] for position in body.allPositions])
               

    def detectOrbit(self):
        sun = None
        for body in self.bodies:
            if body.name == "Sun":
                sun = body
        if sun == None:
            print("no sun, no orbital period.")
        else:
            for body in self.bodies:
                if body.allPositions[-2][1] < sun.position[1] and body.allPositions[-1][1] >= sun.position[1] and body.orbitFound == False and body != sun:
                    body.orbitPeriod = self.currentTime
                    body.orbitFound = True
                    print(f"orbit period of {body.name} = {body.orbitPeriod/(365*24*60*60)} earth years")               
    
    def totalPotentialEnergy(self):
        self.potentialEnergy = 0
        for i in range(len(self.bodies)-1):
            for j in range(len(self.bodies)-1):
                if i != j:
                    relativePos = self.bodies[j].position - self.bodies[i].position
                    distance = ((relativePos[0])**2 + (relativePos[1])**2)**(1/2)
                    partPotentialEnergy = - (self.G * self.bodies[i].mass * self.bodies[j].mass)/distance
                    self.potentialEnergy = self.potentialEnergy + 1/2 * partPotentialEnergy

    def totalKineticEnergy(self):
        self.kineticEnergy = 0
        for body in self.bodies:
            body.calcKineticEnergy()
            self.kineticEnergy = self.kineticEnergy + body.kineticEnergy

    def energyFile(self, filename="energyData.json"):
        self.totalKineticEnergy()
        self.totalPotentialEnergy()
        totalEnergy = self.kineticEnergy + self.potentialEnergy
        energyData = ({ "time (s)": self.currentTime, "total_energy (j)": totalEnergy})

        with open(filename, 'a') as file: 
            json.dump(energyData, file)
            file.write("\n") 
        

    def run(self):
        if os.path.exists("energyData.json"):
            os.remove("energyData.json")
        for body in self.bodies:
            body.newForce(self.bodies)
            body.calcAcceleration()
        while self.currentTime <= self.runTime:
            self.step()
            self.detectOrbit()
            if self.currentTime % (5000*self.timeStep) == 0:
                self.energyFile()
        self.separateCoords()

    
    def animate(self):
        fig, ax = plt.subplots()
        scatters = ax.scatter([], [], s= 5, zorder=3,)  
        trails = {body.name: ax.plot([], [], label=body.name, linewidth=1, zorder=2)[0] for body in self.bodies}
        ax.set(xlim=[-5*10**12, 5*10**12], ylim=[-5*10**12, 5*10**12], xlabel='x (metres)', ylabel='y (metres)')
        #ax.set(xlim=[-5*10**11, 5*10**11], ylim=[-5*10**11, 5*10**11], xlabel='x (metres)', ylabel='y (metres)')
        ax.legend()

        def update(frame):
            x_data = []
            y_data = []
            for body in self.bodies:
                x_data.append(body.xCoords[frame])
                y_data.append(body.yCoords[frame])
                trails[body.name].set_data(body.xCoords[:frame+1], body.yCoords[:frame+1])
            scatters.set_offsets(np.column_stack((x_data, y_data)))
            return scatters

        ani = animation.FuncAnimation(fig, update, frames=len(self.bodies[0].xCoords), interval=5, repeat=False)
        plt.show()
        


sim = simulation()
sim.run()
sim.animate()

