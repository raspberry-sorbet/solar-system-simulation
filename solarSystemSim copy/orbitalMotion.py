bu
import os
import json
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import math



class body():
    """
    Represents a celestial body in the simulation with attributes like mass,
    orbital radius, position, velocity, and force.
    Attributes:
    self.allAccelerations: list of numpy arrays representing all the body's acceleration values during the simulation.
    self.allPositions: list of numpy arrays representing all body's positions during the simulation.
    self.force: sum of all forces exerted onto body by other bodies in the system
    self.simOrbitPeriod: orbital period of body, calculated from the body's movement during the simulation.
    self.actualOrbitPeriod: true orbital period of body, taken from a pre-written file ("NASAorbitalPeriods.json")
    """
    def __init__(self, name, mass, orbitalRad, colour, G, sunMass, jsonFile ="NASAorbitalPeriods.json"):
        self.name = name
        self.mass = mass
        self.colour = colour
        self.orbitalRad = orbitalRad
        self.G = G
        #ensures that division by 0 does not occur
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
        self.simOrbitPeriod = None

        with open(jsonFile) as f:
            parameters = json.load(f) 
            #since the sun has no orbital period, the file does not have a value for it.
            self.actualOrbitPeriod = parameters.get(self.name, None)

    
    def newForce(self, bodies):
        """
        Calculates the total force acting on the
        body due to all other bodies in the system.
        """
        self.force = np.zeros(2) 
        for body in bodies:
            #ensures that the force between the body and itself is not calculated
            if body != self:
                relativePos = body.position - self.position
                relativeDistance = ((body.position[0] - self.position[0])**2 + (body.position[1] - self.position[1])**2)**(1/2)
                forceMagnitude= (self.G*body.mass*self.mass)/(relativeDistance)**2
                forceDirection = relativePos/relativeDistance
                partForce = forceMagnitude * forceDirection
                self.force = self.force + partForce

    def calcAcceleration(self):
        """updates the body's acceleration based on the total force acting on it."""
        self.allAccelerations.append(self.force/self.mass)
    
    def newPositionBeeman(self, timeStep):
        """
        updates the body's position using Beeman's algorithm, 
        based on the body's previous velocity, position and, both 
        previous and current accelerations
        """
        self.position = self.position + self.velocity * timeStep + (1/6)*(4*self.allAccelerations[-1] - self.allAccelerations[-2])*timeStep**2

    def newVelocityBeeman(self, timeStep):
        """
        updates the body's velocity using Beeman's algorithm,
        based on the body's previous velocity and, both previous 
        and current accelerations
        """  
        self.velocity = self.velocity + (1/6)*(2*self.allAccelerations[-1]+5*self.allAccelerations[-2] - self.allAccelerations[-3])*timeStep

    def newVelocityEC(self, timeStep):
        """
        updates the body's velocity using the Euler-Cromer method,
        based on the body's acceleration and previous velocity
        """
        self.velocity = self.velocity + self.allAccelerations[-1]*timeStep

    def newPositionEC(self, timeStep):
        """
        updates the body's position using the Euler-Cromer method, 
        based on the body's velocity and the previous position
        """
        self.position = self.position + self.velocity * timeStep

    def newPositionDE(self, timeStep):
        """
        updates the body's position using the Direct-Euler method, 
        based on the body's previous position and velocity
        """
        self.position = self.position + self.velocity*timeStep
    
    def newVelocityDE(self, timeStep):
        """
        updates the body's velocity using the Direct-Euler method, 
        based on the body's acceleration and previous velocity 
        """
        self.velocity = self.velocity + self.allAccelerations[-1]*timeStep

    def calcKineticEnergy(self):
        """calculates the body's kinetic energy"""
        return 1/2 * self.mass * ((self.velocity[0])**2 + (self.velocity[1])**2)
    
    

class simulation():
    """
    represents a system of celestial bodies interacting under
    gravitational forces, simulating their motion over time.
    Allows runTime to be changed if necessary.
    Attributes:
    self.runTime: the total duration for which the simulation will run
    self.G: gravitational constant
    self.timeStep: the time interval between each set of recalculations in the simulation.
    self.bodies: array of all bodies in the simulation
    self.allECEnergies: array of total energies of the simulation at each timestep, with the simulation having been run using the Euler-Cromer method.
    self.allBeemanEnergies: array of total energies of the simulation at each timestep, with the simulation having been run using Beeman's algorithm.
    self.allDEEnergies: array of total energies of the simulation at each timestep, with the simulation having been run using the Direct-Euler method.
    self.allTimeSteps: array of time values at each timestep during the simulation.
    self.alignedTime: array of timesteps at which a certain amount of planets aligned
    """
    def __init__(self, runTime = None, jsonFile = 'planetInfo.json'):
        with open(jsonFile) as f:
            parameters = json.load(f)
        self.timeStep = parameters["timestep"]
        if runTime == None:
            self.runTime = parameters["runtime"]
        else:
            self.runTime = runTime
        self.G = parameters["gravitational_constant"]
        sunMass = parameters["mass of sun"]
        self.currentTime = 0
        self.allECEnergies = []
        self.allBeemanEnergies = []
        self.allDEEnergies = []
        self.allTimeSteps = []
        self.alignedTimes = []
        self.bodies = []
        #creates an array of objects belonging to the body class
        for data in parameters["bodies"]:
            newBody = body(data["name"], data["mass"], data["orbital_radius"], data["colour"], self.G, sunMass)
            self.bodies.append(newBody)


    def stepBeeman(self):
        """
        updates the bodies' positions, velocities, and accelerations using 
        Beeman's algorithm, and appends the new positions to arrays
        """
        # since the values of newForce() is dependant on all bodies'
        # positions, using two different loops ensures that all bodies'
        # positions are updated before newForce() is called again.
        for body in self.bodies:
            body.newPositionBeeman(self.timeStep)
            body.allPositions.append(body.position.copy())
        for body in self.bodies:
            body.newForce(self.bodies)
            body.calcAcceleration()
            body.newVelocityBeeman(self.timeStep)
           
        

    def stepEC(self):
        """
        updates the bodies' positions, velocities, and accelerations using the 
        Euler-Cromer method, and appends the new positions to arrays
        """
        for body in self.bodies:
            body.newVelocityEC(self.timeStep)
            body.newPositionEC(self.timeStep)
            body.allPositions.append(body.position.copy())
            body.newForce(self.bodies)
            body.calcAcceleration()
            
    
    
    def stepDE(self):
        """
        updates the bodies' positions, velocities, and accelerations using the 
        Direct-Euler method, and appends the new positions to arrays
        """
        for body in self.bodies:
            body.newPositionDE(self.timeStep)
            body.newVelocityDE(self.timeStep)
            body.allPositions.append(body.position.copy())
            body.newForce(self.bodies)
            body.calcAcceleration()
            
   
            

    
    def separateCoords(self):
        """separates the body's position vectors into arrays of x and y coordinates."""
        for body in self.bodies:
            body.xCoords = np.array([position[0] for position in body.allPositions])
            body.yCoords = np.array([position[1] for position in body.allPositions])
    
    def sunDetector(self):
        """checks if a sun is present in the system"""
        for body in self.bodies:
            if body.name == "Sun":
                return body  
    
       

    def detectOrbit(self):
        """detects when a body completes its orbit around the sun, by comparing the Sun's and body's y coordinates."""
        sun = self.sunDetector()
        #checks if any bodies have passed the sun from below to above
        for body in self.bodies:
            if body != sun and body.orbitFound == False:
                if body.allPositions[-2][1] < sun.position[1] and body.allPositions[-1][1] >= sun.position[1]:
                    #converts orbital period from seconds into earth years
                    body.simOrbitPeriod = self.currentTime/((365*24*60*60))
                    body.orbitFound = True

    def printSimOrbits(self):
        """displays all simulated orbit periods"""
        for body in self.bodies:
            print(f"orbit period of {body.name} = {body.simOrbitPeriod} earth years")               
    
    def totalPotentialEnergy(self):
        """
        calculates potential energy of the system, based on the
        gravitational potential energy between each pair of bodies
        """
        self.potentialEnergy = 0
        for i in range(len(self.bodies)-1):
            for j in range(len(self.bodies)-1):
                #ensures that the force between the body and itself is not calculated
                if i != j:
                    relativePos = self.bodies[j].position - self.bodies[i].position
                    distance = ((relativePos[0])**2 + (relativePos[1])**2)**(1/2)
                    partPotentialEnergy = - (self.G * self.bodies[i].mass * self.bodies[j].mass)/distance
                    #force is halved, to account for all pairs being counted twice
                    self.potentialEnergy = self.potentialEnergy + 1/2 * partPotentialEnergy

    def totalKineticEnergy(self):
        """calculates the kinetic energy of the system"""
        self.kineticEnergy = 0
        for body in self.bodies:
            self.kineticEnergy = self.kineticEnergy + body.calcKineticEnergy()

    def storeEnergy(self, type, filename="energyData.json"):
        """
        Stores the total energy of the system in separate 
        arrays for Beeman's algorithm, Euler-Cromer, and Direct-Euler methods
        """
        self.totalKineticEnergy()
        self.totalPotentialEnergy()
        if type == "EC":
            self.allECEnergies.append(self.kineticEnergy + self.potentialEnergy)
            self.allTimeSteps.append(self.currentTime)
        elif type == "Beeman":
            self.allBeemanEnergies.append(self.kineticEnergy + self.potentialEnergy)
            self.allTimeSteps.append(self.currentTime)
            #writes the energy of the system to a file, every 5000 timesteps
            if self.currentTime % (5000*self.timeStep) == 0:
                energyData = ({ "time (s)": self.currentTime, "total_energy (j)": (self.kineticEnergy + self.potentialEnergy)})
                with open(filename, 'a') as file: 
                    json.dump(energyData, file)
                    file.write("\n")
        elif type == "DE":
            self.allDEEnergies.append(self.kineticEnergy + self.potentialEnergy)
            self.allTimeSteps.append(self.currentTime)
            

         
    def plotEnergyGraph(self, type):
        """
        Plots a graph of energy values over time, allowing
        comparison between energy conservation of Beeman's, 
        Euler-Cromer, and Direct-Euler methods
        """  
        plt.figure(figsize=(8, 6))
        plt.plot(self.allTimeSteps, self.allBeemanEnergies, label='Beeman', color='blue')
        #adds required plots if requested
        #ensures second plot is added unless only Beeman is requested
        if type != "Beeman":
            plt.plot(self.allTimeSteps, self.allECEnergies, label='Euler-Cromer', color='green')
        if type == "Beeman, EC, DE":
            plt.plot(self.allTimeSteps, self.allDEEnergies, label='Direct-Euler', color='red') 
            
        else:
            #sets appropriate limits for first two graphs
            plt.ylim(min(self.allBeemanEnergies)*0.999, max(self.allBeemanEnergies)*1.001)
            plt.xlim(0, self.allTimeSteps[-1])
        plt.xlabel('time (s)')
        plt.ylabel('total energy of system (J)') 
        plt.title("Energy Conservation Comparision")
        plt.legend()
        plt.grid()
        plt.show()


    def compareOrbitalPeriods(self):
        """"
        compares the simulated orbital period with the 
        actual orbital period for each planet in the system
        """
        for body in self.bodies:
            if body.name != "Sun":
                print(f"Simulated orbital period for {body.name}: {body.simOrbitPeriod} earth years")
                print(f"Actual orbital period for {body.name}: {body.actualOrbitPeriod} earth years")
                print(f"Difference: {abs(body.simOrbitPeriod - body.actualOrbitPeriod)} earth years \n")

    
    def calcAngle(self, body):
        """calculates and normalises a body's angle, in relation to the sun"""
        sun = self.sunDetector()
        relativePosition = body.position - sun.position
        planetAngle = math.atan2(relativePosition[1], relativePosition[0])
        if planetAngle < 0:
            planetAngle = planetAngle + 2*math.pi
        return planetAngle


    def calcMeanAngle(self):
        """calculates the mean angle of the first five planets in the system, treating the sun as the origin."""
        self.planets = []
        self.numCountedPlanets = 5
        #creates an array of only planets
        for body in self.bodies:
            if body.name != "Sun":
                self.planets.append(body)
        #sorts planets by orbital radius, to ensure that only innermost planets are used
        self.planets.sort(key=lambda body: body.orbitalRad)
        #ensures that enough planets exist
        if len(self.planets) >= self.numCountedPlanets:
            angleTotal = 0
            for i in range(self.numCountedPlanets):
                angleTotal = angleTotal + self.calcAngle(self.planets[i])
            return angleTotal/self.numCountedPlanets
 
    def detectAlignment(self):
        """
        detects when the first five planets fall within a 
        specified angular threshold of their mean angle
        """
        meanAngle = self.calcMeanAngle()
        alignmentCounter = 0
        alignmentGap = 7 * (math.pi/180)
        for i in range(self.numCountedPlanets):
            #checks if planet is within threshhold 
            if abs(self.calcAngle(self.planets[i]) - meanAngle) <= alignmentGap/2:
                alignmentCounter = alignmentCounter + 1
        #checks that all necessary planets have aligned
        if alignmentCounter == self.numCountedPlanets:
            self.alignedTimes.append(self.currentTime)


    def printAlignments(self):
        """displays times at which the planets aligned"""
        for time in self.alignedTimes:
            print(f"planetary alignment occured at {time} seconds")
        print(f"{len(self.alignedTimes)} planetary alignments, including the {self.numCountedPlanets} innermost planets, occured within {self.runTime} seconds")



    def run(self, type):
        """
        Runs the simulation by calling the appropriate step method 
        for the requested integration scheme at each timestep, 
        while also calling the methods for the experiments during the simulation loop. 
        This continues until the simulation duration reaches the specified length.

        Parameters:
        type: indicates which integration scheme is to be used.
        """
        self.allTimeSteps = []
        self.currentTime = 0
        #deletes file with energy data, if it already exists (if the simulation has already been run)
        if os.path.exists("energyData.json"):
            os.remove("energyData.json")
        for body in self.bodies:
            #ensures that previous acceleration values exist
            body.newForce(self.bodies)
            body.calcAcceleration()
        while self.currentTime <= self.runTime:
            if type == "Beeman":
                self.detectAlignment()
                self.stepBeeman()
            elif type == "EC":
                self.stepEC()
            elif type == "DE":
                self.stepDE()
            self.storeEnergy(type)
            self.detectOrbit()
            self.currentTime = self.currentTime + self.timeStep
        self.separateCoords()

    

    def animate(self):
        """
        Animates the motion of the planets in the simulation, 
        showing their positions and orbital trails over time.
        """
        fig, ax = plt.subplots(figsize=(8, 6))
        scatters = ax.scatter([], [], s= 5, zorder=3,)  
        trails = {body.name: ax.plot([], [], label=body.name, linewidth=1, zorder=2)[0] for body in self.bodies}
        plt.title("Planetary Motion Animation")
        ax.set(xlim=[-5*10**12, 5*10**12], ylim=[-5*10**12, 5*10**12], xlabel='x (metres)', ylabel='y (metres)')
        ax.legend()

        def update(frame):
            """
            updates the animation at each timestep by adjusting the scatter plot 
            to show the bodies' new positions and extending their orbital trails.       
            """
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
        

