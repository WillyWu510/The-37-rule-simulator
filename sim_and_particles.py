import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import animation
from itertools import combinations

class Ball:
    def __init__(self, x, y, vx, vy,radius, score):
        """Initialize the position, velocity, radius, and the score of the ball."""
        self.position = np.array((x , y))
        self.velocity = np.array((vx,vy))
        self.radius=radius
        self.score=score
    #For convenience, map the component of the position and the velocity
    #onto their attribute x, y, vx, vy.
    @property
    def x(self):
        return self.position[0]
    @x.setter
    def x(self,value):
        self.r[0] = value
    @property
    def y(self):
        return self.position[1]
    @y.setter
    def y(self, value):
        self.position[1] = value
    @property
    def vx(self):
        return self.velocity[0]
    @vx.setter
    def vx(self, value):
        self.velocity[0] = value
    @property
    def vy(self):
        return self.velocity[1]
    @vy.setter
    def vy(self, value):
        self.velocity[1] = value   
    
    def overlaps(self, other):
        """check whether the particle overlaps with others."""
        return np.hypot(*(self.r - other.r)) < self.radius + other.radius
    def move(self, dt):
        """Move the Particle's position forward in time by dt."""
        
        self.position += self.velocity * dt
        # Make the Particles bounce off the walls
        # We set the size of box to be 10x10
        if self.x - self.radius < 0:
            self.x = self.radius
            self.vx = -self.vx
        if self.x + self.radius > 10:
            self.x = 1-self.radius
            self.vx = -self.vx
        if self.y - self.radius < 0:
            self.y = self.radius
            self.vy = -self.vy
        if self.y + self.radius > 10:
            self.y = 1-self.radius
            self.vy = -self.vy
class Simulation:
    """
       Simulation of two-dimensional elastic collision of particles.
       We set the size of box be 10x10 
        
    """
    def __init__(self, n,score):
        """
        Initialize the simulation of n particles.
        """
    def init_particles(self,n,score):
        """Initialize the position, velocity, radius, score of particles."""
    def collision(self):
        """Handle the collision of particles."""
    def animation(self):
        """Visualize the simulation"""