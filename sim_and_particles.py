import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import animation
from itertools import combinations

class Particle:
    """A class representing a two-dimensional particle."""

    def __init__(self, x, y, vx, vy, radius=0.01, styles=None , score=0):
        """Initialize the particle's position, velocity, and radius.

        Any key-value pairs passed in the styles dictionary will be passed
        as arguments to Matplotlib's Circle patch constructor.

        """
        
        self.position = np.array((x, y))
        self.velocity = np.array((vx, vy))
        self.radius = radius
        self.score = score
        self.styles = styles
        if not self.styles:
            # Default circle styles
            self.styles = {'edgecolor': 'b', 'fill': False}
        self.particles_met=[]
        self.choice = -1
        self.expected_value = 0
    # For convenience, map the components of the particle's position and
    # velocity vector onto the attributes x, y, vx and vy.
    @property
    def x(self):
        return self.position[0]
    @x.setter
    def x(self, value):
        self.position[0] = value
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
        """Does the circle of this Particle overlap that of other?"""

        return np.hypot(*(self.position - other.position)) < self.radius + other.radius

    def draw(self, ax):
        """Add this Particle's Circle patch to the Matplotlib Axes ax."""

        circle = Circle(xy=self.position, radius=self.radius, **self.styles)
        ax.add_patch(circle)
        return circle

    def advance(self, dt):
        """Advance the Particle's position forward in time by dt."""

        self.position += self.velocity * dt

        # Make the Particles bounce off the walls
        if self.x - self.radius < 0:
            self.x = self.radius
            self.vx = -self.vx
        if self.x + self.radius > 10:
            self.x = 10-self.radius
            self.vx = -self.vx
        if self.y - self.radius < 0:
            self.y = self.radius
            self.vy = -self.vy
        if self.y + self.radius > 10:
            self.y = 10-self.radius
            self.vy = -self.vy
    def have_met(self,other):
        return other.score in self.particles_met
    
    def have_chosen(self):
        return self.choice != -1
    
class Simulation:
    """A class for a simple hard-circle molecular dynamics simulation.

    The simulation is carried out on a square domain: 0 <= x < 10, 0 <= y < 10.

    """

    def __init__(self, n, radius=0.01, styles=None, l_score = np.arange(1,1000), strategy = 0):
        """Initialize the simulation with n Particles with radii radius.

        radius can be a single value or a sequence with n values.

        Any key-value pairs passed in the styles dictionary will be passed
        as arguments to Matplotlib's Circle patch constructor when drawing
        the Particles.

        """
        self.fig, self.ax = plt.subplots()
        self.strategy = strategy
        self.init_particles(n, radius, styles , score=l_score)

    def init_particles(self, n, radius, styles=None , score = np.arange(1,1000)):
        """Initialize the n Particles of the simulation.

        Positions and velocities are chosen randomly; radius can be a single
        value or a sequence with n values.

        """
        try:
            iterator = iter(radius)
            assert n == len(radius)
        except TypeError:
            # r isn't iterable: turn it into a generator that returns the
            # same value n times.
            def r_gen(n, radius):
                for i in range(n):
                    yield radius
            radius = r_gen(n, radius)
        #np.random.shuffle(score)
        self.n = n+1
        self.particles = []
        domain_length=10
        
        x, y = 0.1 + (1- 2*0.1) * np.random.uniform(0,domain_length),0.1 + (1- 2*0.1) *   np.random.uniform(0,domain_length)
        vr = np.random.random() + 3.05
        vphi = 2*np.pi * np.random.random()
        vx, vy = vr * np.cos(vphi), vr * np.sin(vphi)
        self.You = Particle(x, y, vx, vy, 0.1, {'edgecolor': 'r','facecolor':'r', 'fill': True})
        self.particles.append(self.You)
        for i, rad in enumerate(radius):
            # Try to find a random initial position for this particle.
            while True:
                # Choose x, y so that the Particle is entirely inside the
                # domain of the simulation.
                
                x, y = rad + (1- 2*rad) * np.random.uniform(0,domain_length),rad + (1- 2*rad) * np.random.uniform(0,domain_length)
                # Choose a random velocity (within some reasonable range of
                # values) for the Particle.
                vr = np.random.random() + 3.05
                vphi = 2*np.pi * np.random.random()
                vx, vy = vr * np.cos(vphi), vr * np.sin(vphi)
                particle = Particle(x, y, vx, vy, rad, styles, score[i])
                # Check that the Particle doesn't overlap one that's already
                # been placed.
                if score[i] == n:
                    self.prince = particle
                    self.prince.styles = {'edgecolor': 'g','facecolor':'g', 'fill': True}
                for p2 in self.particles:
                    if p2.overlaps(particle):
                        break
                else:
                    self.particles.append(particle)
                    break
  
    def handle_collisions(self):
        """Detect and handle any collisions between the Particles.

        When two Particles collide, they do so elastically: their velocities
        change such that both energy and momentum are conserved.

        """

        def change_velocities(p1, p2):
            """
            Particles p1 and p2 have collided elastically: update their
            velocities.

            """

            m1, m2 = p1.radius**2, p2.radius**2
            M = m1 + m2
            r1, r2 = p1.position, p2.position
            d = np.linalg.norm(r1 - r2)**2
            v1, v2 = p1.velocity, p2.velocity
            u1 = v1 - 2*m2 / M * np.dot(v1-v2, r1-r2) / d * (r1 - r2)
            u2 = v2 - 2*m1 / M * np.dot(v2-v1, r2-r1) / d * (r2 - r1)
            p1.velocity = u1
            p2.velocity = u2

        # We're going to need a sequence of all of the pairs of particles when
        # we are detecting collisions. combinations generates pairs of indexes
        # into the self.particles list of Particles on the fly.
        pairs = combinations(range(self.n), 2)
        for i,j in pairs:
            if self.particles[i].overlaps(self.particles[j]):
                change_velocities(self.particles[i], self.particles[j])
                
                #Determin whether you have met the one before
                if self.You==self.particles[i]:
                    if not self.You.have_met(self.particles[j]) and len(self.You.particles_met) < self.strategy:
                        self.You.particles_met.append(self.particles[j].score)
                elif self.You==self.particles[j] :
                    if not self.You.have_met(self.particles[i]) and len(self.You.particles_met) < self.strategy:
                        self.You.particles_met.append(self.particles[i].score)
                if len(self.You.particles_met)!=0:
                    self.You.expected_value=max(self.You.particles_met)        
                #Choose the right one
                if self.You==self.particles[i] and len(self.You.particles_met) >= self.strategy and not self.You.have_chosen():
                    if self.particles[j].score > self.You.expected_value:
                        self.You.choice = self.particles[j].score
                elif self.You==self.particles[j] and len(self.You.particles_met) >= self.strategy and not self.You.have_chosen():
                    if self.particles[i].score > self.You.expected_value:
                        self.You.choice = self.particles[j].score
                    
    def advance_animation(self, dt):
        """Advance the animation by dt, returning the updated Circles list."""

        for i, p in enumerate(self.particles):
            p.advance(dt)
            self.circles[i].center = p.position
        self.handle_collisions()
        return self.circles

    def advance(self, dt):
        """Advance the animation by dt."""
        for i, p in enumerate(self.particles):
            p.advance(dt)
        self.handle_collisions()

    def init(self):
        """Initialize the Matplotlib animation."""

        self.circles = []
        for particle in self.particles:
            self.circles.append(particle.draw(self.ax))
        return self.circles

    
    def animate(self, i):
        """The function passed to Matplotlib's FuncAnimation routine."""

        self.advance_animation(0.1)
        if self.You.have_chosen():
            print(self.You.particles_met)
            print("I choose you!!!")
            plt.close(self.fig)
        elif len(self.You.particles_met)==self.n-1:
            print("No, I can't make a choice. I'm now a frog...Ribbit.")
            plt.close(self.fig)
        return self.circles

    def do_animation(self, save=False):
        """Set up and carry out the animation of the molecular dynamics.

        To save the animation as a MP4 movie, set save=True.
        """

        #fig, self.ax = plt.subplots()
        
        
        for s in ['top','bottom','left','right']:
            self.ax.spines[s].set_linewidth(2)
        self.ax.set_aspect('equal', 'box')
        self.ax.set_xlim(0, 10)
        self.ax.set_ylim(0, 10)
        self.ax.xaxis.set_ticks([])
        self.ax.yaxis.set_ticks([])
        anim = animation.FuncAnimation(self.fig, self.animate, init_func=self.init,
                               frames=800, interval=2, blit=True)
        
        
        if save:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=100, bitrate=1800)
            anim.save('collision.mp4', writer=writer)
        else:
            plt.show()
        if self.You.choice!=-1:
            print(self.You.choice)
            if self.You.choice == self.prince.score:
                print("Yes! You're my Mr.right.")
            else:
                print("Yuck, I kiss a real frog.")

if __name__ == '__main__':
    nparticles = 100
    strategy = int(input("please enter an integer:"))

    #radii = np.random.random(nparticles)*0.03+0.02
    radii=0.1
    styles = {'edgecolor': 'C0', 'linewidth': 2, 'fill': False}
    sim = Simulation(nparticles, radii, styles, np.arange(1,nparticles+1), strategy)
    sim.do_animation(save=False)
