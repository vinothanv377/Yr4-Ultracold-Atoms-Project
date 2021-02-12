import numpy as np

class Particle: 
    """A class repersenting the 3D partilce"""
    
    def __init__(self, x, y, z, vx, vy, vz, radius=0.01, styles=None):
        
        self.r = np.array((x, y, z))
        self.v = np.array((vx, vy, vz))
        self.radius = radius
    
    @property
    def x(self):
        return self.r[0]
    @x.setter
    def x(self, value):
        self.r[0] = value
    @property
    def y(self):
        return self.r[1]
    @y.setter
    def y(self, value):
        self.r[1] = value
    @property
    def z(self):
        return self.r[2]
    @z.setter
    def z(self, value):
        self.r[2] = value
    
    
    @property
    def vx(self):
        return self.v[0]
    @vx.setter
    def vx(self, value):
        self.v[0] = value
    @property
    def vy(self):
        return self.v[1]
    @vy.setter
    def vy(self, value):
        self.v[1] = value
    @property
    def vz(self):
        return self.v[2]
    @vz.setter
    def vz(self, value):
        self.v[2] = value
        
test_Particle = Particle(1, 2, 3, 4, 5, 6)
print(test_Particle.y)
test_Particle.y = 10
print(test_Particle.y)
print(test_Particle.r)

        
