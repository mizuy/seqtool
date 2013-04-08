
from math import sqrt

def rpm(grav, radius=0.071):
    return sqrt((grav*9.8)/radius)*60/(2*3.14)

print rpm(6000)
print rpm(13000)
