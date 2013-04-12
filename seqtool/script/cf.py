
from math import sqrt

def rpm(grav, radius=0.071):
    return sqrt((grav*9.8)/radius)*60/(2*3.14)

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(prog='rpm', description='calc rpm')
    parser.add_argument('gravity', nargs=1, type=float, help="gravity")
    parser.add_argument('radius', nargs='?', type=float, help="radius")
    args = parser.parse_args()

    gravity = args.gravity[0]
    radius = args.radius or 7.1

    print '{} rpm. where gravity is {}, radius is {}cm'.format(rpm(gravity, radius/100.), gravity, radius)

if __name__=='__main__':
    main()
