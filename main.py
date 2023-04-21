import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    print("Hello")
    e1 = Engine(b=2, a=3)
    e1.add_stuff()
    print(e1.add_stuff)
    print(e1.c)
    print(e1.a)














class MaterialProperties():
    def __init__(self, YoungsModulus, PoissonsRatio):
        Eps = YoungsModulus
        Mu = PoissonsRatio
        E11 = Eps*(1-Mu)/((1+Mu)*(1-2*Mu))
        E12 = Eps*Mu/((1+Mu)*(1-2*Mu))
        E44 = Eps/(2*(1-Mu^2))
        self.D = np.array([[E11, E12, E12, 0, 0, 0],
                           [E12, E11, E12, 0, 0, 0],
                           [E12, E12, E11, 0, 0, 0],
                           [0, 0, 0, E44, 0, 0],
                           [0, 0, 0, 0, E44, 0],
                           [0, 0, 0, 0, 0, E44]])
        print(Eps); print(Mu)

if __name__ == '__main__':
    main()
