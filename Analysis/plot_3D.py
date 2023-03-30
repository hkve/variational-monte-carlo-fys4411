import numpy as np
import matplotlib.pyplot as plt
import cpp_utils

def main():
    df = cpp_utils.rLoad("test_pos")

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')    
    ax.scatter(df.x, df.y, df.z)

    l = 2
    ax.set(xlim=(-l,l), ylim=(-l,l), zlim=(-l,l))
    plt.show()
if __name__ == '__main__':
    main()