import numpy as np
import matplotlib.pyplot as plt
import plot_utils
import cpp_utils

def plot_acceptance(filename="acceptance", stepLength_range=(0.1, 10, 10), save=False, impo=False):
    N = 10
    n = 100
    stepLengths = np.linspace(*stepLength_range)
    alphas = [0.2, 0.3, 0.4, 0.6, 0.7, 0.8]
    xlabel = {
        False: r"$\delta$",
        True: r"$\Delta t$",
    }
    extension = {
        False: "met",
        True: "impo"
    }

    filename = f"{filename}_{extension[impo]}"
    total = len(stepLengths)*len(alphas)
    if not cpp_utils.dataPath(filename).exists():
        for i, stepLength in enumerate(stepLengths):
            for j, alpha in enumerate(alphas):
                cpp_utils.vmcRun(D=3, N=N, logMet=15, logEq=10, filename=filename, stepLength=stepLength, alpha=alpha, importance=impo)
                print(f"{i*len(alphas)+j+1}/{total}")

    df = cpp_utils.vmcLoad(filename + ".txt")

    fig, ax = plt.subplots()
    for alpha in alphas:
        A = df[df["WF1"]==alpha]["Accept_ratio"]
        ax.plot(stepLengths, A, label=rf"$\alpha$ = {alpha}")

    ax.set(xlabel=xlabel[impo], ylabel="$A_r$")
    ax.legend()
    if save:
        plot_utils.save(filename + "_plot")
    plt.show()

def plot_compare_algos(filename="compare_algs"):
    N = 30
    Ms = np.linspace(14, 22, 100)

    if not cpp_utils.dataPath(filename + ".txt").exists():
        for i, M in enumerate(Ms):
            cpp_utils.vmcRun(D=3, N=N, logMet=M, logEq=10, filename=filename, stepLength=1.2, alpha=0.4, importance=False)
            cpp_utils.vmcRun(D=3, N=N, logMet=M, logEq=10, filename=filename, stepLength=0.8, alpha=0.4, importance=True)
            print(f"{i+1}/{len(Ms)}")

    df = cpp_utils.vmcLoad(filename + ".txt")

    df_m = df[df["Imposampling"] == 0]
    df_mh = df[df["Imposampling"] == 1]
    print(df_m["Accept_ratio"], df_mh["Accept_ratio"])
    
    fig, ax = plt.subplots()
    ax.plot(Ms, df_m.Energy, label="M")
    ax.plot(Ms, df_mh.Energy, label="MH")

    ax.legend()
    plt.show()
if __name__ == '__main__':
    # plot_acceptance(stepLength_range=(0.05, 5, 100), impo=False, save=True)
    # plot_acceptance(filename="test", stepLength_range=(0.05, 5, 100), impo=True, save=False)

    plot_compare_algos(filename="compare_algos_two")