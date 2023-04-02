import numpy as np
import matplotlib.pyplot as plt
import plot_utils
import cpp_utils
import time 
import matplotlib as mpl

cmap = plot_utils.cmap 

def plot_energy_per_particle(filename="energy_per_particle", D=3, save=False):
    Ns = np.array([1, 10, 50, 100])
    alphas = np.linspace(0.1, 1.1, 11)

    stepLength = 0.55

    # start time measurement
    start = time.time()
    filename = f"{filename}_{D}D.txt"
    if not cpp_utils.dataPath(filename).exists():
        total = len(Ns)*len(alphas)
        for i, N in enumerate(Ns):
            for j, alpha in enumerate(alphas):
                cpp_utils.interactRun(D=D, filename=filename, N=N, alpha=alpha, stepLength=stepLength, logMet=18, logEq=14)
                print(f"Done N = {N}, alpha = {alpha} {i*len(alphas)+j+1}/{total}...")
                print(f"Time elapsed: {time.time() - start}")
                print("Expected time if linear: ", (time.time() - start)/(i*len(alphas)+j+1)*total / 60, "min")
 
    df = cpp_utils.vmcLoad(filename=filename)

    c = plot_utils.colors
    fig, ax = plt.subplots()
    for i, N in enumerate(Ns):
        df_N = df[ df.Particles == N ]
        E, E_std, alpha, MCCs = df_N["Energy"].to_numpy(), df_N["Energy_var"].to_numpy(), df_N["WF1"].to_numpy(), df_N["Metro-steps"].to_numpy()
        ax.errorbar(alpha, E/N, np.sqrt(E_std)/(MCCs), c=c[i], label=f"{N =}", marker="o")
        print(f"Minimum energy at alpha = {alphas[np.argmin(E)]}")

    ax.legend(ncol=4, bbox_to_anchor=(1.05, 1.15))
    ax.set(xlabel=r"$\alpha$", ylabel=r"$\langle E_L\rangle/N$")
    if save:
        plot_utils.save(filename.replace(".txt",f"_plot"))
    plt.show()

if __name__ == "__main__":
    plot_energy_per_particle(filename="energy_per_particle_interact", save=True)