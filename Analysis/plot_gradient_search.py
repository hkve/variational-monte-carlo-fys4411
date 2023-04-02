import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plot_utils
import cpp_utils
import seaborn as sns
import time

import matplotlib as mpl
cmap = plot_utils.cmap 
from matplotlib.lines import Line2D


def plot_alpha_search(filename="gradientSearch", D=3, beta=1.0, alpha_range=(0.5, 0.5, 1), save=False, interacting=False):
    alphas = np.linspace(*alpha_range)
    Ns = [10, 50, 100]
    stepLengths = [0.55]#, 0.1, 0.5, 1.0]
    epsilon = 0.01
    lr= 0.002
    logMet = 20 # 2 ^ 18 = 262144
    logEq = 20 # 2 ^ 16 = 65536


    filename = f"{filename}_{D}D.txt" # prolly a good idea to add the dimension
    if not cpp_utils.dataPath(filename).exists(): # If the file does not exist, run the gradient search code
        total = len(Ns)*len(alphas)
        time_start = time.time()
        count = 0
        for N in Ns:
            for stepLength in stepLengths:
                for alpha in alphas:
                    cpp_utils.gradientRun(logMet=logMet, lr=lr ,logEq=logEq, D=D, N=N, alpha=alpha, stepLength=stepLength, epsilon=epsilon, filename=filename, beta=beta, interacting=interacting)
                    print(f"Gradient run done for alpha = {alpha} and N = {N}...")
                    count += 1
                    print(f"====>> Progress: {count}/{total} ({count/total*100:.2f}%)")
                    print(f"====>>====>>Time elapsed: {time.time() - time_start:.2f} s")
                    print(f"====>>====>>====>>Time remaining: {(time.time() - time_start)/count*(total-count)/60:.2f} min")


    info = f"_D={D}_N={Ns}_stepLength={stepLengths[0]}_met={logMet}_eq={logEq}_eps={epsilon}"

    df = cpp_utils.gradientLoad(filename=f"../Data/{filename}") # Load the data
    df_detailed = cpp_utils.gradientLoad(filename=f"../Data/detailed_{filename}")

    ax = sns.lineplot(data=df_detailed, x="Epoch", y="Alpha", hue="Particles", legend=True, palette=plot_utils.cmap)
    ax.get_legend().remove()
    plt.xlabel("Epoch")
    plt.ylabel(r"$\alpha$")
    # add legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, title="Particles")


    #norm = plt.Normalize(df_detailed["Alpha_0"].min(), df_detailed["Alpha_0"].max())
    #sm = plt.cm.ScalarMappable(cmap=plot_utils.cmap, norm=norm)
    #sm.set_array([])
    #plt.colorbar(sm, label=r"$\alpha_0$", cmap=plot_utils.cmap)

    if save:
        plot_utils.save("alpha_search" + info)
    plt.show()

    return df, df_detailed, info


def plot_energy_var(filename, df_detailed, info, save=False):
    # Plot the energy variance

    ax = sns.lineplot(data=df_detailed, x="Alpha", y=df_detailed["Energy_var"] / df_detailed["Particles"], hue="Particles", legend=True, palette=plot_utils.cmap)

    #norm = plt.Normalize(df_detailed["Alpha_0"].min(), df_detailed["Alpha_0"].max())
    #sm = plt.cm.ScalarMappable(cmap=plot_utils.cmap, norm=norm)
    #sm.set_array([])
    #plt.colorbar(sm, label=r"$\alpha_0$", cmap=plot_utils.cmap)

    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"Var$(\langle E_L \rangle)/N$")


    if save:
        plot_utils.save(filename + info)
    plt.show()


def plot_energy_per_particle(filename="GD_energy_per_particle", D=3, interacting=False, save=False):
    Ns = np.array([1, 10, 30, 50])
    alphas = [0.51]
    epsilon = 0.01 # this does not need to be super small. This is a tolerance with respect to the gradient, but notice 
                     # that this gets multiplied by the learning rate, so the paramaeter update at the end is what matters
    lr = 0.01 # this is scaled by the number of particles as epsilon = 1/(sqrt(N))
    stepLength = 1.2
    logMet = 2 # 2^19 = 524288
    logEq = 2 # 2 ^ 14 = 16384

    # start time measurement
    start = time.time()
    filename = f"{filename}_{D}D_interac{interacting}.txt"
    if not cpp_utils.dataPath(filename).exists():
        total = len(Ns)*len(alphas)
        for i, N in enumerate(Ns):
            for j, alpha in enumerate(alphas):
                cpp_utils.gradientRun(logMet=logMet, logEq=logEq, D=D, epsilon=epsilon, filename=filename, N=N, alpha=alpha, stepLength=stepLength, interacting=interacting, lr = lr)
                print(f"Done N = {N}, alpha = {alpha} {i*len(alphas)+j+1}/{total}...")
                print(f"Time elapsed: {time.time() - start}")
                print("Expected time if linear: ", (time.time() - start)/(i*len(alphas)+j+1)*total / 60, "min")

    df_detailed = cpp_utils.gradientLoad(filename=f"../Data/detailed_{filename}")

    c = plot_utils.colors
    fig, ax = plt.subplots()
    for i, N in enumerate(Ns):
        df_detailed_N = df_detailed[ df_detailed.Particles == N ]
        E, E_std, alpha, MCCs = df_detailed_N["Energy"].to_numpy(), df_detailed_N["Energy_var"].to_numpy(), df_detailed_N["WF1"].to_numpy(), df_detailed_N["Metro-steps"].to_numpy()
        ax.errorbar(alpha, E/N, np.sqrt(E_std)/(MCCs), c=c[i], label=f"{N =}")
        #print(f"Minimum energy at alpha = {alphas[np.argmin(E)]}")

    ax.legend(ncol=4, bbox_to_anchor=(1.05, 1.15))
    ax.set(xlabel=r"$\alpha$", ylabel=r"$\langle E_L \rangle [\hbar \omega]$")

    if save:
        plot_utils.save(filename.replace(".txt",f"_plot"))
    plt.show()

if __name__ == "__main__":
    #df, df_detailed, info = plot_alpha_search(filename="alpha_search_ho_int",D=3, save=True, beta=1, interacting=True)

    #plot_energy_var("energy_var_alpha_search_ho_int", df_detailed, info, save=True)


    # plot alpha search for ho and eo in same plot
    #df_ho, df_detailed_ho, info_ho = plot_alpha_search(filename="alpha_search_ho_int",D=3, save=True, beta=1, interacting=True)
    #df_eo, df_detailed_eo, info_eo = plot_alpha_search(filename="alpha_search_eo_int",D=3, save=True, beta=1, interacting=True)

    # plot energy variance for ho and eo in same plot
    #plot_energy_var("energy_var_alpha_search_ho_int", df_detailed_ho, info_ho, save=True)
    #plot_energy_var("energy_var_alpha_search_eo_int", df_detailed_eo, info_eo, save=True)

    # plot alpha search for ho and eo in same plot
    df_ho, df_detailed_ho, info_ho = plot_alpha_search(filename="alpha_search_ho_int",D=3, save=True, beta=1, interacting=True)
    df_eo, df_detailed_eo, info_eo = plot_alpha_search(filename="alpha_search_eo_int",D=3, save=True, beta=1, interacting=True)

    # plot alpha and epoch for ho and eo in same plot (with different colors)
    ax1 = sns.lineplot(data=df_detailed_ho, x="Epoch", y="Alpha", hue="Particles", legend=True, palette=plot_utils.cmap)

    # make this dashed
    ax2 = sns.lineplot(data=df_detailed_eo, x="Epoch", y="Alpha", hue="Particles", legend=False, palette=plot_utils.cmap, linestyle="--")

    # add legend
    handles, labels = ax1.get_legend_handles_labels()

    ax1.legend(handles=handles + [Line2D([0], [0], color="gray", linestyle='--')] + [Line2D([0], [0], color="gray", linestyle='-')], labels=labels + ["EO"] + ["HO"], title="Particles")
    # add to legend a continuous line as saying it is for the ho case without overriting the previous legend.

    # add labels
    plt.xlabel("Epoch")
    plt.ylabel(r"$\alpha$")

    # save and show
    plot_utils.save("alpha_and_epoch_ho_eo")
    plt.show()

    # do the same thig but with the energy variance per particle

    ax1 = sns.lineplot(data=df_detailed_ho, x="Epoch", y="Energy_var", hue="Particles", legend=True, palette=plot_utils.cmap)

    # make this dashed
    ax2 = sns.lineplot(data=df_detailed_eo, x="Epoch", y="Energy_var", hue="Particles", legend=False, palette=plot_utils.cmap, linestyle="--")

    # add legend
    handles, labels = ax1.get_legend_handles_labels()

    ax1.legend(handles=handles + [Line2D([0], [0], color="gray", linestyle='--')] + [Line2D([0], [0], color="gray", linestyle='-')], labels=labels + ["EO"] + ["HO"], title="Particles")
    # add to legend a continuous line as saying it is for the ho case without overriting the previous legend.

    # add labels
    plt.xlabel("Epoch")
    plt.ylabel(r"Log Var$(\langle E_L \rangle)/N$")
    plt.yscale("log")

    # save and show
    plot_utils.save("energy_var_and_epoch_ho_eo")
    plt.show()

    







