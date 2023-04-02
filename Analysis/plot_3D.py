import numpy as np
import matplotlib.pyplot as plt
import cpp_utils
import plot_utils
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

cmap = plot_utils.cmap 

def runSimulation(D=3, logMet=20, logEq=16, beta=1.0, alpha=0.5, stepLength=0.55, importance=False, analytical=True, GD=False, filename="test.txt", detailed=False):
    Ns = [10,50,100]

    if "eliptical" in filename:
        alphas = [0.49754, 0.48895, 0.48273] # these are the alphas found in the grid search
    elif "spherical" in filename:
        alphas = [0.49769, 0.48927, 0.48273]
    else:
        exit("Filename must contain either 'eliptical' or 'spherical'")


    if not cpp_utils.dataPath(filename).exists():
        total = len(Ns)
        k = 1
        for N, alpha in zip(Ns, alphas):
            if "Part" in filename:
                previous_N = int(filename.split("_")[1])
                filename = filename.replace(f"Part_{previous_N}_", f"Part_{N}_")
            else:
                filename = f"Part_{N}_" + filename
            cpp_utils.parallelinteractRun(D=D, N=N, logMet=logMet, logEq=logEq, beta=beta, alpha=alpha, stepLength=stepLength, importance=importance, analytical=analytical, GD=GD, filename=filename, detailed=detailed)
            print(f"Done {k}/{total}...")
            k += 1



def plot(filename="parallel_eliptical"):
    num_of_walkers = 10


    for N in [100]:
        if "Part" in filename:
            previous_N = int(filename.split("_")[1])
            filename = filename.replace(f"Part_{previous_N}_", f"Part_{N}_")
        else:
            filename = f"Part_{N}_" + filename

        df = pd.DataFrame()
        for i in range(num_of_walkers):
            # concatenate dataframes
            df_i = cpp_utils.rLoad(f"{filename}_{i}")
            df = pd.concat([df, df_i], ignore_index=True)

        plt.plot(df.x, df.y, df.z)

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')    
        ax.scatter(df.x, df.y, df.z)

        l = 2.5
        ax.set(xlim=(-l,l), ylim=(-l,l), zlim=(-l,l))
        # title
        potential = "eliptical" if "eliptical" in filename else "spherical"
        ax.set_title(f"{potential} potential of {N*num_of_walkers} interacting parts")
        plt.show()

        # make 3 histograms. One for each dimension of the position

        fig, axs = plt.subplots(3, 1)
        axs[0].hist(df.x, bins=20)
        axs[1].hist(df.y, bins=20)
        axs[2].hist(df.z, bins=20)

        # set all histogram x axis to be between -l and l
        l = 2.5
        for ax in axs:
            ax.set(xlim=(-l,l))

        # make a 3d color map volume plot that shows the density of the walkers
  
        # label the axes
        axs[0].set(xlabel='x', ylabel='Frequency')
        axs[1].set(xlabel='y', ylabel='Frequency')
        axs[2].set(xlabel='z', ylabel='Frequency')

        # space more evenly
        fig.tight_layout()
        plt.show()


def plot_overlay(filename1="parallel_eliptical", filename2="parallel_spherical", save=False):
    """
    for 100 particles, overlap the histrogram and 3d scatter plot of two different potentials
    """
    num_of_walkers = 10

    potential1 = filename1.removeprefix("parallel_")
    potential1 = potential1[0].upper() + potential1[1:]
    potential1 = potential1.replace("_", " ")


    potential2 = filename2.removeprefix("parallel_")
    potential2 = potential2[0].upper() + potential2[1:]
    potential2 = potential2.replace("_", " ")
   

    for N in [100]:
        if "Part" in filename1:
            previous_N = int(filename1.split("_")[1])
            filename1 = filename1.replace(f"Part_{previous_N}_", f"Part_{N}_")
        else:
            filename1 = f"Part_{N}_" + filename1
        
        if "Part" in filename2:
            previous_N = int(filename2.split("_")[1])
            filename2 = filename2.replace(f"Part_{previous_N}_", f"Part_{N}_")
        else:
            filename2 = f"Part_{N}_" + filename2


        df = pd.DataFrame()
        for i in range(num_of_walkers):
            df_i_1 = cpp_utils.rLoad(f"{filename1}_{i}")
            df_i_1["potential"] = potential1

            df_i_2 = cpp_utils.rLoad(f"{filename2}_{i}")
            df_i_2["potential"] = potential2

            # concatenate the two dataframes but saying from which potential they are from
            df_i = pd.concat([df_i_1, df_i_2], ignore_index=True)
            df = pd.concat([df, df_i], ignore_index=True)


        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        # legend the scatter plot

        types = np.unique(df.potential)

        c_array = np.array([0 if potential == types[0] else 1 for potential in df.potential])
        scatter = ax.scatter(df.x, df.y, df.z, c=c_array, cmap=cmap)
        plt.legend(handles=scatter.legend_elements()[0], labels=[potential1, potential2])

        info = f"overlay_N{N*num_of_walkers}_{potential1}_and_{potential2}_scatter"

        if save:
            plt.savefig(f"figs/{info}.pdf")
        plt.show()

        # show only the xz plane in a 2d scatter plot
        fig, ax = plt.subplots()
        scatter = ax.scatter(df.x, df.z, c=c_array, cmap=cmap)
    
        # get the x and z values for each potential
        x1 = df.x[df.potential == potential1]
        z1 = df.z[df.potential == potential1]
        x2 = df.x[df.potential == potential2]
        z2 = df.z[df.potential == potential2]

        # the standard deviation is the radius of the elipsis, the mean is the center
        # get the mean and standard deviation of the x and z positive values
        x1_mean = np.mean(x1)
        z1_mean = np.mean(z1)
        x2_mean = np.mean(x2)
        z2_mean = np.mean(z2)

        x1_std = np.std(x1)
        z1_std = np.std(z1)
        x2_std = np.std(x2)
        z2_std = np.std(z2)


        # get the x and z values for the elipsis
        x1_ellipse = np.linspace(x1_mean - 2*x1_std, x1_mean + 2*x1_std, 100) 
        z1_ellipse = np.linspace(z1_mean - 2*z1_std, z1_mean + 2*z1_std, 100) # not used because the elipsis is symmetric
        x2_ellipse = np.linspace(x2_mean - 2*x2_std, x2_mean + 2*x2_std, 100)
        z2_ellipse = np.linspace(z2_mean - 2*z2_std, z2_mean + 2*z2_std, 100) # not used because the elipsis is symmetric

        # draw the elipsis given these mean and standard deviations
        ax.plot(x1_ellipse, z1_mean + np.sqrt(1 - (x1_ellipse - x1_mean)**2/x1_std**2)*z1_std, color="black")
        ax.plot(x1_ellipse, z1_mean - np.sqrt(1 - (x1_ellipse - x1_mean)**2/x1_std**2)*z1_std, color="black")
        ax.plot(x2_ellipse, z2_mean + np.sqrt(1 - (x2_ellipse - x2_mean)**2/x2_std**2)*z2_std, color="orange")
        ax.plot(x2_ellipse, z2_mean - np.sqrt(1 - (x2_ellipse - x2_mean)**2/x2_std**2)*z2_std, color="orange")


        # legend the scatter plot and the elipsis in the same legend
        handles, labels = scatter.legend_elements()
        handles.append(Line2D([0], [0], color="black", lw=2))
        handles.append(Line2D([0], [0], color="orange", lw=2))

        print(labels, labels[0], labels[1])

        labels = []
        labels.append(potential1)
        labels.append(potential2)
        labels.append(potential1)
        labels.append(potential2)
        plt.legend(handles=handles, labels=labels)


        # label the axes
        ax.set(xlabel='x', ylabel='z')
        # set limits
        l = 2.5
        ax.set(xlim=(-l,l), ylim=(-l,l))
        
        info = f"overlay_N{N*num_of_walkers}_{potential1}_and_{potential2}_xz"
        if save:
            plt.savefig(f"figs/{info}.pdf")
        plt.show()


        # make 3 histograms. One for each dimension of the position. Overlay the two potentials
        fig, axs = plt.subplots(3, 1)
        axs[0].hist(df.x[df.potential == potential2], bins=30, label=potential2)
        axs[0].hist(df.x[df.potential == potential1], bins=30, label=potential1)
        axs[1].hist(df.y[df.potential == potential2], bins=30, label=potential2)
        axs[1].hist(df.y[df.potential == potential1], bins=30, label=potential1)
        axs[2].hist(df.z[df.potential == potential2], bins=30, label=potential2)
        axs[2].hist(df.z[df.potential == potential1], bins=30, label=potential1)
        

        # label the axes
        axs[0].set(xlabel='x', ylabel='')
        axs[1].set(xlabel='y', ylabel='Frequency')
        axs[2].set(xlabel='z', ylabel='')

        # space more evenly without tight layout
        fig.subplots_adjust(hspace=0.55)

        # set all histogram x axis to be between -l and l
        l = 2.5
        for ax in axs:
            ax.set(xlim=(-l,l))

        # use only one legend, outside the plot, on top
        ax.legend(ncol=2, bbox_to_anchor=(0.97, 4.7), loc='upper right', fontsize=14)

        info = f"overlay_N{N*num_of_walkers}_{potential1}_and_{potential2}_hist"
        if save:
            plt.savefig(f"figs/{info}.pdf")
        plt.show()



if __name__ == '__main__':

    #eliptical
    #runSimulation(logMet=22, logEq=18, beta=2.82843, stepLength=0.55, filename="parallel_eliptical", detailed=True)

    #eliptical no interaction !important: there is no interaction flag in the c++ code and i will not add it so this was hand-made
    # if this is run again, it will not work because the interactionTerm will not be 0
    #runSimulation(logMet=22, logEq=18, beta=2.82843, stepLength=0.55, filename="parallel_eliptical_no_interaction", detailed=True)


    #spherical
    #runSimulation(logMet=22, logEq=18, beta=1.0, stepLength=0.55, filename="parallel_spherical", detailed=True)
    #plot(filename="parallel_spherical")
    #plot(filename="parallel_eliptical")

    plot_overlay(filename1="parallel_eliptical", filename2="parallel_spherical", save=True)

    plot_overlay(filename1="parallel_eliptical", filename2="parallel_eliptical_no_interaction", save=True)



