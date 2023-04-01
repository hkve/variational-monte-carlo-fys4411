# Variational Monte Carlo solver for FYS4411

Hello! Welcome to our project in Variational Monte Carlo!
In this project we have tried to simulate a Bose-Einstein system of particles in a magnetic trap, and find the lowest energy state for this system. To do this, we have used tools such as Variational Monte Carlo, Gradient Decent and the Blocking Method. This file will show how to navigate our repository and use our programs.

## The project report
Our project report can be found in Project1.pdf. This gives a further explenation about the project, our results and discussions.

## Compiling the project
To compile the different programs run `compile_project`. This will make excecutables for every program we have. It is recommended that you use a Mac for this.

## The different excecutables and how to run them
There are different excecutables built for different purposes. If you write `./<excecutable>` you will be given a list of different variables who needs a number. Write what number you want these variables to have in the order they are explained in the code. These will also be explained in each section below.
### vmc
The vmc excecutable works as a skelleton for the other codes. In itself it runs a Monte-Carlo simulation to find the lowest energy state for a given starting alpha.\
The excecutables in order are:\
#dims: How many dimensions do you want to the system to have and calculate in?\
#particles: How many particles do you want to  have in your system?\
#log2(metropolis steps): This is how many metropolis steps you want in your system. It will be the exponential of 2, so if this number is 3, the program will run $2^3=8$ metropolis steps.\
#log2(@-steps): This will determine the amount of Monte-Carlo steps the program will run before doing the Metropolis-steps. It will be the exponential of 2, so if this number is 3, the program will run $2^3=8$ equilibration steps before the Metropolis-steps begin.\
omega: This is the frequency of the trap. See the project report for more explenation.\
alpha: The starting alpha used to calculate the energy. See the project report for more explenation.\
stepLenght: Determines how far a particle is moved per Monte-Carlo cycle. \
Importantce sampling?: If this is set to 1, then the program will use importance sampling in the calculations. If it is set to 0, it will not.\
analytical?: If this is set to 1, the system will calculate the energy analytically. If it is set to 0, it will use a centered difference approximation to determine the energy. It is recommended to set this to 1 for shorter run time.\
gradientDescent?: If this is set to 1, the program will use the gradient descent to find the find the lowest alpha more efficiently.\
filename: The file the program will write out results to. If nothing is written, the program prints out the results in the terminal instead.


### interact
This excecutable works like vmc, but is also contains particle interaction in the model.\
The excecutables in order are:\
#dims: How many dimensions do you want to the system to have and calculate in?\
#particles: How many particles do you want to  have in your system?\
#log2(metropolis steps): This is how many metropolis steps you want in your system. It will be the exponential of 2, so if this number is 3, the program will run $2^3=8$ metropolis steps.\
#log2(@-steps): This will determine the amount of Monte-Carlo steps the program will run before doing the Metropolis-steps. It will be the exponential of 2, so if this number is 3, the program will run $2^3=8$ equilibration steps before the Metropolis-steps begin.\
omega: This is the frequency of the trap. See the project report for more explenation.\
alpha: The starting alpha used to calculate the energy. See the project report for more explenation.\
stepLenght: Determines how far a particle is moved per Monte-Carlo cycle. \
Importantce sampling?: If this is set to 1, then the program will use importance sampling in the calculations. If it is set to 0, it will not.\
analytical?: If this is set to 1, the system will calculate the energy analytically. If it is set to 0, it will use a centered difference approximation to determine the energy. It is recommended to set this to 1 for shorter run time.\
gradientDescent?: If this is set to 1, the program will use the gradient descent to find the find the lowest alpha more efficiently.\
filename: The file the program will write out results to. If nothing is written.\
detailed: If set to 1 the program will give the coordinates of the particles to a file and the information in the middle of the mc calculation to be used in the blocking.

### parallellinteract
This excecutable does the same as interact, but uses parallelization to speed things up a bit.
The excecutables in order are:
#dims: How many dimensions do you want to the system to have and calculate in?\
#particles: How many particles do you want to  have in your system?\
#log2(metropolis steps): This is how many metropolis steps you want in your system. It will be the exponential of 2, so if this number is 3, the program will run $2^3=8$ metropolis steps.\
#log2(@-steps): This will determine the amount of Monte-Carlo steps the program will run before doing the Metropolis-steps. It will be the exponential of 2, so if this number is 3, the program will run $2^3=8$ equilibration steps before the Metropolis-steps begin.\
omega: This is the frequency of the trap. See the project report for more explenation.\
alpha: The starting alpha used to calculate the energy. See the project report for more explenation.\
stepLenght: Determines how far a particle is moved per Monte-Carlo cycle. \
Importantce sampling?: If this is set to 1, then the program will use importance sampling in the calculations. If it is set to 0, it will not.\
analytical?: If this is set to 1, the system will calculate the energy analytically. If it is set to 0, it will use a centered difference approximation to determine the energy. It is recommended to set this to 1 for shorter run time.\
gradientDescent?: If this is set to 1, the program will use the gradient descent to find the find the lowest alpha more efficiently.\
filename: The file the program will write out results to. If nothing is written.\
detailed: If set to 1 the program will give the coordinates of the particles to a file and the information in the middle of the mc calculation to be used in the blocking.

### timing
The timing program is used for calculating the timing for the program under different conditions, such as different amount of particles or whether you use the analytical or numerical solver.\
The excecutables in order are:\
#dims: How many dimensions do you want to the system to have and calculate in?\
#particles: How many particles do you want to  have in your system?\
#log2(metropolis steps): This is how many metropolis steps you want in your system. It will be the exponential of 2, so if this number is 3, the program will run $2^3=8$ metropolis steps.\
#log2(@-steps): This will determine the amount of Monte-Carlo steps the program will run before doing the Metropolis-steps. It will be the exponential of 2, so if this number is 3, the program will run $2^3=8$ equilibration steps before the Metropolis-steps begin.\
omega: This is the frequency of the trap. See the project report for more explenation.\
alpha: The starting alpha used to calculate the energy. See the project report for more explenation.\
stepLenght: Determines how far a particle is moved per Monte-Carlo cycle. \
analytical?: If this is set to 1, the system will calculate the energy analytically. If it is set to 0, it will use a centered difference approximation to determine the energy. It is recommended to set this to 1 for shorter run time.\
filename: The file the program will write out results to. If nothing is written, the program prints out the results in the terminal instead.

### gradient
The gradient excecutable uses the gradient descent to calculate the best alpha for the system.\
The excutables in order are:\
#dims: How many dimensions do you want to the system to have and calculate in?\
#particles: How many particles do you want to  have in your system?\
#log2(metropolis steps): This is how many metropolis steps you want in your system. It will be the exponential of 2, so if this number is 3, the program will run $2^3=8$ metropolis steps.\
#log2(@-steps): This will determine the amount of Monte-Carlo steps the program will run before doing the Metropolis-steps. It will be the exponential of 2, so if this number is 3, the program will run $2^3=8$ equilibration steps before the Metropolis-steps begin.\
omega: This is the frequency of the trap. See the project report for more explenation.\
alpha: The starting alpha used to calculate the energy. See the project report for more explenation.\
stepLenght: Determines how far a particle is moved per Monte-Carlo cycle. \
Importantce sampling?: If this is set to 1, then the program will use importance sampling in the calculations. If it is set to 0, it will not.\
analytical?: If this is set to 1, the system will calculate the energy analytically. If it is set to 0, it will use a centered difference approximation to determine the energy. It is recommended to set this to 1 for shorter run time.\
learning rate: Decides a variable $l$ to be used in the calculation of the learning rate of the particles. The rest of the finction is given as $\frac{l}{ln(N)+1}$ hvere $N$ is the number of particles.\
epsilon: Decides the tolerance for the gradient descent.\
interaction?: If this is set to 1, the interaction gaussian will be used. If it is set to 0, it will not use interaction.\
filename: The file the program will write out results to. If nothing is written, the program prints out the results in the terminal instead.

## The python scripts
In the Analysis folder you will see some python scripts. Some of them are used to run the excecutables made in c++ and analyse the data in different ways.
### plot_derivative_timing.py
This is used to plot the timing excecutable and for different particles for both analytical and numerical calculation of the energy. It will then plot the time the caclualations took over the amount of particles calculated for both the analytical and numerical.


# From before. Should be deleted afterwards.
ah
## Compiling and running the project
The recommend way to compile this project is by using CMake to create a Makefile that you can then run. You can install CMake through one of the Linux package managers, e.g., `apt install cmake`, `pacman -S cmake`, etc. For Mac you can install using `brew install cmake`. Other ways of installing are shown here: [https://cmake.org/install/](https://cmake.org/install/).

### Compiling the project using CMake
In a Linux/Mac terminal this can be done by the following commands
```bash
# Create build-directory
mkdir build

# Move into the build-directory
cd build

# Run CMake to create a Makefile
cmake ../

# Make the Makefile using two threads
make -j2

# Move the executable to the top-directory
mv vmc ..
```
Or, simply run the script `compile_project` via
```bash
./compile_project
```
and the same set of commands are done for you. Now the project can be run by executing
```bash
./vmc
```
in the top-directory.

#### Cleaning the directory
Run `make clean` in the top-directory to remove the executable `vmc` and the `build`-directory.

#### Windows
Compilation of the project using Windows should work using CMake as it is OS-independent, but `make` does not work on Windows so the `compile_project`-script will not work.

## Completing the missing parts ##
Here follows a suggestion for how you can work to complete the missing parts of the code:
- Start by implementing the `SimpleGaussian` wave function: Write the `evaluate` function. Assume for now that the number of particles is always one, and the number of dimensions is always one. Next, compute the Laplacian analytically, and implement the `computeDoubleDerivative` function.
- Secondly, use the `Random` class (or your own favorite random number generator, should you have one) to implement the missing part of the `setupRandomUniformInitialState` function.
- Next, implement the metropolisStep function in the System class. Implement also the small missing part of the runMetropolisSteps function.
- Now, the last big thing needed is to implement the energy calculation. This is done by the `Hamiltonian` sub-class `HarmonicOscillator`. Here you will have to use the Laplacian you calculated for the wave function earlier.
- Now the code should be functioning and you should see (somewhat) reasonable results. Try to set the oscillator frequency to 1 and calculate analytically the energy of the oscillator. Recall the form of the ground state wave function of the harmonic oscillator, and set the parameter `alpha` accordingly. What is the resulting energy?
- If this energy is NOT correct, the last bit missing is to take a look at the `computeAverages` function in the `Sampler` class. What is missing here?

