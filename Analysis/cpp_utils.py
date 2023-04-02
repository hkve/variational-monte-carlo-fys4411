import subprocess
import pathlib as pl
import pandas as pd
import numpy as np

import os

def rootPath():
    cur_path = pl.Path(__file__)
    root_path = cur_path

    while root_path.name != "variational-monte-carlo-fys4411":
        root_path = root_path.parent

    return root_path

def vmcPath():
    vmc_path = rootPath() / pl.Path("build/vmc")
    return vmc_path

def gradientPath():
    gradient_path = rootPath() / pl.Path("build/gradient")
    return gradient_path

def timingPath():
    timing_path = rootPath() / pl.Path("build/timing")
    return timing_path

def interactPath():
    filename_interact = rootPath() / pl.Path(f"build/interact")
    return filename_interact

def parallelinteractPath():
    filename_interact = rootPath() / pl.Path(f"build/parallelinteract")
    return filename_interact

def dataPath(filename):
    filename_path = rootPath() / pl.Path(f"Data/{filename}")
    return filename_path


def vmcRun(D=3, N=10, logMet=20, logEq=16, omega=1.0, alpha=0.5, stepLength=0.1, importance=False, analytical=True, timing=False, GD=False, filename="test.txt"):
    vmc_path = vmcPath()
    filename_path = dataPath(filename)

    assert vmc_path.exists(), f"I cannot find {vmc_path} :((, are you sure you have compiled?"
    args = [
        vmc_path,
        D,
        N,
        logMet,
        logEq,
        omega,
        alpha,
        stepLength,
        int(importance),
        int(analytical),
        int(GD),
        filename_path,
    ]

    if not filename:
        args.pop()

    args_run = [str(arg) for arg in args]

    subprocess.run(args_run)

def parallelinteractRun(D=3, N=10, logMet=20, logEq=16, beta=1.0, alpha=0.5, stepLength=0.55, importance=False, analytical=True, GD=False, filename="test.txt", detailed=False):
    interact_path = parallelinteractPath()
    filename_path = dataPath(filename)
    #2^21 = 2097152

    assert interact_path.exists(), f"I cannot find {interact_path} :((, are you sure you have compiled?"
    args = [
        interact_path,
        D,
        N,
        logMet,
        logEq,
        alpha,
        beta,
        stepLength,
        int(importance),
        int(analytical),
        int(GD),
        filename_path,
        int(detailed)
    ]

    if not filename:
        args.pop()
    
    args_run = [str(arg) for arg in args]

    subprocess.run(args_run)


def gradientRun(D=3, N=10, logMet=16, logEq=14, alpha=0.5, stepLength=0.1, epsilon=0.01, lr= 0.01, importance=False, analytical=True, interacting=False, beta=1.0,  filename="gradientSearch.txt"):
    """
    This funcitons will run the gradient search for the best alpha parameter, as asked in the project description.
    Notice the regular VMC can run gradient descent, but this function will run the gradient search.
    """
    gradient_path = gradientPath()
    filename_path = dataPath(filename)

    assert gradient_path.exists(), f"I cannot find {gradient_path} :((, are you sure you have compiled?"
    args = [
        gradient_path,
        D,
        N,
        logMet,
        logEq,
        alpha,
        beta,
        stepLength,
        int(importance),
        int(analytical),
        lr,
        epsilon,
        int(interacting),
        filename_path,
    ]

    if not filename:
        args.pop()

    args_run = [str(arg) for arg in args]

    subprocess.run(args_run)


def timingRun(D=3, N=10, logMet=6, logEq=5, omega=1.0, alpha=0.5, stepLength=0.1, analytical=True, filename="timing.txt"):
    timing_path = timingPath()
    filename_path = dataPath(filename)

    assert timing_path.exists(), f"I cannot find {timing_path} :((, are you sure you have compiled?"
    args = [
        timing_path,
        D,
        N,
        logMet,
        logEq,
        omega,
        alpha,
        stepLength,
        int(analytical),
        filename_path,
    ]

    if not filename:
        args.pop()

    args_run = [str(arg) for arg in args]

    subprocess.run(args_run)

def interactRun(D=3, N=10, logMet=20, logEq=16, beta=1.0, alpha=0.5, stepLength=0.1, importance=False, analytical=True, timing=False, GD=False, filename="test.txt",detailed=False):
    interact_path = interactPath()
    filename_path = dataPath(filename)

    assert interact_path.exists(), f"I cannot find {interact_path} :((, are you sure you have compiled?"
    args = [
        interact_path,
        D,
        N,
        logMet,
        logEq,
        alpha,
        beta,
        stepLength,
        int(importance),
        int(analytical),
        int(GD),
        filename_path,
        int(detailed)
    ]

    if not filename:
        args.pop()

    args_run = [str(arg) for arg in args]

    subprocess.run(args_run)


def vmcLoad(filename):
    filename_path = dataPath(filename)

    df = pd.read_csv(filename_path, delim_whitespace=True)

    int_cols = ["Dimensions", "Particles" ,"Metro-steps", "Analytical", "Imposampling"]
    numeric_cols = [col for col in df.columns if col not in int_cols]
    
    for col in numeric_cols:
        df[col] = df[col].astype(float)

    return df

def gradientLoad(filename):
    filename_path = dataPath(filename)

    df = pd.read_csv(filename_path, delim_whitespace=True)

    int_cols = ["Dimensions", "Particles" ,"Metro-steps", "Analytical", "Imposampling"]
    numeric_cols = [col for col in df.columns if col not in int_cols]

    for col in numeric_cols:
        df[col] = df[col].astype(float)

    return df

def binaryLoad(filename):
    filename = dataPath(filename)
    with open(filename, "rb") as f:
        f.seek(0, 2)
        num_doubles = f.tell() // 8
        f.seek(0)

        return np.fromfile(f, dtype=np.float64, count=num_doubles)
    
def rLoad(filename):
    filename = dataPath(filename + "_Rs.txt")

    return pd.read_csv(filename, delim_whitespace=True)