import numpy as np
import cpp_utils

def block(x):
    n = len(x)
    d = int(np.log2(n))
    s, gamma = np.zeros(d), np.zeros(d)
    mu = np.mean(x)

    for i in np.arange(0,d):
        n = len(x)
        gamma[i] = (n)**(-1)*np.sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        s[i] = np.var(x)
        x = 0.5*(x[0::2] + x[1::2])
   
    # generate the test observator M_k from the theorem
    M = (np.cumsum( ((gamma/s)**2*2**np.arange(1,d+1)[::-1])[::-1] )  )[::-1]

    q = np.array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in np.arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")

    return mu, s[k]/2**(d-k) # I think this is the /M

if __name__ == "__main__":
    """
    X = cpp_utils.binaryLoad("Part_10_parallel_eliptical_9_blocking_samples.dat")
    print(X)

    mu_X, error_X = np.mean(X), np.std(X)/np.sqrt(len(X))
    print("No blocking")
    print(f"{mu_X = }, {error_X = }")

    mu_block_X, var_block_X = block(X)
    err_block_X = np.sqrt(var_block_X)
    print("Blocking")
    print(f"{mu_block_X = }, err_block = {err_block_X}")
    """
    parts = [10, 50, 100]
    walkers = 8

    for part in parts:
        print(f"#### Part {part}")
        print(">>>>>blocking error of the eliptical with interaction")
        X = np.array([])
        for i in range(walkers):
            X = np.append(X, cpp_utils.binaryLoad(f"Part_{part}_parallel_eliptical_{i}_blocking_samples.dat"))

        mu_X, error_X = np.mean(X), np.std(X)/np.sqrt(len(X))
        #print(">>>>>>> No blocking")
        #print(f"{mu_X/part = }, {error_X/part = }")
        mu_block_X, var_block_X = block(X)
        err_block_X = np.sqrt(var_block_X)
        print(">>>>>>> Blocking")
        print(f"{mu_block_X/part = }, err_block/part = {err_block_X/part}")


        print(">>>>>blocking error of the eliptical without interaction")
        X = np.array([])
        for i in range(walkers):
            X = np.append(X, cpp_utils.binaryLoad(f"Part_{part}_parallel_eliptical_no_interaction_{i}_blocking_samples.dat"))


        mu_X, error_X = np.mean(X), np.std(X)/np.sqrt(len(X))
        #print(">>>>>>> No blocking")
        #print(f"{mu_X/part = }, {error_X/part = }")
    #
        mu_block_X, var_block_X = block(X)
        err_block_X = np.sqrt(var_block_X)
        print(">>>>>>> Blocking")
        print(f"{mu_block_X/part = }, err_block/part = {err_block_X/part}")

        print(">>>>> blocking error of the harmonic oscillator with interaction")

        X = np.array([])
        for i in range(walkers):
            X = np.append(X, cpp_utils.binaryLoad(f"Part_{part}_parallel_spherical_{i}_blocking_samples.dat"))


        mu_X, error_X = np.mean(X), np.std(X)/np.sqrt(len(X))
        #print(">>>>>>> No blocking")
        #print(f"{mu_X/part = }, {error_X/part = }")

        mu_block_X, var_block_X = block(X)
        err_block_X = np.sqrt(var_block_X)
        print(">>>>>>> Blocking")
        print(f"{mu_block_X/part = }, err_block/part = {err_block_X/part}")
