def double_factorial(n):
    """computes the double factorial
    
    args
    ----
    n: integer greater than -1
    """
    assert n >= -1, "n >= -1"
    l = np.arange(1,n+1)
    if n == 0 or n == -1:
        return 1
    elif n%2 == 1:
        return np.prod(l[::2])
    else:
        return np.prod(l[1::2])