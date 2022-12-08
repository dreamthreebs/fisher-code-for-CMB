import numpy as np
import numba as nb
import time

def py_sum(a):
    Sum = 0
    for i in range(len(a)):
        Sum += a[i]
    return Sum

@nb.jit()
def nb_sum(a):
    Sum = 0
    for i in np.arange(len(a)):
        Sum += a[i]
    return Sum


#a = np.linspace(0,100,100)
a = np.linspace(0,100,10**6)
start=time.time()
np.sum(a)
end=time.time()
print('np.sum= ',end-start)
start=time.time()
nb_sum(a)
end=time.time()
print('nb.sum= ',end-start)


