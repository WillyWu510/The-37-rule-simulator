import numpy as np
import matplotlib.pyplot as plt 
def f(M,N):
    total = 0
    for i in range(M,N):
        total+=1/i
    return total

N=int(input())
x=np.arange(1,N+1)
constant1=np.array([1 for i in x])
f=np.array([f(i,N) for i in x])
idx = np.argwhere(np.diff(np.sign(f - constant1))).flatten()

plt.plot(x,f,'.',markeredgewidth = 0.01)
plt.plot(x,constant1,markeredgewidth = 0.01)
plt.plot(x[idx], f[idx], 'ro',label=f"M={x[idx][0]}")
plt.xlabel("M")
plt.ylabel(f"f(M,{N})")
plt.legend()
