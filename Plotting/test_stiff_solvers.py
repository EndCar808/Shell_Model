import numpy as np
from numba import njit
from rkstiff import grids
from rkstiff import if34, etd35, if4, etd4, etd5

#### ---------- System Parameters
N     = 22
nu    = 5e-7
k0    = 0.0625
lam   = 2.0
delta = 0.5
alpha = 1.5

t0 = 0.0
t = t0
dt = 1e-4
T  = 1000

forcing_shell = 1
forcing_scale = 0.005

Printing = False

save_every=1e5

#### ---------- Function Defs
def init_k(N, k0, lmabda):
    k = np.zeros((N + 4, ), dtype=np.float64)
    for i in range(N + 4):
        if i >= 2 and  i < N + 2:
            k[i] = k0 * lmabda**(i - 2)
    return k

def get_IC(k, N, u0):

    u = np.zeros((N + 4, ), dtype=np.complex128)
    amp = np.zeros((N + 4, ), dtype=np.float64)
    phi = np.zeros((N + 4, ), dtype=np.float64)

    if u0 == "N_SCALING":
        for i in range(N + 4):
            if i >= 2 and  i < N + 2:
                amp[i]   = 1.0 / (k[i]**alpha) / np.sqrt(75)
                phi[i] = (i - 1)**2
                u[i]  = amp[i] * np.exp(1j * phi[i])
    
    return u, amp, phi

@njit
def NL(u):
    non = np.zeros((u.shape), dtype=np.complex128)
    for i in range(u.shape[0] - 4):
        n = i + 2
        non[n] = 1j * k[n] * np.conjugate(u[n + 1] * u[n + 2] - (delta/lam) * (u[n - 1] * u[n + 1]) - ((1.0 - delta)/lam**2) * (u[n - 2] * u[n - 1]))
        if n == forcing_shell + 2:
              non[n] += forcing_scale * (1.0 + 1j)
    return non

@njit
def NL_AO(a):
    non = np.zeros((a.shape), dtype=np.float64)
    for i in range(a.shape[0] - 4):
        n = i + 2
        non[n] = k[n] * (a[n + 1] * a[n + 2] * np.sin(phi[n + 1] + phi[n + 2] + phi[n]) - (delta/lam) * (a[n - 1] * a[n + 1] * np.sin(phi[n - 1] + phi[n + 1] + phi[n])) - ((1.0 - delta)/lam**2) * (a[n - 2] * a[n - 1] * np.sin(phi[n - 1] + phi[n - 2] + phi[n])))
        if n == forcing_shell + 2:
              non[n] += np.real(forcing_scale * (1.0 + 1j) * np.exp(phi[n]))
    return non

def print_update(u, h, iters):
    if np.mod(iters, save_every) == 0.0:
        print("Iter: {}/{} t: {:g} \t Enrg: {:g}".format(iters, int(T/dt), h, np.sum(np.absolute(u)**2)*0.5))

def run_IF34(u0, dt, nl, l):
    solver = if34.IF34(linop=l,NLfunc=nl,epsilon=1e-8)

    h = dt
    u = u0.copy()
    iters = 1
    t = 0.0
    while t < T:
        u, h, h_suggest = solver.step(u, h)
        print_update(u, h, iters)
        t += h
        # use suggested step 
        h = h_suggest
        iters+=1

def run_ETD35(u0, dt, nl, l):
    solver = etd35.ETD35(linop=l,NLfunc=nl,epsilon=1e-8)

    h = dt
    u = u0.copy()
    iters = 1
    t = 0.0
    while t < T:
        u, h, h_suggest = solver.step(u, h)
        print_update(u, h, iters)
        t += h
        # use suggested step 
        h = h_suggest
        iters+=1


def run_IF4(u0, dt, nl, l):
    solver = if4.IF4(linop=l,NLfunc=nl)

    h = dt
    u = u0.copy()
    iters = 1
    t = 0.0
    while t < T:
        u = solver.step(u, h)
        print_update(u, h, iters)
        t += h
        # use suggested step 
        iters+=1

def run_ETD4(u0, dt, nl, l):
    solver = etd4.ETD4(linop=l,NLfunc=nl)

    h = dt
    u = u0.copy()
    iters = 1
    t = 0.0
    while t < T:
        u = solver.step(u, h)
        print_update(u, h, iters)
        t += h
        # use suggested step 
        iters+=1


def run_ETD5(u0, dt, nl, l):
    solver = etd5.ETD5(linop=l,NLfunc=nl)

    h = dt
    u = u0.copy()
    iters = 1
    t = 0.0
    while t < T:
        u = solver.step(u, h)
        print_update(u, h, iters)
        t += h
        # use suggested step 
        iters+=1



#### ---------- Setup  
k     = init_k(N, k0, lam)
k_sqr = k**2
L     = -nu * k_sqr

u0, amp, phi = get_IC(k, N, "N_SCALING")

####---------------- Test Full Model
# run_IF34(u0, dt, NL, L)

# run_ETD35(u0, dt, NL, L)

# run_IF4(u0, dt, NL, L)

# run_ETD4(u0, dt, NL, L)

# run_ETD5(u0, dt, NL, L)

####---------------- Test AO
# run_IF34(amp, dt, NL_AO, L)

# run_ETD35(amp, dt, NL_AO, L)

run_IF4(amp, dt, NL_AO, L)

# run_ETD4(amp, dt, NL_AO, L)

# run_ETD5(amp, dt, NL_AO, L)
