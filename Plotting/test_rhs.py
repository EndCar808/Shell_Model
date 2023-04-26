import numpy as np
from numba import njit


def init_k(N, k0, lmabda):
	k = np.zeros((N + 4, ), dtype=np.float64)
	for i in range(N + 4):
		if i >= 2 and  i < N + 2:
			k[i] = k0 * lmabda**(i - 2)
	return k


def get_IC(k, N, u0, turb_type):

	u = np.zeros((N + 4, ), dtype=np.complex128)

	if u0 == "N_SCALING":
		for i in range(N):
			n = i + 2
			amp   = 1.0 / (k[n]**alpha)
			phase = (i - 1)**2
			u[n]  = amp * np.exp(1j * phase)
	
	return u

@njit
def nonlin(non, u, k, N, delta, lam):
	for i in range(N):
		n = i + 2
		non[n] = 1j * k[n] * np.conjugate(u[n + 1] * u[n + 2] - (delta/lam) * (u[n - 1] * u[n + 1]) - ((1.0 - delta)/lam**2) * (u[n - 2] * u[n - 1]))
	

def IntFactRK4(non, u, k, k_sqr, N, dt, nu, delta, lam):

	exp_dt_m  = np.exp(-dt * nu * k_sqr)
	exp_dt2_m = np.exp(-(dt/2.0) * nu * k_sqr)
	exp_dt_p  = np.exp(dt * nu * k_sqr)
	exp_dt2_p = np.exp((dt/2.0) * nu * k_sqr)

	u_tmp = u

	rk1 = nonlin(non, u_tmp, k, N, delta, lam)
	
	u_tmp = exp_dt2_p * u + (dt/2.0) * exp_dt2_p * rk1
	rk2   = nonlin(non, u_tmp, k, N, delta, lam)
	
	u_tmp = exp_dt2_p * u + (dt/2.0) * rk2
	rk3   = nonlin(non, u_tmp, k, N, delta, lam)
	
	u_tmp = exp_dt_p * u + (dt) * exp_dt2_p * rk3
	rk4   = nonlin(non, u_tmp, k, N, delta, lam)


	u = exp_dt_p * u + dt * (1.0/6.0) * (exp_dt_m * rk1) + dt * (1.0/3.0) * (exp_dt2_m * rk2) + dt * (1.0/3.0) * (exp_dt2_m * rk3) + dt * (1.0/6.0) * rk4

N     = 25
nu    = 5e-7
k0    = 0.05
lam   = 2.0
delta = 0.5
alpha = 1.5

t0 = 0.0
dt = 1e-4
T  = 100


## Get k
k = init_k(N, k0, lam)

## Get IC
u = get_IC(k, N, "N_SCALING", "FULL")
for i in range(N + 4):
	print("k[{}]: {}\tu[{}]: {:1.16f} {:1.16f}".format(i, k[i], i, np.real(u[i]), np.imag(u[i])))


non = np.zeros((N + 4, ), dtype=np.complex128)

nonlin(non, u, k, N, delta, lam)