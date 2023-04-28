import numpy as np
from numba import njit
import sys


N     = 25
nu    = 5e-7
k0    = 0.05
lam   = 2.0
delta = 0.5
alpha = 1.5

t0 = 0.0
t = t0
dt = 1e-4
T  = 1000

forcing_shell = 0
forcing_scale = 0.1

rk1 = np.zeros((N + 4, ), dtype=np.complex128)
rk2 = np.zeros((N + 4, ), dtype=np.complex128)
rk3 = np.zeros((N + 4, ), dtype=np.complex128)
rk4 = np.zeros((N + 4, ), dtype=np.complex128)
u_tmp = np.zeros((N + 4, ), dtype=np.complex128)

a_rk1 = np.zeros((N + 4, ), dtype=np.float64)
a_rk2 = np.zeros((N + 4, ), dtype=np.float64)
a_rk3 = np.zeros((N + 4, ), dtype=np.float64)
a_rk4 = np.zeros((N + 4, ), dtype=np.float64)
a_tmp = np.zeros((N + 4, ), dtype=np.float64)
p_rk1 = np.zeros((N + 4, ), dtype=np.float64)
p_rk2 = np.zeros((N + 4, ), dtype=np.float64)
p_rk3 = np.zeros((N + 4, ), dtype=np.float64)
p_rk4 = np.zeros((N + 4, ), dtype=np.float64)
p_tmp = np.zeros((N + 4, ), dtype=np.float64)


Printing = False

save_every=1e4




##########################################################################################
######################### 					FUNCTION DEFS
##########################################################################################
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
def nonlin(non, u, k, N, delta, lam):
	for i in range(N):
		n = i + 2
		non[n] = 1j * k[n] * np.conjugate(u[n + 1] * u[n + 2] - (delta/lam) * (u[n - 1] * u[n + 1]) - ((1.0 - delta)/lam**2) * (u[n - 2] * u[n - 1]))
		if n == forcing_shell + 2:
			  non[n] += forcing_scale * (1.0 + 1j)
	return non

@njit
def nonlin_POAO(non_AO, non_PO, a, p, k, N, delta, lam, Mod_Type):
	for i in range(N):
		n = i + 2
		
		## Nonlinear Term
		if Mod_Type == "AMP_PHASE":
			non_AO[n] = k[n] * ((a[n + 1] * a[n + 2] * np.sin(p[n + 1] + p[n + 2] + p[n]))
			   			      - (a[n - 1] * a[n + 1] * np.sin(p[n - 1] + p[n + 1] + p[n])) * (delta/lam)
							  - (a[n - 2] * a[n - 1] * np.sin(p[n - 2] + p[n - 1] + p[n])) * ((1.0 - delta)/lam**2))

			non_PO[n] = (k[n] / a[n]) * ((a[n + 1] * a[n + 2] * np.cos(p[n + 1] + p[n + 2] + p[n]))
			            			   - (a[n - 1] * a[n + 1] * np.cos(p[n - 1] + p[n + 1] + p[n])) * (delta/lam) 
			 						   - (a[n - 2] * a[n - 1] * np.cos(p[n - 2] + p[n - 1] + p[n])) * ((1.0 - delta)/lam**2))
		elif Mod_Type == "AO":
			non_AO[n] = k[n] * ((a[n + 1] * a[n + 2] * np.sin(p[n + 1] + p[n + 2] + p[n]))
			   			      - (a[n - 1] * a[n + 1] * np.sin(p[n - 1] + p[n + 1] + p[n])) * (delta/lam)
							  - (a[n - 2] * a[n - 1] * np.sin(p[n - 2] + p[n - 1] + p[n])) * ((1.0 - delta)/lam**2))
			
			non_PO[n] = 0.0
		elif Mod_Type == "PO":
			non_AO[n] = 0.0
			
			non_PO[n] = (k[n] / a[n]) * ((a[n + 1] * a[n + 2] * np.cos(p[n + 1] + p[n + 2] + p[n]))
			            			   - (a[n - 1] * a[n + 1] * np.cos(p[n - 1] + p[n + 1] + p[n])) * (delta/lam) 
			 						   - (a[n - 2] * a[n - 1] * np.cos(p[n - 2] + p[n - 1] + p[n])) * ((1.0 - delta)/lam**2))
	
		## Add Forcing
		if n == forcing_shell + 2:
			if Mod_Type == "AMP_PHASE" or Mod_Type == "AO":
			  non_AO[n] += np.real(forcing_scale * (1.0 + 1j) * np.exp(-1j * p[n])) 		##(np.sqrt(2.0) * np.cos(np.pi/4.0 - p[n]))
			elif Mod_Type == "AMP_PHASE":
			  non_PO[n] += np.imag(forcing_scale * (1.0 + 1j) * np.exp(-1j * p[n])) / a[n]  ##(np.sqrt(2.0) * np.sin(np.pi/4.0 - p[n])) / a[n]

	return non_AO, non_PO
	
@njit
def get_exp(dt, nu, k_sqr):
	return np.exp(dt * -nu * k_sqr)






##-----------------------------------------
## START SOLVER
##-----------------------------------------
model_type  = "AO"


if model_type == "AO":
	solver_type = "INT_FAC"
elif model_type == "PO":
	solver_type = "RK4"


##-----------------------------------------
## Get k
##-----------------------------------------
k = init_k(N, k0, lam)
k_sqr = k**2

##-----------------------------------------
## Get IC
##-----------------------------------------
u, amp, phi = get_IC(k, N, "N_SCALING")
if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
	if Printing:
		for i in range(N + 4):
			print("k[{}]: {}\tu[{}]: {:1.16f} {:1.16f}\ta[{}]: {:1.16f}\tp[{}]: {:1.16f}".format(i, k[i], i, np.real(amp[i] * np.exp(1j * phi[i])), np.imag(amp[i] * np.exp(1j * phi[i])), i, amp[i], i, np.angle(np.exp(1j * phi[i]))))
		print("\n")
else:
	if Printing:
		for i in range(N + 4):
			print("k[{}]: {}\tu[{}]: {:1.16f} {:1.16f}\ta[{}]: {:1.16f}\tp[{}]: {:1.16f}".format(i, k[i], i, np.real(u[i]), np.imag(u[i]), i, np.absolute(u[i]), i, np.angle(u[i])))
		print("\n")



##-----------------------------------------
## Check Nonlinear Term
##-----------------------------------------
if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
	a_rk1, p_rk1 = nonlin_POAO(a_rk1, p_rk1, amp, phi, k, N, delta, lam, model_type)
	print(a_rk1, p_rk1)
	if Printing:
		for i in range(N + 4):
			print("k[{}]: {}\tRK_TEST[{}]: {:1.16f} {:1.16f}\ta[{}]: {:1.16f}\tp[{}]: {:1.16f}\tta: {:1.16f}\ttp: {:1.16f}".format(i, k[i], i, np.real(a_rk1[i]* np.exp(1j * p_rk1[i])), np.imag(a_rk1[i]* np.exp(1j * p_rk1[i])), i, np.absolute(a_rk1[i]* np.exp(1j * p_rk1[i])), i, np.angle(a_rk1[i]* np.exp(1j * p_rk1[i])), a_rk1[i], p_rk1[i]))
		print("\n")
else:
	rk1 = nonlin(rk1, u, k, N, delta, lam)
	if Printing:
		for i in range(N + 4):
			print("k[{}]: {}\tRK_TEST[{}]: {:1.16f} {:1.16f}\ta[{}]: {:1.16f}\tp[{}]: {:1.16f}".format(i, k[i], i, np.real(rk1[i]), np.imag(rk1[i]), i, np.absolute(rk1[i]), i, np.angle(rk1[i])))
		print("\n")






iters=1
while t <= T:

	##-----------------------------------------
	## Solver
	##-----------------------------------------
	if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
		u_tmp[:] = amp[:] * np.exp(1j * phi[:])
		a_tmp[:] = amp[:]
		p_tmp[:] = phi[:]
	else:
		u_tmp[:] = u[:]
	if Printing:
		for i in range(N + 4):
			print("in[{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}".format(i, np.real(u_tmp[i]), np.imag(u_tmp[i]), np.absolute(u_tmp[i]), np.angle(u_tmp[i])))
		print("\n\n")

	##-----------------------------------------
	## INTEGRATING FACTOR RK4
	##-----------------------------------------
	if solver_type == "INT_FAC":
			
		exp_dt  = get_exp(dt, nu, k_sqr)
		exp_dt2 = get_exp((dt/2.0), nu, k_sqr)

		##-------- STAGE 1
		if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
			a_rk1, p_rk1 = nonlin_POAO(a_rk1, p_rk1, a_tmp, p_tmp, k, N, delta, lam, model_type)
			a_tmp = exp_dt2 * amp + (dt/2.0) * exp_dt2 * a_rk1
			p_tmp = phi + (dt/2.0) * p_rk1
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("k[{}]: {}\tRK1[{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}\tRK_tmp[{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}".format(i, k[i], 
											i - 2, np.real(a_rk1[i]* np.exp(1j * p_rk1[i])), np.imag(a_rk1[i]* np.exp(1j * p_rk1[i])), np.absolute(a_rk1[i]* np.exp(1j * p_rk1[i])), np.angle(a_rk1[i]* np.exp(1j * p_rk1[i])),
											i - 2, np.real(a_tmp[i]* np.exp(1j * p_tmp[i])), np.imag(a_tmp[i]* np.exp(1j * p_tmp[i])), np.absolute(a_tmp[i]* np.exp(1j * p_tmp[i])), np.angle(a_tmp[i]* np.exp(1j * p_tmp[i]))))
				print("\n")
		else:
			rk1 = nonlin(rk1, u_tmp, k, N, delta, lam)
			u_tmp = exp_dt2 * u + (dt/2.0) * exp_dt2 * rk1
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("RK1[{}]: {:1.16f} {:1.16f}\t\tRK_tmp[{}]: {:1.16f} {:1.16f}".format(i - 2, np.real(rk1[i]), np.imag(rk1[i]), i - 2, np.real(u_tmp[i]), np.imag(u_tmp[i])))
				print("\n")

		##-------- STAGE 2
		if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
			a_rk2, p_rk2 = nonlin_POAO(a_rk2, p_rk2, a_tmp, p_tmp, k, N, delta, lam, model_type)
			a_tmp = exp_dt2 * amp + (dt/2.0) * a_rk2
			p_tmp = phi + (dt/2.0) * p_rk2
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("k[{}]: {}\tRK2[{}]: {:1.16f} {:1.16f}\ta[{}]: {:1.16f}\tp[{}]: {:1.16f}".format(i - 2, k[i], i - 2, np.real(a_rk2[i]* np.exp(1j * p_rk2[i])), np.imag(a_rk2[i]* np.exp(1j * p_rk2[i])), i - 2, a_rk2[i], i - 2, p_rk2[i]))
				print("\n")
		else:
			rk2    = nonlin(rk2, u_tmp, k, N, delta, lam)
			u_tmp = exp_dt2 * u + (dt/2.0) * rk2
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("RK2[{}]: {:1.16f} {:1.16f}\t\tRK_tmp[{}]: {:1.16f} {:1.16f}".format(i - 2, np.real(rk2[i]), np.imag(rk2[i]), i - 2, np.real(u_tmp[i]), np.imag(u_tmp[i])))
				print("\n")

		##-------- STAGE 3
		if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
			a_rk3, p_rk3 = nonlin_POAO(a_rk3, p_rk3, a_tmp, p_tmp, k, N, delta, lam, model_type)
			a_tmp = exp_dt * amp + dt * exp_dt2 * a_rk3
			p_tmp = phi + dt * p_rk3
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("k[{}]: {}\tRK3[{}]: {:1.16f} {:1.16f}\ta[{}]: {:1.16f}\tp[{}]: {:1.16f}".format(i - 2, k[i], i - 2, np.real(a_rk3[i]* np.exp(1j * p_rk3[i])), np.imag(a_rk3[i]* np.exp(1j * p_rk3[i])), i - 2, a_rk3[i], i - 2, p_rk3[i]))
				print("\n")
		else:
			rk3    = nonlin(rk3, u_tmp, k, N, delta, lam)
			u_tmp = exp_dt * u + dt * exp_dt2 * rk3
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("RK3[{}]: {:1.16f} {:1.16f}\t\tRK_tmp[{}]: {:1.16f} {:1.16f}".format(i - 2, np.real(rk3[i]), np.imag(rk3[i]), i - 2, np.real(u_tmp[i]), np.imag(u_tmp[i])))
				print("\n")

		##-------- STAGE 4
		if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
			a_rk4, p_rk4 = nonlin_POAO(a_rk4, p_rk4, a_tmp, p_tmp, k, N, delta, lam, model_type)
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("k[{}]: {}\tRK4[{}]: {:1.16f} {:1.16f}\ta[{}]: {:1.16f}\tp[{}]: {:1.16f}".format(i - 2, k[i], i - 2, np.real(a_rk4[i]* np.exp(1j * p_rk4[i])), np.imag(a_rk4[i]* np.exp(1j * p_rk4[i])), i - 2, a_rk4[i], i - 2, p_rk4[i]))
				print("\n")
		else:
			rk4    = nonlin(rk4, u_tmp, k, N, delta, lam)
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("RK4[{}]: {:1.16f} {:1.16f}".format(i - 2, np.real(rk4[i]), np.imag(rk4[i])))
				print("\n")

		##-------- Update
		if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
			amp = exp_dt * amp + dt * (1.0/6.0) * (exp_dt * a_rk1) + dt * (1.0/3.0) * (exp_dt2 * a_rk2) + dt * (1.0/3.0) * (exp_dt2 * a_rk3) + dt * (1.0/6.0) * a_rk4
			phi = phi + dt * (1.0/6.0) * p_rk1 + dt * (1.0/3.0) * p_rk2 + dt * (1.0/3.0) * p_rk3 + dt * (1.0/6.0) * p_rk4
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("INTFACunew[{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}".format(i - 2, np.real(amp[i] * np.exp(1j * phi[i])), np.imag(amp[i] * np.exp(1j * phi[i])), amp[i], phi[i]))
				print("\n")
		else:
			u = exp_dt * u + dt * (1.0/6.0) * (exp_dt * rk1) + dt * (1.0/3.0) * (exp_dt2 * rk2) + dt * (1.0/3.0) * (exp_dt2 * rk3) + dt * (1.0/6.0) * rk4
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("INTFACunew[{}]: {:1.16f} {:1.16f}".format(i - 2, np.real(u[i]), np.imag(u[i])))
				print("\n")
	

	##-----------------------------------------
	## RK4
	##-----------------------------------------
	elif solver_type == "RK4":
		##-------- STAGE 1
		if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
			a_rk1, p_rk1 = nonlin_POAO(a_rk1, p_rk1, a_tmp, p_tmp, k, N, delta, lam, model_type)
			if model_type == "AMP_PHASE" or model_type == "AO":
				a_rk1 -= nu * k * k * amp
			a_tmp = amp + (dt/2.0) * a_rk1
			p_tmp = phi + (dt/2.0) * p_rk1
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("k[{}]: {}\tRK1[{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}\tta: {:1.16f}\ttp: {:1.16f}\tRK_tmp[{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}".format(i, k[i], 
											i - 2, np.real(a_rk1[i]* np.exp(1j * p_rk1[i])), np.imag(a_rk1[i]* np.exp(1j * p_rk1[i])), np.absolute(a_rk1[i]* np.exp(1j * p_rk1[i])), np.angle(a_rk1[i]* np.exp(1j * p_rk1[i])), a_rk1[i], p_rk1[i],
											i - 2, np.real(a_tmp[i]* np.exp(1j * p_tmp[i])), np.imag(a_tmp[i]* np.exp(1j * p_tmp[i])), np.absolute(a_tmp[i]* np.exp(1j * p_tmp[i])), np.angle(a_tmp[i]* np.exp(1j * p_tmp[i]))))
				print("\n")
		else:
			rk1 = nonlin(rk1, u_tmp, k, N, delta, lam)
			rk1 -= nu * k * k * u 
			u_tmp = u + (dt/2.0) * (rk1)
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("k: {}\tnu: {:g}\tRK1[{}]: {:1.16f} {:1.16f}\t\tRK_tmp[{}]: {:1.16f} {:1.16f}".format(k[i], nu, i - 2, np.real(rk1[i]), np.imag(rk1[i]), i - 2, np.real(u_tmp[i]), np.imag(u_tmp[i])))
				print("\n")

		##-------- STAGE 2
		if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
			a_rk2, p_rk2 = nonlin_POAO(a_rk2, p_rk2, a_tmp, p_tmp, k, N, delta, lam, model_type)
			if model_type == "AMP_PHASE" or model_type == "AO":
				a_rk2 -= nu * k * k * amp
			a_tmp = amp + (dt/2.0) * a_rk2
			p_tmp = phi + (dt/2.0) * p_rk2
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("k[{}]: {}\tRK2[{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}\tta: {:1.16f}\ttp: {:1.16f}\tRK_tmp[{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}".format(i, k[i], 
											i - 2, np.real(a_rk2[i]* np.exp(1j * p_rk2[i])), np.imag(a_rk2[i]* np.exp(1j * p_rk2[i])), np.absolute(a_rk2[i]* np.exp(1j * p_rk2[i])), np.angle(a_rk2[i]* np.exp(1j * p_rk2[i])), a_rk2[i], p_rk2[i],
											i - 2, np.real(a_tmp[i]* np.exp(1j * p_tmp[i])), np.imag(a_tmp[i]* np.exp(1j * p_tmp[i])), np.absolute(a_tmp[i]* np.exp(1j * p_tmp[i])), np.angle(a_tmp[i]* np.exp(1j * p_tmp[i]))))
				print("\n")
		else:
			rk2    = nonlin(rk2, u_tmp, k, N, delta, lam)
			rk2 -= nu * k * k * u 
			u_tmp = u + (dt/2.0) * (rk2)
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("RK2[{}]: {:1.16f} {:1.16f}\t\tRK_tmp[{}]: {:1.16f} {:1.16f}".format(i - 2, np.real(rk2[i]), np.imag(rk2[i]), i - 2, np.real(u_tmp[i]), np.imag(u_tmp[i])))
				print("\n")

		##-------- STAGE 3
		if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
			a_rk3, p_rk3 = nonlin_POAO(a_rk3, p_rk3, a_tmp, p_tmp, k, N, delta, lam, model_type)
			if model_type == "AMP_PHASE" or model_type == "AO":
				a_rk3 -= nu * k * k * amp
			a_tmp = amp + dt * a_rk3
			p_tmp = phi + dt * p_rk3
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("k[{}]: {}\tRK3[{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}\tta: {:1.16f}\ttp: {:1.16f}\tRK_tmp[{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}".format(i, k[i], 
											i - 2, np.real(a_rk3[i]* np.exp(1j * p_rk3[i])), np.imag(a_rk3[i]* np.exp(1j * p_rk3[i])), np.absolute(a_rk3[i]* np.exp(1j * p_rk3[i])), np.angle(a_rk3[i]* np.exp(1j * p_rk3[i])), a_rk3[i], p_rk3[i],
											i - 2, np.real(a_tmp[i]* np.exp(1j * p_tmp[i])), np.imag(a_tmp[i]* np.exp(1j * p_tmp[i])), np.absolute(a_tmp[i]* np.exp(1j * p_tmp[i])), np.angle(a_tmp[i]* np.exp(1j * p_tmp[i]))))
				print("\n")
		else:
			rk3    = nonlin(rk3, u_tmp, k, N, delta, lam)
			rk3 -= nu * k * k * u 
			u_tmp = u + (dt) * (rk3)
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("RK3[{}]: {:1.16f} {:1.16f}\t\tRK_tmp[{}]: {:1.16f} {:1.16f}".format(i - 2, np.real(rk3[i]), np.imag(rk3[i]), i - 2, np.real(u_tmp[i]), np.imag(u_tmp[i])))
				print("\n")

		##-------- STAGE 4
		if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
			a_rk4, p_rk4 = nonlin_POAO(a_rk4, p_rk4, a_tmp, p_tmp, k, N, delta, lam, model_type)
			if model_type == "AMP_PHASE" or model_type == "AO":
				a_rk4 -= nu * k * k * amp
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("k[{}]: {}\tRK4[{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}\tta: {:1.16f}\ttp: {:1.16f}".format(i, k[i], 
											i - 2, np.real(a_rk4[i]* np.exp(1j * p_rk4[i])), np.imag(a_rk4[i]* np.exp(1j * p_rk4[i])), np.absolute(a_rk4[i]* np.exp(1j * p_rk4[i])), np.angle(a_rk4[i]* np.exp(1j * p_rk4[i])), a_rk4[i], p_rk4[i],))
				print("\n")
		else:
			rk4    = nonlin(rk4, u_tmp, k, N, delta, lam)
			rk4 -= nu * k * k * u
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("RK4[{}]: {:1.16f} {:1.16f}".format(i - 2, np.real(rk4[i]), np.imag(rk4[i])))
				print("\n")

		##-------- Update
		if model_type == "AMP_PHASE" or model_type == "AO" or model_type == "PO":
			amp = amp + dt * (1.0/6.0) * (a_rk1) + dt * (1.0/3.0) * (a_rk2) + dt * (1.0/3.0) * (a_rk3) + dt * (1.0/6.0) * a_rk4
			phi = phi + (dt * (1.0/6.0) * p_rk1 + dt * (1.0/3.0) * p_rk2 + dt * (1.0/3.0) * p_rk3 + dt * (1.0/6.0) * p_rk4)
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("RK4unew[{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}".format(i - 2, np.real(amp[i] * np.exp(1j * phi[i])), np.imag(amp[i] * np.exp(1j * phi[i])), amp[i], phi[i]))
				print("\n")
		else:
			u = u + dt * (1.0/6.0) * (rk1) + dt * (1.0/3.0) * rk2 + dt * (1.0/3.0) * rk3 + dt * (1.0/6.0) * rk4
			if Printing:
				for i in range(N + 4):
					if i >= 2 and i < N + 2:
						print("RK4unew[{}]: {:1.16f} {:1.16f}".format(i - 2, np.real(u[i]), np.imag(u[i])))
				print("\n")
	else:
		break

	##-----------------------------------------
	## Print Update to screen
	##-----------------------------------------
	if np.mod(iters, save_every) == 0.0:
		print("Iter: {}/{}\tEnrg: {:g}".format(iters, int(T/dt), np.sum(np.absolute(u)**2)*0.5))

	##-----------------------------------------
	## Update Time
	##-----------------------------------------
	iters+=1 
	t = iters*dt





























