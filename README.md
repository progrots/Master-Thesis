# Project description
This project entails simulation of an extended Klausmeier model, given by

	|u_t = u_{xx} + a - u - uvs
	|v_t = Dv_{xx} - mv + uvs(1-bv),	x\in\mathbb{R},
	|s_t = \frac{v-s}{\tau_s} 

and continuation of steady states in an approximation of the fast behaviour
of the extended Klausmeier model, given by				

	|u_t = u_{xx} + a - u - uvs_0  
	|v_t = Dv_{xx} - mv + uvs_0(1-bv).
___________________________________________________________________________________
## Simulation scripts:
	sim_Kplus_ode	(simulation of extended Klausmeier without diffusion terms)
	sim_Kplus_pde	(simulation of extended Klausmeier)
	sim_fast_pde	(simulation of approximation extended Klausmeier)

## Functions used in simulation:
	Klausmeier_plus_ode (extended Klausmeier without diffusion terms)
	Klausmeier_plus_pde (extended Klausmeier)
	Klausmeier_fast_pde (approximation extended Klausmeier)
	JacKlausmeier_plus_pde (jacobian of Klausmeier_plus_pde)
	JacKlausmeier_fast_pde (jacobian of Klausmeier_fast_pde)

## Continuation script:
	continuation_main

## Functions used in continuation:
	cont_Klausmeier.m, which uses
	findzeros.m, which uses
	newton.m

## Other scripts:
FindTuring (script that finds the location of saddle-node and Turing bifurcations in parameter a,
given the other model parameters).
