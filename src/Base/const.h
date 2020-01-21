#ifndef CONST_H
#define CONST_H

/*
 * Options
 *
 */
#define GAMMA 1.4

#define NDIM 3

#define RK 3	// 3rd order RK
#define ff_mach 1.2
#define deg_angle_of_attack 0.0

/*
 * not options
 */
#define VAR_DENSITY 0
#define VAR_MOMENTUM 1
#define VAR_MOMENTUMX 1
#define VAR_MOMENTUMY 2
#define VAR_MOMENTUMZ 3
#define VAR_DENSITY_ENERGY (VAR_MOMENTUMX+NDIM)
 
#define NVAR (VAR_DENSITY_ENERGY+1)

#define MG_UP 0
#define MG_DOWN 1

#define COMPUTE_STEP 0
#define COMPUTE_FLUX_EDGE 1
#define UPDATE 2
#define INDIRECT_RW 3
#define TIME_STEP 4
#define RESTRICT 5
#define PROLONG 6
#define NUM_KERNELS (DOWN+1)

#define MESH_FVCORR 0
#define MESH_M6_WING 2
#define MESH_LA_CASCADE 3
#define MESH_ROTOR_37 4

#endif
