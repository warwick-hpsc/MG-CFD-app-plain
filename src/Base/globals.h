#ifndef GLOBALS_H
#define GLOBALS_H

extern config conf;

extern MeshName::MeshName mesh_name;

extern int levels;
extern int level;
extern int current_kernel;

extern double ff_variable[NVAR];
extern double3 ff_flux_contribution_momentum_x;
extern double3 ff_flux_contribution_momentum_y;
extern double3 ff_flux_contribution_momentum_z;
extern double3 ff_flux_contribution_density_energy;

#ifdef PAPI
extern int papi_monitoring_state;
#endif
#ifdef TIME
extern int timer_monitoring_state;
#endif

#endif
