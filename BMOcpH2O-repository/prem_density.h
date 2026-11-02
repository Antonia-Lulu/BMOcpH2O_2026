/* PREM mantle density table and BMO rho_nv(z) profile helpers. */

#ifndef PREM_DENSITY_H
#define PREM_DENSITY_H

struct All_variables;

/* Linear interpolation on density_PREM table [kg/m^3]. */
float prem_density_kg_m3(float depth_km);

/* Build BMO_rho_nv_prof on grid nodes; anchor at MO_depth_row grid depth.
 * Call after MO_depth_row is set (setup_plume_problem). */
void init_BMO_rho_nv_prof(struct All_variables *E);

/* rho_nv at dimensionless depth z_ctr (linear interp on prof or constant ref). */
float lookup_BMO_rho_nv(struct All_variables *E, float z_ctr);

/* Du et al. (2019) density from rho_nv and composition [kg/m^3]. */
float BMO_du2019_rho_kg_m3(float rho_nv, float XFe, float XH2O,
                           float n_spin, float alpha_w);

/* Du+2019 density [kg/m^3] at element-row row_k (same logic as compofuncbuoy). */
float BMO_rho_grid_row(struct All_variables *E, int row_k,
                              float *out_rho_nv, float *out_XFe, float *out_XH2O);

/* Depth-dependent liquidus pressure slope mP(z) and integral from z0_init. */
float BMO_get_mP_at_z(struct All_variables *E, float z);
float BMO_integral_mP_fromz0(struct All_variables *E, float z);

#endif /* PREM_DENSITY_H */
