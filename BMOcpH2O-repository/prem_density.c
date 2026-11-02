/* PREM mantle density vs depth, and BMO non-volatile density profile
 * rho_nv(z) = BMO_rho_solid_ref * PREM(z) / PREM(z_MO_grid), anchored at the MO grid row. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "prem_density.h"

static const float prem_depth_km_tab[] = {
    0.0f, 100.0f, 200.0f, 300.0f, 400.0f, 500.0f, 600.0f, 700.0f,
    800.0f, 900.0f, 1000.0f, 1100.0f, 1200.0f, 1300.0f, 1400.0f, 1500.0f,
    1600.0f, 1700.0f, 1800.0f, 1900.0f, 2000.0f, 2100.0f, 2200.0f, 2300.0f,
    2400.0f, 2500.0f, 2600.0f, 2700.0f, 2800.0f, 2891.0f
};

static const float prem_rho_kg_m3_tab[] = {
    3383.0f, 3373.0f, 3362.0f, 3484.0f, 3543.0f, 3850.0f, 3976.0f, 4399.0f,
    4461.0f, 4521.0f, 4580.0f, 4638.0f, 4695.0f, 4751.0f, 4806.0f, 4860.0f,
    4913.0f, 4966.0f, 5018.0f, 5070.0f, 5121.0f, 5171.0f, 5222.0f, 5272.0f,
    5322.0f, 5372.0f, 5421.0f, 5471.0f, 5521.0f, 5566.0f
};

static const int prem_ntab = 30;

float prem_density_kg_m3(float depth_km)
{
    int i;
    float d, t;

    if (depth_km <= prem_depth_km_tab[0])
        return prem_rho_kg_m3_tab[0];
    if (depth_km >= prem_depth_km_tab[prem_ntab - 1])
        return prem_rho_kg_m3_tab[prem_ntab - 1];

    for (i = 0; i < prem_ntab - 1; i++) {
        if (depth_km >= prem_depth_km_tab[i] && depth_km <= prem_depth_km_tab[i + 1]) {
            d = prem_depth_km_tab[i + 1] - prem_depth_km_tab[i];
            if (fabs(d) < 1.0e-12f)
                return prem_rho_kg_m3_tab[i];
            t = (depth_km - prem_depth_km_tab[i]) / d;
            return prem_rho_kg_m3_tab[i] * (1.0f - t)
                 + prem_rho_kg_m3_tab[i + 1] * t;
        }
    }
    return prem_rho_kg_m3_tab[prem_ntab - 1];
}

float BMO_du2019_rho_kg_m3(float rho_nv, float XFe, float XH2O,
                           float n_spin, float alpha_w)
{
    float coeff = 0.28f + 0.04f * n_spin;
    float rho_anh = rho_nv * (1.0f + XFe * coeff);
    /* XH2O is mass fraction; Du+2019 water term uses alpha_w * XH2O directly. */
    return rho_anh * (1.0f + alpha_w * XH2O);
}

void init_BMO_rho_nv_prof(E)
    struct All_variables *E;
{
    int k, k0, noz, lev;
    float depth_km, rho_prem_k, rho_ref, z_MO_km, z_MO_nom_km, rho_prem_MO;
    float coeff, rho_nv_MO, rho_upper_MO, rho_melt_MO;
    float XFe_upper, XH2O_upper, XFe_melt, XH2O_melt;

    if (E->control.BMO_density_mode != 1 || E->control.BMO_rho_nv_prem == 0)
        return;

    if (E->control.BMO_rho_nv_prof != NULL)
        return;

    lev = E->mesh.levmax;
    noz = E->lmesh.NOZ[lev];

    /* Anchor at the actual MO grid row, not the nominal input MO_depth [km].
     * PREM at the anchor is used only as a normaliser for relative scaling. */
    k0 = E->control.MO_depth_row;
    if (k0 < 1)     k0 = 1;
    if (k0 > noz)   k0 = noz;
    z_MO_nom_km = E->control.MO_depth_initial * E->data.layer_km;
    z_MO_km = E->X[2][k0] * E->data.layer_km;
    rho_prem_MO = prem_density_kg_m3(z_MO_km);
    if (rho_prem_MO < 1.0f) {
        fprintf(stderr, "[BMO_rho_nv] ERROR: invalid PREM anchor density at z_MO=%g km\n",
                (double)z_MO_km);
        return;
    }

    E->control.BMO_rho_nv_prof = (float *)calloc(noz + 1, sizeof(float));
    if (E->control.BMO_rho_nv_prof == NULL) {
        fprintf(stderr, "[BMO_rho_nv] ERROR: cannot allocate rho_nv profile (noz=%d)\n", noz);
        return;
    }

    rho_ref = E->control.BMO_rho_solid_ref;
    E->control.BMO_rho_nv_anchor = rho_ref;
    E->control.BMO_prem_rho_anchor = rho_prem_MO;

    for (k = 1; k <= noz; k++) {
        depth_km = E->X[2][k] * E->data.layer_km;
        rho_prem_k = prem_density_kg_m3(depth_km);
        E->control.BMO_rho_nv_prof[k] = rho_ref * rho_prem_k / rho_prem_MO;
    }
    /* Absolute anchor: rho_nv at MO row is the input ref, not PREM(z_MO). */
    E->control.BMO_rho_nv_prof[k0] = rho_ref;

    if (E->parallel.me != 0)
        return;

    coeff = 0.28f + 0.04f * E->control.BMO_n_spin;
    rho_nv_MO = E->control.BMO_rho_nv_prof[k0];
    {
        float KFe = E->control.BMO_partition_D;
        float xi_l = E->control.BMO_xi_l;
        float denom;
        if (xi_l >= 1.0f - 1.0e-6f)
            XFe_upper = 1.0f;
        else {
            denom = 1.0f - (1.0f - KFe) * xi_l;
            if (fabsf(denom) < 1.0e-12f)
                XFe_upper = 1.0f;
            else
                XFe_upper = KFe * xi_l / denom;
        }
    }
    XH2O_upper = E->control.BMO_Dw * E->control.BMO_Cw_l;
    XFe_melt = E->control.BMO_xi_l;
    XH2O_melt = E->control.BMO_Cw_l;

    rho_upper_MO = BMO_du2019_rho_kg_m3(rho_nv_MO, XFe_upper, XH2O_upper,
                                        E->control.BMO_n_spin, E->control.BMO_alpha_water);
    rho_melt_MO = BMO_du2019_rho_kg_m3(rho_nv_MO, XFe_melt, XH2O_melt,
                                       E->control.BMO_n_spin, E->control.BMO_alpha_water);

    fprintf(stderr,
        "[BMO_rho_nv] PREM scaling ON: z_MO_nom=%.4f (%.2f km)"
        "  z_MO_grid row=%d (%.2f km)"
        "  rho_nv_anchor=%.1f  PREM(z_MO_grid)=%.1f kg/m^3"
        "  (PREM shape only; anchor abs value is rho_ref)\n",
        (double)E->control.MO_depth_initial, (double)z_MO_nom_km,
        k0, (double)z_MO_km,
        (double)rho_ref, (double)rho_prem_MO);
    fprintf(stderr,
        "[BMO_rho_nv] rho_nv: surface=%.1f  MO_row=%.1f  CMB=%.1f kg/m^3\n",
        (double)E->control.BMO_rho_nv_prof[1],
        (double)E->control.BMO_rho_nv_prof[k0],
        (double)E->control.BMO_rho_nv_prof[noz]);
    fprintf(stderr,
        "[BMO_rho_nv] interface (rho_nv continuous):"
        "  rho_upper(MO-)=%.1f  rho_melt(MO+)=%.1f  dRho=%.1f kg/m^3"
        "  (XFe_sol=%.4g  xi_l=%.4g)\n",
        (double)rho_upper_MO, (double)rho_melt_MO,
        (double)(rho_melt_MO - rho_upper_MO),
        (double)XFe_upper, (double)XFe_melt);
    (void)coeff;
}

float lookup_BMO_rho_nv(E, z_ctr)
    struct All_variables *E;
    float z_ctr;
{
    int k, noz, lev;
    float *prof, *zcoord;
    float z1, z2, r1, r2, t;

    if (E->control.BMO_density_mode != 1
        || E->control.BMO_rho_nv_prem == 0
        || E->control.BMO_rho_nv_prof == NULL)
        return E->control.BMO_rho_solid_ref;

    lev = E->mesh.levmax;
    noz = E->lmesh.NOZ[lev];
    prof = E->control.BMO_rho_nv_prof;
    zcoord = E->X[2];

    if (z_ctr <= zcoord[1])
        return prof[1];
    if (z_ctr >= zcoord[noz])
        return prof[noz];

    for (k = 1; k < noz; k++) {
        z1 = zcoord[k];
        z2 = zcoord[k + 1];
        if (z_ctr >= z1 && z_ctr <= z2) {
            r1 = prof[k];
            r2 = prof[k + 1];
            if (fabs(z2 - z1) < 1.0e-12f)
                return r1;
            t = (z_ctr - z1) / (z2 - z1);
            return r1 * (1.0f - t) + r2 * t;
        }
    }
    return prof[noz];
}

float BMO_rho_grid_row(E, row_k, out_rho_nv, out_XFe, out_XH2O)
    struct All_variables *E;
    int row_k;
    float *out_rho_nv;
    float *out_XFe;
    float *out_XH2O;
{
    int lev, noz, elz, ielem, ifunc, ifunc_Fe, ifunc_H2O;
    float z_ctr, XFe, XH2O, rho_nv, rho;

    ifunc_Fe  = -1;
    ifunc_H2O = -1;
    lev = E->mesh.levmax;
    noz = E->lmesh.NOZ[lev];
    elz = noz - 1;
    if (row_k < 1)     row_k = 1;
    if (row_k > elz)   row_k = elz;

    ielem = row_k;
    z_ctr = 0.5f * (E->XL[2][row_k] + E->XL[2][row_k + 1]);

    for (ifunc = 0; ifunc < E->tracers.ntracfunc; ifunc++) {
        if (E->tracers.tracfunc[ifunc] == 80 && E->tracers.felem[ifunc] != NULL)
            ifunc_Fe = ifunc;
        if (E->tracers.tracfunc[ifunc] == 81 && E->tracers.felem[ifunc] != NULL)
            ifunc_H2O = ifunc;
    }

    XFe  = (ifunc_Fe  >= 0) ? E->tracers.felem[ifunc_Fe ][ielem]
                            : E->control.BMO_XFe_ref;
    XH2O = (ifunc_H2O >= 0) ? E->tracers.felem[ifunc_H2O][ielem] : 0.0f;

    rho_nv = lookup_BMO_rho_nv(E, z_ctr);
    rho = BMO_du2019_rho_kg_m3(rho_nv, XFe, XH2O,
                               E->control.BMO_n_spin, E->control.BMO_alpha_water);

    if (out_rho_nv != NULL) *out_rho_nv = rho_nv;
    if (out_XFe    != NULL) *out_XFe    = XFe;
    if (out_XH2O   != NULL) *out_XH2O   = XH2O;
    return rho;
}
