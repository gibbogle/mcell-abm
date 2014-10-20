#ifndef GLOBAL_H
#define GLOBAL_H

#include <QtGui>

#define MAX_CONC 4
#define MAX_CELLS 10000
#define N_CELLINFO 7

#include "hex.h"

namespace Global
{
//    extern int data1;
//    extern int data2;

//    extern int MAX_CHEMO;
    extern int Ncirc, Nlong;
    extern double DELTA_T;
//    extern double dfraction;
    extern int nt_vtk;
    extern int istep;
    extern bool leftb;

//    extern int nvars_used;
//    extern int GUI_to_DLL_index[32];
//    extern int DLL_to_GUI_index[32];
//    extern QString var_string[32];

//    extern double *FACS_data;
//    extern int nFACS_cells;
//    extern int nFACS_dim;

//    extern double *histo_data;
//    extern int nhisto_boxes;
//    extern int nhisto_dim;
//    extern double histo_vmax[3*32];
//    extern int histo_celltype;

    extern int summaryData[100];
//    extern int i_hypoxia_cutoff;
//    extern int i_growth_cutoff;

    extern double concData[4000];
    extern int conc_nc;
    extern double conc_dx;

//    extern double volProb[100];
//    extern int vol_nv;
//    extern double vol_v0;
//    extern double vol_dv;
//    extern double oxyProb[100];
//    extern int oxy_nv;
//    extern double oxy_dv;

    extern int cell_list[N_CELLINFO*MAX_CELLS];
    extern int ncell_list;

    extern HEXAHEDRON *hex_list;
    extern int Ncirc;
    extern int Nlong;
    extern int Nhex;

    extern bool showingVTK;
    extern bool recordingVTK;

} // namespace Global

#endif // GLOBAL_H
