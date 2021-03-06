#include "global.h"

// Note that in the Fortran DLL the constituent numbering starts at 1
// WARNING: the constituent indices are hard-wired - must be maintained consistent with
// cells-abm.dll

namespace Global
{
    int data1=0;
    int data2;

    int MAX_CHEMO;
    double DELTA_T;
    double dfraction;
    int nt_vtk;
    int istep;
    bool leftb;

    int nvars_used;
    int GUI_to_DLL_index[32];
    int DLL_to_GUI_index[32];
    QString var_string[32];

    double *FACS_data=NULL;
    int nFACS_cells=0;
    int nFACS_dim=0;

    double *histo_data=NULL;
    int nhisto_boxes=20;
    int nhisto_dim=0;
    double histo_vmax[3*32];
    int histo_celltype=0;

    int summaryData[100];
    int i_hypoxia_cutoff;
    int i_growth_cutoff;

    double concData[4000];
    int conc_nc;
    double conc_dx;

    double volProb[100];
    int vol_nv;
    double vol_v0;
    double vol_dv;
    double oxyProb[100];
    int oxy_nv;
    double oxy_dv;

    int cell_list[N_CELLINFO*MAX_CELLS];
    int ncell_list;

    HEXAHEDRON *hex_list;
    int Ncirc;
    int Nlong;
    int Nhex;

    bool showingVTK = false;
    bool recordingVTK = false;


} // namespace Global
