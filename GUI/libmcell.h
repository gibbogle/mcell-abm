#ifndef LIBMCELL_H
#define LIBMCELL_H

#define DllExport __declspec( dllexport )

#ifdef __cplusplus
extern "C" {
#endif

#include "hex.h"

void execute(int *, char *, int *,char *, int *);
void simulate_step(int *);
void terminate_run(int *);
void get_dimensions(int *, int *, int *, double *);
void get_scene(int *, HEXAHEDRON *);
void get_summary(int *);
void get_concdata(int *, double *, double *);
void CellML_load_file(char *, int, int *, int *);
void CellML_get_component_info(int, char *, int *);
void CellML_get_variable_name(int, int, char *);
void CellML_get_variable_value(int, int, double *);
void CellML_set_variable_value(int, int, double);

#ifdef __cplusplus
}
#endif

#endif // LIBMCELL_H
