#ifndef HEX_H
#define HEX_H

struct point3d_str {
    double x[3];
};
typedef point3d_str POINT3D;

struct hexahedron {
    double centre[3];
    POINT3D vertex[8];
    double width[3];
    double vx[3];
    double vy[3];
    double vz[3];
};
typedef hexahedron HEXAHEDRON;

#endif // HEX_H
