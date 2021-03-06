// Structs (Data types) used in the program

// Bool is the shortest possible data type (ideally should be 1-bit)
typedef char bool;
// Byte is the shortest value possible. Must be able to contain -1, 0, and 1.
typedef char byte;

// Double-precision complex number. This data structure is shared with the
// FORTRAN functions (PARDISO, ARPACK and BLAS).
typedef struct
{
    double re;
    double i;
}
doublecomplex;

// String List element. Used to store configuration parameters.
struct sl
{
    char *v;
    struct sl *n;
};
typedef struct sl strlist;

// Point in 2-dimensional space.
typedef struct
{
    double x;
    double y;
}
point2d;

// Field (electric, polarization etc.) in 2-dimensional space
typedef struct
{
    doublecomplex x;
    doublecomplex y;
}
field2d;

// For the geometrical structs, the first n elements must be the n nodes number (ints),
// and the (n+1)th must be the domain number (int). All these elements must be int. (see mesh_reading.c)

// Edge (2-nodes) element. Contains the 2 nodes number and the domain number.
typedef struct
{
    int p[2]; // DO NOT MODIFY
    //int p2; // DO NOT MODIFY
    int domain; // DO NOT MODIFY
    double length;
}
edge;

// Triangular (3-nodes) element. Contains the 3 nodes number, the domain number,
// the edges number (ints), and the orientation of the edges (+1 or -1).
typedef struct
{
    int p[3]; // DO NOT MODIFY
    int domain; // DO NOT MODIFY
    int e[3];
    byte o[3];
    double area;

}
triangle;
