/**
* DISPHOTN Dispersive Metamaterial Photonic Bands Calculator
* Version 0.1, March 2011
*
* This program is released under the GNU General Public License Version 3.
*
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>


// DEFINES
#define MAJOR 0 // Version Number
#define MINOR 1
#define NEARZERO_TOLERANCE (1e-4) // Tolerance used for geometrical comparisons
#define NEARZERO_STRICT_TOLERANCE (1e-14) // Tolerance used to correct rounding errors
#define MAX_SEQ_LEN 300 // Maximum sequence length in configuration file
#define PI (3.141592653589793115997963468544185161590576171875)
#define TWOPI (6.28318530717958623199592693708837032318115234375)
#define PENALTY (0.7)  // Default value for interface penalty
#define NUM_PROCESSORS 1 // Default number of processors, used by PARDISO.
#define ARPACK_MAX_ITERATIONS 100 // Maximum number of iterations for ARPACK.
#define TRUE (-1) // These are used in conjunction with the custom "bool" type.
#define FALSE (0)


// INCLUDES
#include "structs.c"
#include "utility.c"
#include "file_io.c"
#include "geometry_manipulation.c"
#include "matrix_manipulation.c"
#include "matrix_building.c"
#include "eigensolver.c"

int main( int argc, char *argv[] )
{

    // VARIABLE DECLARATIONS

    int i, j, k, l, *temp;
    char *dummy;
    point2d pt1, pt2, pt3;

    // File names
    char filename[256], meshfile[256], bandsfile[256], modesfile[256], hfieldfile[256], gridfile[256];

    // Configuration
    strlist *variables, *values;
    strlist *pvar, *pval;

    // Mesh information
    int npoints, nedges, ntries, nbedges;
    edge *edges, *bedges;
    triangle *tries;

    // Geometric Information
    point2d *a1, *a2;
    point2d *points, point_zero;
    point_zero.x = 0;
    point_zero.y = 0;
    bool xperiodic = FALSE, yperiodic = FALSE;
    byte *boundary_orientations_rl, *boundary_orientations_tb;
    bool *on_PEC;

    // Boundary Information
    int nbound[4];
    int *boundary[4]; // Left, Bottom, Right, Top
    int *boundary_info;
    point2d *p1, *p2, *par1, *par2;

    // Interface Information
    int ndomains, *interface_info;
    double *omegap_info;

    // Matrix Information
    doublecomplex *a,*b;
    int n, nnza, nnzb, *ia, *ja, *ib, *jb;
    int *keep_indices, *keep_indicesB, nKeep, nKeepB;
    doublecomplex *m, *f;
    int tot_size_a, tot_size_b;
    int ka, kb;
    int *iak, *jak, *ibk, *jbk, nnzak, nnzbk;
    doublecomplex *ak, *bk;

    // Interpolation grid
    int xgrid = 21, ygrid = 21, gridn;
    point2d *grid;
    doublecomplex *field, *hfield;
    field2d *modes;

    // Metal and dielectric information
    int *metal_domains, num_metal_domains, nomega0, nomegap;
    int *dielectric_domains, num_dielectric_domains, nepsilon;
    double *omegap, *omega0, *epsilon, eps;
    nomega0 = nomegap = num_metal_domains = num_dielectric_domains = nepsilon = 0;

    // K-points
    int kvals;
    double *kvx, *kvy;
    doublecomplex phase;

    // Solver Information
    int nev = 0, ncv = 0, solver = 1, processors = NUM_PROCESSORS;
    int iparam[11];
    doublecomplex *eigs, *eigv, *bands;
    double penalty = PENALTY;
    doublecomplex sigma;
    sigma.re = 0;
    sigma.i = 0;

    // Memory Allocation
    a1 = malloc(sizeof(point2d));
    a2 = malloc(sizeof(point2d));
    *a1 = *a2 = point_zero;
    kvx = malloc(MAX_SEQ_LEN*sizeof(double));
    kvy = malloc(MAX_SEQ_LEN*sizeof(double));

    // Sets all the filenames as empty strings
    bandsfile[0] = '\0';
    modesfile[0] = '\0';
    meshfile[0] = '\0';
    gridfile[0] = '\0';

    // Prints title
    printf("DISPHOTN Dispersive Photonic Bands v.%d.%d\n", MAJOR, MINOR);


    if (argc > 1)
    {
        // Reads configuration filename from the command line
        strcpy(filename, argv[1]);

    }
    else
    {
        // Inputs configuration filename
        printf("Insert configuration filename:\n");
        scanf("%s", filename);
    }

    // **************** READING CONFIGURATION **************
    printf("Reading configuration file %s\n", filename);


    // The configuration file is read and the information is stored in two lists.
    if(!read_configuration(filename, &variables, &values))
    {
        printf("ERROR: Error reading configuration file.\n");
        return 0;
    }

    // The list-parsing pointers pvar and pval are positioned at the beginning of the lists.
    pvar = variables;
    pval = values;

    // At each cycle, a parameter is analyzed and stored in memory.
    while (pvar != NULL)
    {
        // Identifies the parameter.
        if (pvar->v[0] == '#')
        {
            // The line is a comment: This should not happen, since comments should be already filtered
            // by the read_configuration_file function.
            printf("ERROR: Error in processing configuration file.\n");
            return 0;
        }
        else if (strcmp(pvar->v,"meshfile") == 0)
        {
            // Mesh file name. Required.
            strcpy(meshfile,pval->v);
        }
        else if (strcmp(pvar->v,"nev") == 0)
        {
            // Number of eigenvalues to compute. Required.
            nev = atoi(pval->v);
        }
        else if (strcmp(pvar->v,"ncv") == 0)
        {
            // Number of Arnoldi vectors, used by ARPACK.
            ncv = atoi(pval->v);
        }
        else if (strcmp(pvar->v,"processors") == 0)
        {
            // Number of processors to use, used by PARDISO.
            processors = atoi(pval->v);
            processors = processors < 1 ? 1 : processors;
        }
        else if ( strcmp(pvar->v,"sigmar") == 0)
        {
            // Real part of the eigenvalue shift.
            sigma.re = strtod(pval->v,&dummy);
            sigma.re = square(sigma.re * TWOPI);
        }
        else if (strcmp(pvar->v,"sigmai") == 0)
        {
            // Imaginary part of the eigenvalue shift.
            sigma.i = strtod(pval->v,&dummy);
            sigma.i = square(sigma.i * TWOPI);
        }
        else if (strcmp(pvar->v, "metal-domains") == 0)
        {
            // Sequence of the metallic domains. The same domain can appear
            // multiple times in the sequence; in this case, a multi-pole Lorentz
            // model is used.
            metal_domains = malloc(MAX_SEQ_LEN*sizeof(int));
            read_sequence(pval->v, metal_domains, &num_metal_domains, "%d");
        }
        else if (strcmp(pvar->v, "omegap") == 0)
        {
            // Sequence of the plasma frequencies for the metallic domains.
            // The values must correspond (in the same order) to the domains in the
            // metal-domains sequence.
            omegap = malloc(MAX_SEQ_LEN*sizeof(double));
            read_sequence(pval->v, omegap, &nomegap, "%lg");
        }
        else if (strcmp(pvar->v, "omega0") == 0)
        {
            // Sequence of the resonance frequencies for the metallic domains.
            // The values must correspond (in the same order) to the domains in the
            // metal-domains sequence.
            omega0 = malloc(MAX_SEQ_LEN*sizeof(double));
            read_sequence(pval->v, omega0, &nomega0, "%lg");
        }
        else if (strcmp(pvar->v, "dielectric-domains") == 0)
        {
            // Sequence of the dielectric domains.
            dielectric_domains = malloc(MAX_SEQ_LEN*sizeof(int));
            read_sequence(pval->v, dielectric_domains, &num_dielectric_domains, "%d");
        }
        else if (strcmp(pvar->v, "epsilon") == 0)
        {
            // Sequence of the relative permittivities for the domains
            // specified in the dielectric-domains sequence.
            epsilon = malloc(MAX_SEQ_LEN*sizeof(double));
            read_sequence(pval->v, epsilon, &nepsilon, "%lg");
        }
        else if (strcmp(pvar->v,"solver") == 0)
        {
            // PARDISO solver type. 0 = iterative, 1 = direct.
            solver = atoi(pval->v);
        }
        else if (strcmp(pvar->v,"bandsfile") == 0)
        {
            // Where to save the eigenfrequencies.
            strcpy(bandsfile,pval->v);
        }
        else if (strcmp(pvar->v,"modesfile") == 0)
        {
            // Where to save the interpolated electric field
            strcpy(modesfile,pval->v);
        }
        else if (strcmp(pvar->v,"hfieldfile") == 0)
        {
            // Where to save the magnetic field.
            strcpy(hfieldfile,pval->v);
        }
        else if (strcmp(pvar->v,"xperiodic") == 0)
        {
            // Periodicity along the a1 direction.
            xperiodic = -1;
        }
        else if (strcmp(pvar->v,"yperiodic") == 0)
        {
            // Periodicity along the a2 direction.
            yperiodic = -1;
        }
        else if (strcmp(pvar->v,"a1") == 0)
        {
            // First Bravais lattice vector. Required.
            sscanf(pval->v, "%lg %lg", &a1->x, &a1->y);
        }
        else if (strcmp(pvar->v,"a2") == 0)
        {
            // Second Bravais lattice vector. Required.
            sscanf(pval->v, "%lg %lg", &a2->x, &a2->y);
        }
        else if (strcmp(pvar->v,"kx") == 0)
        {
            // Sequence of the x components of the k-vector, normalized.
            read_sequence(pval->v, kvx, &kvals, "%lf");
        }
        else if (strcmp(pvar->v,"ky") == 0)
        {
            // Sequence of the y components of the k-vector, normalized.
            read_sequence(pval->v, kvy, &kvals, "%lf");
        }
        else if (strcmp(pvar->v,"xgrid") == 0)
        {
            // Size of the interpolation grid along the a1 direction.
            xgrid = atoi(pval->v);
        }
        else if (strcmp(pvar->v,"ygrid") == 0)
        {
            // Size of the interpolation grid along the a2 direction.
            ygrid = atoi(pval->v);
        }
        else if (strcmp(pvar->v,"penalty") == 0)
        {
            // Penalty coefficient.
            sscanf(pval->v, "%lg", &penalty);
        }
        else if (strcmp(pvar->v,"gridfile") == 0)
        {
            // Optional file containing the coordinate on which the electric field
            // must be interpolated. If set, overrides the xgrid and ygrid parameters.
            strcpy(gridfile,pval->v);
        }

        // Moves to the next list element
        pvar = pvar->n;
        pval = pval->n;

    }

    // Checks parameter values
    if (meshfile[0] == '\0')
    {
        // The mesh file name has not been set.
        printf("ERROR: Mesh file name missing from configuration file. Please set the \"meshfile\" parameter.");
        return 0;

    }

    if (nev <= 0)
    {
        // The number of required eigenvalues is not set or invalid.
        printf("ERROR: Invalid number of required eigenvalues. Please set the \"nev\" parameter to an integer > 0.\n");
        return 0;
    }

    // If the three sequences containing the information on metallic domains have different
    // lengths, they are shortened to the length of the shortest one, and a warning is issued.
    if( num_metal_domains != nomegap || num_metal_domains != nomega0)
    {
        printf("WARNING: The sequences metal-domains, omegap and omega0 have different length.\n");
        num_metal_domains = nomegap = nomega0 = (nomega0 < nomegap ? (nomega0 < num_metal_domains ? nomega0 : num_metal_domains) : (nomegap < num_metal_domains ? nomegap : num_metal_domains));
    }

    // The same thing for dielectric-domains and epsilon
    if( num_dielectric_domains != nepsilon)
    {
        printf("WARNING: The sequences dielectric-domains and epsilon have different length.\n");
        num_dielectric_domains = nepsilon = (num_dielectric_domains < nepsilon ? num_dielectric_domains : nepsilon);
    }

    printf("Metal Domains: %d, Dielectric Domains: %d\n", num_metal_domains, num_dielectric_domains);

    // Checks that no domain is set both as metallic and dielectric.
    for(i = 0; i < num_metal_domains; i++)
        for(j = 0; j < num_dielectric_domains; j++)
        {
            if(dielectric_domains[j] == metal_domains[i])
            {
                printf("ERROR: Domain %d is declared both metal and dielectric\n", metal_domains[i]);
                return 0;

            }
        }


    // (de)normalizes metallic material parameters.
    for(i = 0; i < num_metal_domains; i++)
    {
        omegap[i] *= TWOPI;
        omega0[i] *= TWOPI;
    }

    printf("Done\n");

    // (de)normalizes the value of the k-vector components
    printf("K-vector values: ");
    for (i = 0; i < kvals; i++)
    {
        printf("(%lg, %lg), ", kvx[i], kvy[i]);
        kvx[i] *= PI;
        kvy[i] *= PI;

    }
    printf("\n");

    // ****************** INTERPOLATION GRID *************
    if(gridfile[0] != '\0')
    {
        // If an interpolation grid file is provided, the points are
        // loaded from it.
        if(!read_grid_file_2d(gridfile, &grid, &gridn))
        {
            printf("ERROR: Error reading grid file.\n");
            return 0;

        }

    }
    else
    {
        // If no file is provided, a regular grid of the same shape of the unit
        // cell is used.
        gridn = xgrid*ygrid;
        grid = malloc(gridn*sizeof(point2d));
        k = 0;
        for(i = 0; i < xgrid; i++)
        {
            for(j = 0; j < ygrid; j++)
            {
                grid[k].x = a1->x*i/(xgrid - 1) + a2->x*j/(ygrid - 1);
                grid[k].y = a1->y*i/(xgrid - 1) + a2->y*j/(ygrid - 1);
                k++;
            }
        }
    }

    // Reading Mesh Files: the information on the mesh is read from the file
    // and saved in the dynamically allocated arrays.
    printf("Reading Mesh File %s\n", meshfile);
    i = read_mesh_file_2d(meshfile, &points, &npoints, &bedges, &nbedges, &tries, &ntries);
    if (i > 0)
    {
        printf("ERROR: Error %d in read_mesh_file_2d.\n", i);
        return 0;

    }
    printf("File read ended\n");

    printf("Points: %d, Boundary edges: %d, Triangles: %d\n", npoints, nbedges, ntries);

    // ************ EDGE IDENTIFICATION *********
    printf("Identifying edges\n");
    // The mesh file only provides information on nodes and triangles. The following function
    // assigns an unique index to each edge and saves, for each triangle, the indices of the three
    // edges.
    i = identify_edges(&edges, &nedges, tries, &ntries);

    if (i > 0)
    {
        // Some error has occurred. The user should look up in the documentation the error code.
        printf("ERROR: Error %d in edge identification.\n", i);
        return 0;
    }

    printf("Done. Created %d edges\n", nedges);

    // ************ EDGE LENGTHS ************

    // Computes the edge lengths.
    for(i = 0; i < nedges; i++)
    {
        pt1 = points[edges[i].p[0]];
        pt2 = points[edges[i].p[1]];

        edges[i].length = sqrt(square(pt1.x - pt2.x) + square(pt1.y - pt2.y));

    }

    // ************ TRIANGLE AREAS ************

    // Computes the triangle areas using the determinant method.
    for(i = 0; i < ntries; i++)
    {
        pt1 = points[tries[i].p[0]];
        pt2 = points[tries[i].p[1]];
        pt3 = points[tries[i].p[2]];

        tries[i].area = 0.5*(pt2.x*pt3.y + pt3.x*pt1.y + pt1.x*pt2.y - pt2.x*pt1.y - pt3.x*pt2.y - pt1.x*pt3.y);


    }


    // ************ BOUNDARY IDENTIFICATION *************
    // Identifies the edges belonging to the boundaries of the unit cells. To do so,
    // the a1 and a2 vectors are used.
    printf("Identifying boundaries\n");

    // Since in a 2D problem each edge can lie at most on one boundary, we save the index of the
    // boundary for each edge (or -1 for no boundary) in the boundary_info array.
    boundary_info = malloc(nedges*(num_metal_domains + 1)*sizeof(int));

    // The on_PEC array contains a true/false (-1 or 0) value for each edge. The
    // boundary array (created later) contains the indices of the edges on each boundary.

    on_PEC = malloc(nedges*sizeof(bool));

    // Initializes the array
    for (i = 0; i < nedges; i++)
        on_PEC[i] = FALSE;

    // Processes the four boundaries
    for (i = 0; i < 4; i++)
    {
        // k will be set to true/false if the current boundary is a PEC.
        k = FALSE;

        // The process of identifying the edges on a boundary is the same for each
        // boundary. The two points constituting the edge are translated by a vector
        // -par2 and compared to the vector par1. If the three vectors (the two translated
        // points and par1) are collinear, then the edge is on the boundary. The difference
        // between the 4 boundaries lie only in the values of par1 and par2.
        switch (i)
        {
            case 0: // Left boundary
            {
                par1 = a2; // The edge must be collinear to the a2 vector...
                par2 = &point_zero; //... without being translated.

                // The Left boundary is PEC if there is no x-periodicity.
                k = !xperiodic;
            }
            break;
            case 1: // Bottom boundary
            {
                par1 = a1; // The same as above, but now with the a1 vector
                par2 = &point_zero;

                k = !yperiodic;
            }
            break;
            case 2: // Right boundary
            {
                par1 = a2; // The edge must be collinear with the a2 vector...
                par2 = a1; // ...after a translation of -a1.

                k = !xperiodic;
            }
            break;
            case 3: // Top boundary
            {
                par1 = a1; // Same as above.
                par2 = a2;

                k = !yperiodic;
            }
            break;
        }

        // nbound[i] will contain the number of edges on the i-th boundary,
        nbound[i] = 0;
        // while boundary[i] is an array containing the indices of the edges on the i-th boundary
        boundary[i] = malloc(nedges*sizeof(int));
        for (j = 0; j < nedges; j++)
        {

            // Retrieves the coordinates of the edge points.
            p1 = &points[edges[j].p[0]];
            p2 = &points[edges[j].p[1]];

            // Checks the collinearity of the translated points with par1, by computing the determinant
            if (near_zero(det2(vector_subtract(*p1, *par2), *par1)) && near_zero(det2(vector_subtract(*p2, *par2), *par1)))
            {

                // Adds the edge to the boundary
                boundary[i][nbound[i]] = j;
                (nbound[i])++;

                on_PEC[j] = k;
            }
        }

        if(nbound[i] == 0)
        {
            // Checks if there is at least one edge on each boundary
            // If this fails, usually the user specified wrong Bravais lattice vectors

            printf("ERROR: Boundary %d has zero edges on it. Please check the values of a1 and a2 vectors.\n", i);
            return 0;

        }

        // Resizes the boundary array to the actual size
        temp = realloc(boundary[i], nbound[i]*sizeof(int));
        if (temp == NULL)
        {
            printf("ERROR: Realloc failed\n");
            free(boundary[i]);
            return 0;
        }
        boundary[i] = temp;
        temp = NULL;
    }

    // If periodicity is to be imposed, there must be the same number of edges
    // on the two corresponding boundaries.
    if (xperiodic && (nbound[0] != nbound[2]))
    {
        printf("ERROR: Wrong edge count on left-right boundaries.\n");
        return 0;

    }
    if (yperiodic && (nbound[1] != nbound[3]))
    {
        printf("ERROR: Wrong edge count on top-bottom boundaries.\n");
        return 0;
    }
    printf("Done. Edges on boundaries: Left %d, Bottom %d, Right %d, Top %d\n", nbound[0], nbound[1], nbound[2], nbound[3]);

    // ************ PERIODIC CONDITIONS *************
    // If the periodic conditions must be imposed, the program identifies the
    // pairs of corresponding edges on the corresponding boundaries.


    printf("Periodicity: ");
    if(!xperiodic && !yperiodic)
        printf("None\n");
    else
    {
        if(xperiodic)
            printf("x ");
        if(yperiodic)
            printf("y");
        printf("\n");
    }


    if (xperiodic)
    {

        printf("Identifying Left-Right Periodic Conditions\n");
        // The following array identifies the relative orientation between corresponding edges.
        // It will contain only 1 or -1.
        boundary_orientations_rl = (byte*)malloc(nbound[0]*sizeof(byte));

        i = identify_periodic_boundaries(points, edges, boundary, nbound, 2, 0, a1, boundary_orientations_rl);
        if (i >= 0)
        {
            // An error has occurred, the error code is the number of the unpaired edge.
            printf("ERROR: No corrispondence found for edge %d on left boundary.\n", i);
            return 0;
        }

        printf("Done.\n");
    }

    if (yperiodic)
    {
        printf("Identifying Bottom-Top Periodic Conditions\n");

        boundary_orientations_tb = (byte*)malloc(nbound[1]*sizeof(byte));

        i = identify_periodic_boundaries(points, edges, boundary, nbound, 3, 1, a2, boundary_orientations_tb);
        if (i >= 0)
        {
            printf("ERROR: No corrispondence found for edge %d on bottom boundary.\n", i);
            return 0;
        }

        printf("Done.\n");
    }

    // ************* INTERFACES ***************
    // To use the penalty method to impose the interface condition, the edges
    // belonging to the interfaces must be identified.
    printf("Identifying Interfaces\n");

    // Counts the total number of domains (by finding the maximum domain index).
    ndomains = 0;
    for(i = 0; i < ntries; i++)
    {
        if(tries[i].domain > ndomains)
            ndomains = tries[i].domain;
    }
    ndomains++; // Domain 0 is outside the mesh.

    // The omegap_info will contain the values of omegap for each domain.
    // The freespace or dielectric domains will have omegap = 0.
    omegap_info = malloc((ndomains)*sizeof(double));


    for(i = 0; i < ndomains; i++)
    {
        for(j = 0; j < num_metal_domains; j++)
        {
            if(metal_domains[j] == i)
                omegap_info[i] = omegap[j];
            else
                omegap_info[i] = 0;
        }
    }

    // The interface_info array will contain, for each edge, the indices
    // of the two domains bounding it. The interface edges will be the ones with two
    // different domain indices. The interface_info array has twice as elements as the
    // number of edges.
    interface_info = malloc(nedges*2*sizeof(int));

    for(i = 0; i < ntries; i++)
    {
        for(j = 0; j < 3; j++)
        {
            // Each edge is registered with orientation 1 for one of the two bounding triangles,
            // and -1 for the other. The domain index of the "1"-triangle will be stored in the
            // first position in interface_info, the domains index of the "-1"-triangle will be
            // stored in the second position.
            if(tries[i].o[j] == 1)
            {
                interface_info[tries[i].e[j]] = tries[i].domain;
            }
            else
            {

                interface_info[tries[i].e[j] + nedges] = tries[i].domain;
            }
        }
    }

    // Counts the number of edges on the interfaces.
    int ninterface = 0;
    for(i = 0; i < nedges; i++)
    {
        if(interface_info[i] != interface_info[i + nedges])
            ninterface++;
    }


    printf("Done. Edges on Interfaces: %d\n", ninterface);


    // ************* MATRIX CONSTRUCTION **************
    printf("Constructing Matrices\n");

    // During the process of matrix construction, a number of matrix entries
    // (consisting of the row and column indices i,j and the complex value) is
    // created. Many of the entries will correspond to the same position i,j;
    // this will be fixed later.

    // For each triangle, a 3x3 elemental matrix is computed. This accounts for
    // 9 matrix entries. Each Lorentz model requires an expansion of the matrix,
    // so 9*4 more entries are required for each metal domain for each triangle.
    // Moreover, in the A matrix, the entries due to the interface penalty method
    // must be stored as well.
    tot_size_a = 9*ntries*(num_metal_domains*4+1) + ninterface;
    tot_size_b = 9*ntries*(num_metal_domains+1);

    // Allocation of the row indices array (ix), column indices array (jx) and values array (x)
    ia = (int*)malloc(tot_size_a*sizeof(int));
    ja = (int*)malloc(tot_size_a*sizeof(int));
    ib = (int*)malloc(tot_size_b*sizeof(int));
    jb = (int*)malloc(tot_size_b*sizeof(int));
    a = (doublecomplex*)malloc((tot_size_a)*sizeof(doublecomplex));
    b = (doublecomplex*)malloc((tot_size_b)*sizeof(doublecomplex));

    if(ia == NULL || ja == NULL || ib == NULL || jb == NULL || a == NULL || b == NULL)
    {
        // Something went wrong during the allocation.
        printf("ERROR: Allocation Error.\n");
        return 0;
    }

    // Allocates space for the elemental matrices.
    m = (doublecomplex*)malloc(9*sizeof(doublecomplex));
    f = (doublecomplex*)malloc(9*sizeof(doublecomplex));

    // ka and kb are counters that are incremented with each entry created. The values of
    // ka and kb at the end will be the actual number of matrix entries.
    ka = 0;
    kb = 0;

    for (i = 0; i < ntries; i++)
    {
        // The matrix construction is handled by the matrix_construction function,
        // that computes the elemental matrices.
        matrix_construction(points, edges, &tries[i], m, f);

        // If the triangle belongs to a dielectric domain, a non-unitary relative
        // permittivity must be applied. The value of the relative permittivity is
        // retrieved from the epsilon array.
        eps = 1;
        for(j = 0; j < num_dielectric_domains && eps == 1; j++)
        {
            if(tries[i].domain == dielectric_domains[j])
            {
                eps = epsilon[j];
            }
        }

        // The values of the elemental matrices are inserted into the arrays that
        // will form the sparse matrices. See the code of add_to_matrix for more information
        // on the parameters.
        add_to_matrix(ia, ja, a, m, &tries[i], 0, 0, 1, &ka);
        add_to_matrix(ib, jb, b, f, &tries[i], 0, 0, eps, &kb);

        // If the triagle is not into a dielectric domain, it might belong to a metallic
        // domain.
        if(eps == 1)
        {
            for(j = 0; j < num_metal_domains; j++)
            {
                if(tries[i].domain == metal_domains[j])
                {
                    // In this case, the matrix must be expanded with the
                    // P-field part.

                    // Omegap^2
                    add_to_matrix(ia, ja, a, f, &tries[i], 0, 0, omegap[j]*omegap[j], &ka);

                    // E-P
                    add_to_matrix(ia, ja, a, f, &tries[i], 0, (j+1)*nedges, -omegap[j]*omega0[j], &ka);

                    // P-E (f is hermitian)
                    add_to_matrix(ia, ja, a, f, &tries[i], (j+1)*nedges, 0, -omegap[j]*omega0[j], &ka);

                    // P-P
                    add_to_matrix(ia, ja, a, f, &tries[i], (j+1)*nedges, (j+1)*nedges, omega0[j]*omega0[j], &ka);

                    // B matrix
                    add_to_matrix(ib, jb, b, f, &tries[i], (j+1)*nedges, (j+1)*nedges, 1, &kb);
                }
            }
        }

        // Interface conditions.
        // Checks if the edges are on an interface. If this is the case, the penalty
        // method is enforced.
        for(k = 0; k < 3; k++)
        {

            if(interface_info[tries[i].e[k]] != interface_info[tries[i].e[k] + nedges])
            {
                // The edge is on the interface. Apply the correct penalty according
                // to the interface properties.

                if(tries[i].o[k] == 1)
                {

                    double coef = omegap_info[interface_info[tries[i].e[k] + nedges]];
                    coef *= coef*penalty;

                    if(coef != 0)
                    {


                        ja[ka] = ia[ka] = tries[i].e[k];
                        a[ka].re = coef*f[4*k].re;
                        a[ka].i = coef*f[4*k].i;
                        ka++;
                    }
                }
                else
                {

                    double coef = omegap_info[interface_info[tries[i].e[k]]];
                    coef *= penalty;

                    if(coef != 0)
                    {

                        ja[ka] = ia[ka] = tries[i].e[k];
                        a[ka].re = coef*f[4*k].re;
                        a[ka].i = coef*f[4*k].i;
                        ka++;
                    }
                }

            }

        }


    }
    // The sizes of the matrices' arrays is updated to the actual values.
    // The arrays will be resized later (in the assembly_matrix function),
    // so it is not necessary to do it here.
    tot_size_a = ka;
    tot_size_b = kb;

    printf("Assembling matrices\n");

    // The matrix are assembled by summing all the values of the entries with the same i,j position.
    // The matrix entries are first sorted by row and then column, using the sort_matrix_entries
    // function. The assembly_matrix function, then, sums all the consecutive duplicate entries.

    // Assembly matrix A

    // Please see the code for sort_matrix_entries for an explanation of the parameters.
    sort_matrix_entries(ia, ja, a, 0, tot_size_a - 1);

    sort_matrix_entries(ib, jb, b, 0, tot_size_b - 1);

    printf("Assembling matrix A\n");

    // The assembly_matrix function sums all the values sharing the same i,j position and
    // return the number of nonzeros elements in the matrix (i.e. the number of matrix entries.)
    nnza = assembly_matrix(&ia, &ja, &a, tot_size_a);

    if (nnza < 0)
    {
        // If the function returned a negative value, an error occurred.
        printf("ERROR: Error %d in assembly_matrix on matrix A.\n", nnza);
        return 0;
    }

    // Assembly matrix B
    printf("Assembling matrix B\n");

    // Same as above.
    nnzb = assembly_matrix(&ib, &jb, &b, tot_size_b);

    if (nnzb < 0)
    {
        printf("ERROR: Error %d in assembly_matrix on matrix B.\n", nnzb);
        return 0;
    }

    printf("Done.\n");

    // ************* PEC BOUNDARY CONDITIONS **************
    // Two kind of boundary conditions are currently supported:
    // PEC (Perfect electric conduction) and Bloch-periodicity.
    // In the first case, the tangential E-field is enforced to be
    // null; this correspond to the removal of the corresponding DOF.
    // In the second case, the DOF on a boundary are removed by
    // setting them proportional (with a phase factor) to the DOF
    // on the corresponding boundary.

    // The PEC Boundary conditions, being indipendent of the value of the
    // k-vector, are imposed here. The Periodic Boundary Conditions
    // will be imposed later for each value of the k-vector.

    if (!xperiodic || !yperiodic)
    {

        printf("Imposing PEC Boundary Conditions\n");
        // The following functions simply set to zero the values on the rows and columns
        // corresponding to the edges on a PEC. After this process, the matrices will then
        // have empty rows and columns. These will be removed (and the matrices will be shrunk
        // acordingly) after imposing the Periodic Bundary conditions.

        set_PEC_conditions(ia, ja, a, nnza, on_PEC, nedges);

        set_PEC_conditions(ib, jb, b, nnzb, on_PEC, nedges);

        printf("Done\n");

    }

    // If there is no periodicity, there is no point in cycling over different values of the
    // k vector. In this case, the parameters kvals, kx and ky are overridden.
    if(!xperiodic && !yperiodic)
        kvals = 1;

    printf("Number of processors to be used for linear system solving: %d\n", processors);
    printf("Beginning computation for %d k points\n", kvals);

    // Allocates space to store the computation results.
    bands = malloc(kvals*nev*sizeof(doublecomplex));
    modes = malloc(kvals*nev*gridn*sizeof(field2d));
    hfield = malloc(kvals*nev*ntries*sizeof(doublecomplex));


    // Begin cycling over the k-values.
    for (k = 0; k < kvals; k++)
    {

        printf("-------------------------------------------------------------\n k-cycle %d/%d, k value (%g , %g)\n", k + 1, kvals, kvx[k], kvy[k]);

        // Since the matrices must be modified by imposing the periodic boundary
        // conditions, a copy of the original matrices is stored in the ...k arrays
        // at each cycle.

        nnzak = nnza;
        nnzbk = nnzb;

        iak = malloc(nnzak*sizeof(int));
        jak = malloc(nnzak*sizeof(int));
        ak = malloc(nnzak*sizeof(doublecomplex));
        ibk = malloc(nnzbk*sizeof(int));
        jbk = malloc(nnzbk*sizeof(int));
        bk = malloc(nnzbk*sizeof(doublecomplex));

        // Copies the matrices.
        for(i=0; i < nnzak; i++)
        {
            iak[i] = ia[i];
            jak[i] = ja[i];
            ak[i] = a[i];
        }

        for(i=0; i < nnzbk; i++)
        {
            ibk[i] = ib[i];
            jbk[i] = jb[i];
            bk[i] = b[i];
        }

        // Imposition of Periodic Boundary conditions.
        if(xperiodic || yperiodic)
        {
            printf("Imposing periodic boundary conditions: ");

            if (xperiodic)
            {
                // Computes the Bloch phase shift for the periodic conditions.
                printf("x ");
                phase.re = cos(kvx[k]*a1->x + kvy[k]*a1->y);
                phase.i = sin(kvx[k]*a1->x + kvy[k]*a1->y);

                printf("(phase = %g+%gi) ", phase.re, phase.i);

                // The following functions express the values of the field on the edges beloging to a destination boundary as
                // the values of the field on the edges on the source boundary, multiplied by the phase factor and +1 or -1 to
                // take into account the relative orientation.
                // The rows and columns corresponding to the edges on the source boundary are set to zero, and will be removed
                // from the matrix later.
                set_periodic_conditions(iak, jak, ak, nnzak, nedges, nbound[0], boundary[0], boundary[2], boundary_orientations_rl, phase);

                set_periodic_conditions(ibk, jbk, bk, nnzbk, nedges, nbound[0], boundary[0], boundary[2], boundary_orientations_rl, phase);


            }

            if (yperiodic)
            {
                printf("y ");

                //Same as above.
                phase.re = cos(kvx[k]*a2->x + kvy[k]*a2->y);
                phase.i = sin(kvx[k]*a2->x + kvy[k]*a2->y);

                printf("(phase = %g+%gi) ", phase.re, phase.i);


                set_periodic_conditions(iak, jak, ak, nnzak, nedges, nbound[1], boundary[1], boundary[3], boundary_orientations_tb, phase);

                set_periodic_conditions(ibk, jbk, bk, nnzbk, nedges, nbound[1], boundary[1], boundary[3], boundary_orientations_tb, phase);




            }
            printf("\n");

            // Repeats the sorting-summing process done after the matrix construction.
            sort_matrix_entries(iak, jak, ak, 0, nnzak - 1);

            sort_matrix_entries(ibk, jbk, bk, 0, nnzbk - 1);

            nnzak = assembly_matrix(&iak, &jak, &ak, nnzak);

            if (nnzak < 0)
            {
                printf("ERROR: Error %d in assembly_matrix on matrix A.\n", nnzak);
                return 0;
            }

            nnzbk = assembly_matrix(&ibk, &jbk, &bk, nnzbk);

            if (nnzbk < 0)
            {
                printf("ERROR: Error %d in assembly_matrix on matrix B.\n", nnzbk);
                return 0;
            }
        }

        // Removes the empty rows and columns introduced by the boundary conditions (PEC
        // or Periodic) and shrinks the matrix accordingly.
        i = purge_matrix(&iak,&jak,&ak,&nnzak,&nKeep,&keep_indices);

        if (i)
        {
            printf("ERROR: Error %d in purge_matrix on matrix A.\n", i);
            return 0;
        }

        i = purge_matrix(&ibk,&jbk,&bk,&nnzbk,&nKeepB,&keep_indicesB);
        if (i)
        {
            printf("ERROR: Error %d in purge_matrix on matrix A.\n", i);
            return 0;
        }

        // Performs a basic check on the validity of the boundary conditions:
        // after the imposition, the size of the matrices A and B must be the same...
        if (nKeep != nKeepB)
        {
            printf("ERROR: Inconsistency in boundary conditions (edge number). A: %d, B: %d\n", nKeep, nKeepB);
            return 0;
        }
        else
        {
            // ... and, also, the edges not eliminated by the boundary conditions must be the same.
            for (i=0;i<nKeep;i++)
            {
                if (keep_indices[i] != keep_indicesB[i])
                {
                    printf("ERROR: Inconsistency in boundary conditions (edges).\n");
                    return 0;
                }
            }

        }

        // Everything went fine, the matrices are now ready.
        printf("Done. Matrix size: %d, nnza: %d, nnzb: %d\n", nKeep, nnzak, nnzbk);

        n = nKeep;

/*
        // MATRIX WRITING
        printf("Writing matrices \n");
        fid = fopen("/media/Volume/Dropbox/afile","w");

        for(i = 0; i < nnzak; i++)
        {
            if(near_zero_tol(ak[i].i))
                fprintf(fid, "%d %d %.50lg\r\n", iak[i], jak[i], ak[i].re);
            else
                fprintf(fid, "%d %d %.50lg%+.50lgi\r\n", iak[i], jak[i], ak[i].re, ak[i].i);

        }
        fclose(fid);
        fid = fopen("/media/Volume/Dropbox/bfile","w");
        for(i = 0; i < nnzbk; i++)
        {
            if(near_zero_tol(ak[i].i))
                fprintf(fid, "%d %d %.60lg\r\n", ibk[i], jbk[i], bk[i].re);
            else
                fprintf(fid, "%d %d %.60lg%+.60lgi\r\n", ibk[i], jbk[i], bk[i].re, bk[i].i);

        }
        fclose(fid);

        fid = fopen("/media/Volume/Dropbox/lengthsfile","w");
        for(i=0;i<nedges;i++)
            fprintf(fid, "%.60lg\n", edges[i].length);

        fclose(fid);


        fid = fopen("/media/Volume/Dropbox/areasfile","w");
        for(i=0;i<ntries;i++)
            fprintf(fid, "%.60lg\n", tries[i].area);

        fclose(fid);


        // END MATRIX WRITING
*/

        // PARDISO uses the Compressed Row Format.
        printf("Compressing Rows\n");
        compress_rows(iak, jak, ak, nnzak);
        compress_rows(ibk, jbk, bk, nnzbk);

        printf("Starting Eigensolver. Size of the problem: %d\n", n);

        // Allocates space for the current cycle's solutions
        eigs = (doublecomplex*)malloc(nev*sizeof(doublecomplex));
        field = (doublecomplex*)malloc(nedges*sizeof(doublecomplex));

        // Starts the eigensolver.
        i = eigensolver(&iak, &jak, &ak, nnzak, &ibk, &jbk, &bk, nnzbk, n, sigma, nev, ncv, eigs, &eigv, solver, iparam, processors);

        if (i)
        {
            // The eigensolver has returned an error.
            return 0;
        }


        // Frees the memory containing the modified matrices for this cycle.
        free(iak);
        free(jak);
        free(ak);
        free(ibk);
        free(jbk);
        free(bk);

        printf("Done.\n%d converged eigenvalues.\n", iparam[4]);

        if(iparam[4] < nev)
        {
            // ARPACK reported incomplete convergence.
            printf("WARNING: Not all the required eigenvalues have converged.\n");
        }

        // Normalizes eigenvalues and computes the eigenfrequencies.
        for(i=0;i<nev;i++)
        {
            // Eigenvalues should be real.
            eigs[i].re = sqrt((eigs[i].re >= 0) ? eigs[i].re : -eigs[i].re )/TWOPI;
            eigs[i].i = 0;
        }

        // Field interpolation
        printf("Interpolating Field\n");
        for(l = 0; l < nev; l++)
        {
            // The eigv vector contains the field values on the edges not removed during
            // the imposition of boundary conditions. We will now reconstruct the value
            // of the field on the missing edges by putting it to 0 (if they are on PEC)...
            j = 0; // j is the counter that scans the array keep_indices.
            for(i = 0; i < nedges; i++)
            {
                if(j < n && keep_indices[j] == i)
                {
                    field[i] = eigv[l*n + j];

                    j++;
                }
                else
                {
                    field[i].re = field[i].i = 0;
                }
            }

            // ... or by computing it from other edges (if they are on periodic boundaries).
            if(xperiodic)
            {
                phase.re = cos(kvx[k]*a1->x + kvy[k]*a1->y);
                phase.i = sin(kvx[k]*a1->x + kvy[k]*a1->y);

                for(j = 0; j < nbound[0]; j++)
                {

                    field[boundary[0][j]].re = boundary_orientations_rl[j]*(field[boundary[2][j]].re*phase.re - field[boundary[2][j]].i*phase.i);
                    field[boundary[0][j]].i = boundary_orientations_rl[j]*(field[boundary[2][j]].re*phase.i + field[boundary[2][j]].i*phase.re);


                }
            }

            if(yperiodic)
            {
                phase.re = cos(kvx[k]*a2->x + kvy[k]*a2->y);
                phase.i = sin(kvx[k]*a2->x + kvy[k]*a2->y);

                for(j = 0; j < nbound[1]; j++)
                {

                    field[boundary[1][j]].re = boundary_orientations_tb[j]*(field[boundary[3][j]].re*phase.re - field[boundary[3][j]].i*phase.i);
                    field[boundary[1][j]].i = boundary_orientations_tb[j]*(field[boundary[3][j]].re*phase.i + field[boundary[3][j]].i*phase.re);

                }
            }


            // The complete reconstructed field is now interpolated on the provided grid.
            interpolate_2d(grid, gridn, tries, ntries, points, field, &modes[(k*nev + l)*gridn]);

            // Computation of the z component of the magnetic field.
            // Retrieves the pointer to the part of the hfield array in which the
            // current cycle's H-field must be saved.
            doublecomplex *hz = &hfield[(k*nev + l)*ntries];


            // Computes the H field in each triangle. The H field is not interpolated,
            // as the current output can be used with MATLAB's pdeplot function.
            for(i = 0; i < ntries; i++)
            {
                hz[i].re = hz[i].i = 0;
                for(j = 0; j < 3; j++)
                {
                    double coef = edges[tries[i].e[j]].length*tries[i].o[j]/tries[i].area;
                    hz[i].re += coef*field[tries[i].e[j]].re;
                    hz[i].i += coef*field[tries[i].e[j]].i;
                }
            }


        }
        printf("Done.\n");


        // Copies the computed eigenfrequencies in the results array.
        for(i = 0; i < nev; i++)
        {
            bands[k*nev + i] = eigs[i];
        }

        free(eigs);
        free(field);

    }// End kvals cycle

    // ***************** WRITING RESULTS *************

    // In the following, the results (eigenfrequencies, electric field, magnetic field)
    // are written to file, if a valid filename is specified.

    if(bandsfile[0] != '\0')
    {
        printf("Writing Bands to File %s\n", bandsfile);
        i = write_bands(bandsfile, bands, nev, kvals);
        if(i)
        {
            printf("ERROR: Error in opening bands file.\n");
            return 0;
        }
        printf("Done.\n");
    }

    if(modesfile[0] != '\0')
    {
        printf("Writing E-Modes to File %s\n", modesfile);
        i = write_modes(modesfile, modes, grid, gridn, nev, kvals);
        if(i)
        {
            printf("ERROR: Error in opening modes file.\n");
            return 0;
        }
        printf("Done.\n");
    }

    if(hfieldfile[0] != '\0')
    {

        printf("Writing H-Modes to File %s\n", hfieldfile);
        i = write_hfield(hfieldfile, hfield, ntries, nev, kvals);
        if(i)
        {
            printf("ERROR: Error in opening modes file.\n");
            return 0;
        }
        printf("Done.\n");
    }


    // End of the program.
    printf("End.\n");

    // TO DO: Free all
    return -1;

}

