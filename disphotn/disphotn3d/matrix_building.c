void matrix_construction_3d(point3d *points, edge *edges, point3d *edge_versors, tetra *tetr, doublecomplex *m, doublecomplex *f)
{
    // This function handles the construction of the elemental matrices
    // for each tetrahedron. The matrices m and f are defined as:
    // m(i,j) = integral_over_tetrahedron(curl(N_i) . curl(N_j)),
    // f(i,j) = integral_over_tetrahedron(N_i . N_j),
    // where N_i is the vectorial basis function for edge elements.
    // m and f are 6x6 real matrices, but are stored here as complex.
    // INPUT:
    // point3d* points: the array containing the information on the mesh nodes.
    // edge* edges: the array containing the information on the mesh edges.
    // tetra* tetr: a pointer to the tetrahedron over which compute the elemental matrices.
    // doublecomplex* m, doublecomplex* f: arrays of length 9 that will be overwritten with
    //      the elemental matrices.
    // OUTPUT:
    // on exit, the arrays m[] and f[] contain the elemental matrices for the given tetrahedron.

    int i, j;
    double volume, l[6];
    double be[4], ce[4], de[4], cf[6][6];
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;

    // Copies the data from the arrays
    x1 = points[tetr->p[0]].x;
    x2 = points[tetr->p[1]].x;
    x3 = points[tetr->p[2]].x;
    x4 = points[tetr->p[3]].x;
    y1 = points[tetr->p[0]].y;
    y2 = points[tetr->p[1]].y;
    y3 = points[tetr->p[2]].y;
    y4 = points[tetr->p[3]].y;
    z1 = points[tetr->p[0]].z;
    z2 = points[tetr->p[1]].z;
    z3 = points[tetr->p[2]].z;
    z4 = points[tetr->p[3]].z;

    volume = tetr->volume;

    for(i = 0; i < 6; i++)
        l[i] = edges[tetr->e[i]].length;

    // m matrix construction
    for(i = 0; i < 6; i++)
        for(j = 0; j < 6; j++)
        {
            m[6*i + j].re = l[i]*l[j]*l[5-i]*l[5-j]
                      *dot_product_3d(&edge_versors[tetr->e[5-i]], &edge_versors[tetr->e[5-j]])
                      *(tetr->o[5-i]*tetr->o[5-j])/(9*volume);
        }

    // Coefficients computation
    be[0] = y2*z3 - y3*z2 - y2*z4 + y4*z2 + y3*z4 - y4*z3;
    be[1] = y3*z1 - y1*z3 + y1*z4 - y4*z1 - y3*z4 + y4*z3;
    be[2] = y1*z2 - y2*z1 - y1*z4 + y4*z1 + y2*z4 - y4*z2;
    be[3] = y2*z1 - y1*z2 + y1*z3 - y3*z1 - y2*z3 + y3*z2;

    ce[0] = x2*z3 - x3*z2 - x2*z4 + x4*z2 + x3*z4 - x4*z3;
    ce[1] = x3*z1 - x1*z3 + x1*z4 - x4*z1 - x3*z4 + x4*z3;
    ce[2] = x1*z2 - x2*z1 - x1*z4 + x4*z1 + x2*z4 - x4*z2;
    ce[3] = x2*z1 - x1*z2 + x1*z3 - x3*z1 - x2*z3 + x3*z2;

    de[0] = x2*y3 - x3*y2 - x2*y4 + x4*y2 + x3*y4 - x4*y3;
    de[1] = x3*y1 - x1*y3 + x1*y4 - x4*y1 - x3*y4 + x4*y3;
    de[2] = x1*y2 - x2*y1 - x1*y4 + x4*y1 + x2*y4 - x4*y2;
    de[3] = x2*y1 - x1*y2 + x1*y3 - x3*y1 - x2*y3 + x3*y2;

    for(i=0;i<4;i++)
        for(j=0;j<4;j++)
            {
                cf[i][j] = be[i]*be[j] + ce[i]*ce[j] + de[i]*de[j];
            }

    // f matrix construction
    f[0].re = 2*l[0]*l[0]*(cf[1][1] - cf[0][1] + cf[0][0])/(volume*720);
    f[6].re = f[1].re = l[0]*l[1]*(2*cf[1][2] - cf[1][0] - cf[0][2] + cf[0][0])/(volume*720);
    f[12].re = f[2].re = l[0]*l[2]*(2*cf[1][3] - cf[1][0] - cf[0][3] + cf[0][0])/(volume*720);
    f[18].re = f[3].re = l[0]*l[3]*(cf[1][2] - cf[1][1] - 2*cf[0][2] + cf[0][1])/(volume*720);
    f[24].re = f[4].re = l[0]*l[4]*(cf[1][1] - cf[1][3] - cf[0][1] + 2*cf[0][3])/(volume*720);
    f[30].re = f[5].re = l[0]*l[5]*(cf[1][3] - cf[1][2] - cf[0][3] + cf[0][2])/(volume*720);
    f[7].re = 2*l[1]*l[1]*(cf[2][2] - cf[0][2] + cf[0][0])/(volume*720);
    f[13].re = f[8].re = l[1]*l[2]*(2*cf[2][3] - 2*cf[0][2] - 2*cf[0][3] + cf[0][0])/(volume*720);
    f[19].re = f[9].re = l[1]*l[3]*(cf[2][2] - cf[1][2] - cf[0][2] + 2*cf[0][1])/(volume*720);
    f[25].re = f[10].re = l[1]*l[4]*(cf[1][2] - cf[2][3] - cf[0][1] + cf[0][3])/(volume*720);
    f[31].re = f[11].re = l[1]*l[5]*(cf[0][2] - cf[2][2] - 2*cf[0][3] + cf[2][3])/(volume*720);
    f[14].re = 2*l[2]*l[2]*(cf[3][3] - cf[0][3] + cf[0][0])/(volume*720);
    f[20].re = f[15].re = l[2]*l[3]*(cf[2][3] - cf[1][3] - cf[0][2] + cf[0][1])/(volume*720);
    f[26].re = f[16].re = l[2]*l[4]*(cf[1][3] - cf[3][3] - 2*cf[0][1] + cf[0][3])/(volume*720);
    f[32].re = f[17].re = l[2]*l[5]*(cf[3][3] - cf[2][3] - cf[0][3] + 2*cf[0][2])/(volume*720);
    f[21].re = 2*l[3]*l[3]*(cf[2][2] - cf[1][2] + cf[1][1])/(volume*720);
    f[27].re = f[22].re = l[3]*l[4]*(cf[1][2] - 2*cf[2][3] - cf[1][1] + cf[1][3])/(volume*720);
    f[33].re = f[23].re = l[3]*l[5]*(cf[2][3] - cf[2][2] - 2*cf[1][3] + cf[1][2])/(volume*720);
    f[28].re = 2*l[4]*l[4]*(cf[1][1] - cf[1][3] + cf[3][3])/(volume*720);
    f[34].re = f[29].re = l[4]*l[5]*(cf[1][3] - 2*cf[1][2] - cf[3][3] + cf[2][3])/(volume*720);
    f[35].re = 2*l[5]*l[5]*(cf[3][3] - 2*cf[2][3] + cf[2][2])/(volume*720);

    // The matrices are real, so the imaginary part is set to zero.
    for(i=0;i<36;i++)
        m[i].i = f[i].i = 0;

}

void add_to_matrix_3d(int *ia, int *ja, doublecomplex *a, doublecomplex *m, tetra* tetr, int row_offset, int column_offset, double coefficient, int *pos)
{
    // Adds the elemental matrix m to the arrays containing
    // the matrix entries for A, managing row/column offset and a coefficient.
    // The 6x6 matrix m is split and its elements are added to the A matrix
    // arrays in the rows and columns corresponding to the edges of the
    // tetrahedron, after being multiplied by the coefficient and the orientations of the edges.
    // Specifying a nonzero offset allows the construction of extended matrices (for the Lorentz model).
    // INPUT:
    // int *ia, int *ja, doublecomplex *a: the arrays of lenght *pos containing the information on
    //      the matrix A, or empty arrays
    // doublecomplex* m: an array of size 9 containing the elemental matrix to be added
    // tetra* tetr: a pointer to the tetrahedron corresponding to the elemental matrix m
    // int row_offset: the offset on the row indices in which the m matrix is inserted
    // int column_offset: the offset on the column indices in which the m matrix is inserted.
    // double coefficient: a real coefficient that multiplies the matrix m before insertion.
    // int* pos: a pointer to the first empty position in the ia, ja, a arrays. Should be set to 0 before the
    //       first call to this function, and passed unchanged to successive calls. After the last
    //       call, the value of *pos is the size of the arrays.
    // OUTPUT:
    // int *ia, int *ja, doublecomplex *a: the 9 entries from the elemental matrix m, multiplied by
    //       the coefficient and the orientations, are added to these arrays starting from *pos.
    // int* pos: on exit, *pos is the index of the first empty position in ia, ja and a, i.e. the new number
    //       of entries.
    int j, k;
    for (j = 0; j < 6; j++)
        {
            for (k = 0; k < 6; k++)
            {

                int orient = tetr->o[j]*tetr->o[k];

                ia[*pos] = tetr->e[j] + row_offset;
                ja[*pos] = tetr->e[k] + column_offset;

                a[*pos].re = coefficient*orient*m[6*j + k].re;
                a[*pos].i = coefficient*orient*m[6*j + k].i;
                (*pos)++;
            }
        }
}
