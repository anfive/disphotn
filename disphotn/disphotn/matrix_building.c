void matrix_construction(point2d *points, edge *edges, triangle *trie, doublecomplex *m, doublecomplex *f)
{
    // This function handles the construction of the elemental matrices
    // for each triangle. The matrices m and f are defined as:
    // m(i,j) = integral_over_triangle(curl(N_i) . curl(N_j)),
    // f(i,j) = integral_over_triangle(N_i . N_j),
    // where N_i is the vectorial basis function for edge elements.
    // m and f are 3x3 real matrices, but are stored here as complex.
    // INPUT:
    // point2d* points: the array containing the information on the mesh nodes.
    // edge* edges: the array containing the information on the mesh edges.
    // triangle* trie: a pointer to the triangle over which compute the elemental matrices.
    // doublecomplex* m, doublecomplex* f: arrays of length 9 that will be overwritten with
    //      the elemental matrices.
    // OUTPUT:
    // on exit, the arrays m[] and f[] contain the elemental matrices for the given triangle.

    int i;
    double area, ln1, ln2, ln3;
    double x1, x2, x3, y1, y2, y3;

    // Copies the data from the arrays
    x1 = points[trie->p[0]].x;
    x2 = points[trie->p[1]].x;
    x3 = points[trie->p[2]].x;
    y1 = points[trie->p[0]].y;
    y2 = points[trie->p[1]].y;
    y3 = points[trie->p[2]].y;

    area = trie->area;
    ln1 = edges[trie->e[0]].length;
    ln2 = edges[trie->e[1]].length;
    ln3 = edges[trie->e[2]].length;

    // m matrix construction
    m[0].re = ln1*ln1/area;
    m[3].re = m[1].re = ln1*ln2/area;
    m[2].re = m[6].re = ln1*ln3/area;
    m[4].re = ln2*ln2/area;
    m[5].re = m[7].re = ln2*ln3/area;
    m[8].re = ln3*ln3/area;

    // f matrix construction
    f[0].re = (ln1*ln1*(x1*x1 + x2*x2 + x1*(x2 - 3*x3) - 3*x2*x3 + 3*x3*x3 + y1*y1 + y1*y2 + y2*y2 - 3*y1*y3 - 3*y2*y3 + 3*y3*y3))/(24*area);
    f[3].re = f[1].re = -((ln1* ln2*(x1*x1 - x2*x2 + x1*(x2 - 3*x3) + x2*x3 + x3*x3 + y1*y1 + y1*y2 - y2*y2 - 3*y1*y3 + y2*y3 + y3*y3))/(24*area));
    f[2].re = f[6].re = (ln1* ln3*(x1*x1 - x2*x2 + 3*x2*x3 - x3*x3 - x1*(x2 + x3) + y1*y1 - y1*y2 - y2*y2 - y1*y3 + 3*y2*y3 - y3*y3))/(24*area);
    f[4].re = (ln2*ln2*(3*x1*x1 + x2*x2 + x2*x3 + x3*x3 - 3*x1*(x2 + x3) + 3*y1*y1 - 3*y1*y2 + y2*y2 - 3*y1*y3 + y2*y3 + y3*y3))/(24*area);
    f[5].re = f[7].re = -((ln2* ln3*(x1*x1 + x2*x2 + x2*x3 - x3*x3 + x1*(-3*x2 + x3) + y1*y1 - 3*y1*y2 + y2*y2 + y1*y3 + y2*y3 - y3*y3))/(24*area));
    f[8].re = (ln3*ln3*(x1*x1 + 3*x2*x2 - 3*x2*x3 + x3*x3 + x1*(-3*x2 + x3) + y1*y1 - 3*y1*y2 + 3*y2*y2 + y1*y3 - 3*y2*y3 + y3*y3))/(24*area);

    // The matrices are real, so the imaginary part is set to zero.
    for(i=0;i<9;i++)
        m[i].i = f[i].i = 0;
}

void add_to_matrix(int *ia, int *ja, doublecomplex *a, doublecomplex *m, triangle* trie, int row_offset, int column_offset, double coefficient, int *pos)
{
    // Adds the elemental matrix m to the arrays containing
    // the matrix entries for A, managing row/column offset and a coefficient.
    // The 3x3 matrix m is split and its elements are added to the A matrix
    // arrays in the rows and columns corresponding to the edges of the
    // triangle, after being multiplied by the coefficient and the orientations of the edges.
    // Specifying a nonzero offset allows the construction of extended matrices (for the Lorentz model).
    // INPUT:
    // int *ia, int *ja, doublecomplex *a: the arrays of lenght *pos containing the information on
    //      the matrix A, or empty arrays
    // doublecomplex* m: an array of size 9 containing the elemental matrix to be added
    // triangle* trie: a pointer to the triangle corresponding to the elemental matrix m
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
    // Scans the 3 rows...
    for (j = 0; j < 3; j++)
        {
            // ... and the 3 columns of m.
            for (k = 0; k < 3; k++)
            {
                // Computes the orientation coefficient. This can
                // be either 1 or -1 (always 1 on diagonal).
                int orient = trie->o[j]*trie->o[k];

                // Adds a new entry to the row and column corresponding to the
                // current edges, with the provided offset.
                ia[*pos] = trie->e[j] + row_offset;
                ja[*pos] = trie->e[k] + column_offset;

                // Adds a new entry with the value, multiplied by the coefficient
                // and the orientation.
                a[*pos].re = coefficient*orient*m[3*j + k].re;
                a[*pos].i = coefficient*orient*m[3*j + k].i;

                // Updates the index of the first free position
                (*pos)++;
            }
        }
}
