// This file contains the function used to manipulate the geometric information of the problem.

int identify_periodic_boundaries_3d(point3d *points, edge *edges, int *boundary[], int *nbound, int b1, int b2, point3d *a, byte *boundary_orientations)
{
    // Identifies the corresponding edges on two boundaries. Used to apply periodic boundary conditions.
    // INPUT:
    // point3d* points: the array containing the information on the nodes.
    // edge* edges: the array containing the information on the edges.
    // int* boundary[]: an array of pointers. The i-th element of the array points to another array containing the
    //                  number of the edges on the i-th boundary.
    // int* nbound: an array containing the number of elements on each boundary.
    // int b1: the index (within boundary[]) of the source boundary.
    // int b2: the index (within boundary[]) of the destination boundary.
    // point3d* a: the Bravais lattice vector corresponding to the pair of boundaries, expressed as a 3-elements vector.
    // OUTPUT:
    // The function returns 0 on success. A positive return value i indicates that the i-th edge in the Source boundary has
    // no corresponding edge.
    // int* boundary[]: the array corresponding to the Destination boundary (specified by b2) is ordered such that
    //                  each element correspond to the related edge in the Source array.
    // byte* boundary_orientations: an array of the orientations (+1 or -1) of the corresponding edges.

    int *ordered, i, j;
    point3d p1, p2;
    point3d exp1, exp2, dist1, dist2;

    // Auxiliary array that will contain the ordered list of edges on the Destination boundary.
    ordered = (int*)malloc(nbound[b1]*sizeof(int));



    int *bound1, *bound2;
    bound1 = boundary[b1];
    bound2 = boundary[b2];

    for (i = 0; i < nbound[b1]; i++)
    {

        ordered[i] = -1;;

        p1 = points[edges[bound1[i]].p[0]];
        p2 = points[edges[bound1[i]].p[1]];

        // For each edge of the Source boundary, looks into the Destination array
        // looking for a corresponding edge.

        // Computes the expected position of the nodes of the corresponding edge
        exp1 = vector_subtract_3d(p1, *a);
        exp2 = vector_subtract_3d(p2, *a);

        for (j = 0; j < nbound[b2] && ordered[i] == -1; j++)
        {

            p1 = points[edges[bound2[j]].p[0]];
            p2 = points[edges[bound2[j]].p[1]];


            // Computes the difference between the expected position of the corresponding edge
            // and the position of the edge under examination.
            dist1 = vector_subtract_3d(exp1, p1);
            dist2 = vector_subtract_3d(exp2, p2);

            // If this is the corresponding edge, the difference is zero (within rounding errors)
            if (near_zero(dist1.x) && near_zero(dist1.y) && near_zero(dist1.z) && near_zero(dist2.x) && near_zero(dist2.y) && near_zero(dist2.z))
            {
                // The edge is the corresponding one. Its index is saved in the auxiliary array.
                ordered[i] = bound2[j];
                boundary_orientations[i] = 1;
            }

            // Computes the difference again. but this time inverting the edge.
            dist1 = vector_subtract_3d(exp1, p2);
            dist2 = vector_subtract_3d(exp2, p1);

            if (near_zero(dist1.x) && near_zero(dist1.y) && near_zero(dist1.z) && near_zero(dist2.x) && near_zero(dist2.y) && near_zero(dist2.z))
            {
                // The edge is the corresponding one, but inverted. Its index is saved in the auxiliary array,
                // and its orientation is set to -1.
                ordered[i] = bound2[j];
                boundary_orientations[i] = -1;
            }

        }

        if (ordered[i] == -1)
        {
            // No edge in the Destination boundary corresponds to the one on the Source boundary.

            return bound1[i];
        }
    }

    // The correspondence identification is completed.

    // Frees the old (unordered) boundary array.
    free(boundary[b2]);

    // Stores the ordered boundary array.
    boundary[b2] = ordered;

    return -1;
}

int identify_edges_3d(edge **edges, int *nedges, tetra *tetras, int *ntetras)
{
    // Inspects the mesh tetrahedra and assigns an unique edge index to each edge (group of 2 connected nodes).
    // INPUT:
    // tetra *tetras: an array containing the information on the mesh tetrahedra. Only the node fields need to be
    //                  initialized (the first 4 int fields).
    // int *ntetras: the number of elements in the tetras array.
    // OUTPUT:
    // The function returns 0 in case of success. It returns an error code > 0 in case of error.
    // edge **edges: a pointer to an array containing the newly created edge information.
    // int *nedges: the number of edge created
    // tetra *tetras: the edge (int e) and orientation (int o) fields of each tetrahedron contain now the edge
    //                  numbers and their orientations.

    *nedges = 0;
    edge *edg;
    edg = malloc(3*(*ntetras)*sizeof(edge));

    int edge_order[6][2], i, j, k;
    int *p, *e, check = 0;
    bool *o;

    // The convention for numbering edges is the following:
    // Edge 0: from node 1 to node 2.
    edge_order[0][0] = 0; edge_order[0][1] = 1;
    // Edge 1: from node 1 to node 3.
    edge_order[1][0] = 0; edge_order[1][1] = 2;
    // Edge 2: from node 1 to node 4.
    edge_order[2][0] = 0; edge_order[2][1] = 3;
    // Edge 3: from node 2 to node 3.
    edge_order[3][0] = 1; edge_order[3][1] = 2;
    // Edge 4: from node 4 to node 2.
    edge_order[4][0] = 3; edge_order[4][1] = 1;
    // Edge 5: from node 3 to node 4.
    edge_order[5][0] = 2; edge_order[5][1] = 3;


    for (i = 0; i < *ntetras; i++)
    {
        // check is an int variable, initialized to zero. When an edge is found,
        // check is incremented by one. At the end of the cycle, if check < 6 then
        // some edges have not been found and must be created.
        check = 0;

        // Retrieves the pointer to the tetrahedron data.
        p = tetras[i].p;

        e = tetras[i].e;
        o = tetras[i].o;

        // Initializes the edge indices.
        for(j = 0; j < 6; j++)
            e[j] = -1;

        // Scans the array of the existing edges, to
        // see if the edges for this triangle have already been created.
        for (j = 0; j < *nedges && check < 6; j++)
        {
            for (k = 0; k < 6; k++)
            {
                if (e[k] == -1)
                {
                    // The k-th edge has not been found yet.

                    if (p[edge_order[k][0]] == edg[j].p[0] && p[edge_order[k][1]] == edg[j].p[1])
                    {
                        // The j-th edge corresponds to the edge we are looking for. Its index
                        // is saved in the triangle data, with the orientation.

                        e[k] = j;
                        o[k] = 1;
                        check++;
                    }
                    else if (p[edge_order[k][0]] == edg[j].p[1] && p[edge_order[k][1]] == edg[j].p[0])
                    {
                        // We found a corresponding edge, but with opposite orientation.

                        e[k] = j;
                        o[k] = -1;
                        check++;
                    }

                }
            }

        }
        if (check < 6)
        {
            // Some edges are missing
            for (k = 0; k < 6; k++)
            {
                if (e[k] == -1)
                {
                    // If the edge is missing, a new entry in the edg[] array is created
                    // with the missing edge information.
                    edg[*nedges].p[0] = p[edge_order[k][0]];
                    edg[*nedges].p[1] = p[edge_order[k][1]];
                    e[k] = *nedges;
                    o[k] = 1;
                    (*nedges)++;
                }
            }
        }
    }

    // Resizes the edges array to the actual size.
    *edges = realloc(edg, *nedges*sizeof(edge));
    if (*edges == NULL)
    {
        free(edg);
        return -1;
    }

    return 0;


}

double tetra_volume(point3d *p1, point3d *p2, point3d *p3, point3d *p4)
{
    // Computes the volume of a tetrahedron using the determinant method.
    // Pretty ugly, in my opinion.
    // INPUT:
    // point3d *p1, point3d *p2, point3d *p3, point3d *p4: pointers to
    //      the vertices of the tetrahedron.
    // OUTPUT:
    // The function returns the volume of the tetrahedron as a double precision number.
    return( -(p1->x*p2->z*p3->y - p1->x*p2->y*p3->z + p1->y*p2->x*p3->z - p1->y*p2->z*p3->x - p1->z*p2->x*p3->y + p1->z*p2->y*p3->x + p1->x*p2->y*p4->z - p1->x*p2->z*p4->y
            - p1->y*p2->x*p4->z + p1->y*p2->z*p4->x + p1->z*p2->x*p4->y - p1->z*p2->y*p4->x - p1->x*p3->y*p4->z + p1->x*p3->z*p4->y + p1->y*p3->x*p4->z - p1->y*p3->z*p4->x
            - p1->z*p3->x*p4->y + p1->z*p3->y*p4->x + p2->x*p3->y*p4->z - p2->x*p3->z*p4->y - p2->y*p3->x*p4->z + p2->y*p3->z*p4->x + p2->z*p3->x*p4->y - p2->z*p3->y*p4->x)/6);

}


bool is_in_tetra(point3d *p, point3d *p1, point3d *p2, point3d *p3, point3d *p4)
{
    // Determines if a certain point is inside a tetrahedron.
    // INPUT:
    // point3d* p: pointer to a 3D point.
    // point3d* p1, point3d* p2, point3d* p3, point3d* p4: pointers to the vertices of the
    //      tetrahedron.
    // OUTPUT:
    // The function returns true (-1) if the point p is inside the tetrahedron defined
    // by *p1, *p2, *p3 and *p4, 0 otherwise.

    double d, q;

    // To determine if the point is inside the tetrahedron, it computes the signed volumes
    // of the tetrahedra obtained substituting one of the vertices with p. If all the volumes
    // and the signed area of the tetrahedron have the same sign, then the point is inside the tetrahedron.
    // However, due to rounding errors, if the point lies on an edge, the signed volume may be a
    // small but non-zero number, of either sign. So we must allow a certain tolerance on the volume
    // signs.

//    d = (p1->x*p2->z*p3->y - p1->x*p2->y*p3->z + p1->y*p2->x*p3->z - p1->y*p2->z*p3->x - p1->z*p2->x*p3->y + p1->z*p2->y*p3->x + p1->x*p2->y*p4->z - p1->x*p2->z*p4->y
//        - p1->y*p2->x*p4->z + p1->y*p2->z*p4->x + p1->z*p2->x*p4->y - p1->z*p2->y*p4->x - p1->x*p3->y*p4->z + p1->x*p3->z*p4->y + p1->y*p3->x*p4->z - p1->y*p3->z*p4->x
//        - p1->z*p3->x*p4->y + p1->z*p3->y*p4->x + p2->x*p3->y*p4->z - p2->x*p3->z*p4->y - p2->y*p3->x*p4->z + p2->y*p3->z*p4->x + p2->z*p3->x*p4->y - p2->z*p3->y*p4->x);
//    q = (p->x*p2->z*p3->y - p->x*p2->y*p3->z + p->y*p2->x*p3->z - p->y*p2->z*p3->x - p->z*p2->x*p3->y + p->z*p2->y*p3->x + p->x*p2->y*p4->z - p->x*p2->z*p4->y
//        - p->y*p2->x*p4->z + p->y*p2->z*p4->x + p->z*p2->x*p4->y - p->z*p2->y*p4->x - p->x*p3->y*p4->z + p->x*p3->z*p4->y + p->y*p3->x*p4->z - p->y*p3->z*p4->x
//        - p->z*p3->x*p4->y + p->z*p3->y*p4->x + p2->x*p3->y*p4->z - p2->x*p3->z*p4->y - p2->y*p3->x*p4->z + p2->y*p3->z*p4->x + p2->z*p3->x*p4->y - p2->z*p3->y*p4->x);

    d = tetra_volume(p1, p2, p3, p4);
    q = tetra_volume(p, p2, p3, p4);

    if (d*q >= - NEARZERO_STRICT_TOLERANCE)
    {
//        q = (p1->x*p->z*p3->y - p1->x*p->y*p3->z + p1->y*p->x*p3->z - p1->y*p->z*p3->x - p1->z*p->x*p3->y + p1->z*p->y*p3->x + p1->x*p->y*p4->z - p1->x*p->z*p4->y
//            - p1->y*p->x*p4->z + p1->y*p->z*p4->x + p1->z*p->x*p4->y - p1->z*p->y*p4->x - p1->x*p3->y*p4->z + p1->x*p3->z*p4->y + p1->y*p3->x*p4->z - p1->y*p3->z*p4->x
//            - p1->z*p3->x*p4->y + p1->z*p3->y*p4->x + p->x*p3->y*p4->z - p->x*p3->z*p4->y - p->y*p3->x*p4->z + p->y*p3->z*p4->x + p->z*p3->x*p4->y - p->z*p3->y*p4->x);
        q = tetra_volume(p1, p, p3, p4);

        if (d*q >= - NEARZERO_STRICT_TOLERANCE)
        {
//            q = (p1->x*p2->z*p->y - p1->x*p2->y*p->z + p1->y*p2->x*p->z - p1->y*p2->z*p->x - p1->z*p2->x*p->y + p1->z*p2->y*p->x + p1->x*p2->y*p4->z - p1->x*p2->z*p4->y
//                - p1->y*p2->x*p4->z + p1->y*p2->z*p4->x + p1->z*p2->x*p4->y - p1->z*p2->y*p4->x - p1->x*p->y*p4->z + p1->x*p->z*p4->y + p1->y*p->x*p4->z - p1->y*p->z*p4->x
//                - p1->z*p->x*p4->y + p1->z*p->y*p4->x + p2->x*p->y*p4->z - p2->x*p->z*p4->y - p2->y*p->x*p4->z + p2->y*p->z*p4->x + p2->z*p->x*p4->y - p2->z*p->y*p4->x);
            q = tetra_volume(p1, p2, p, p4);

            if (d*q >= - NEARZERO_STRICT_TOLERANCE)
            {
//                q = (p1->x*p2->z*p3->y - p1->x*p2->y*p3->z + p1->y*p2->x*p3->z - p1->y*p2->z*p3->x - p1->z*p2->x*p3->y + p1->z*p2->y*p3->x + p1->x*p2->y*p->z - p1->x*p2->z*p->y
//                    - p1->y*p2->x*p->z + p1->y*p2->z*p->x + p1->z*p2->x*p->y - p1->z*p2->y*p->x - p1->x*p3->y*p->z + p1->x*p3->z*p->y + p1->y*p3->x*p->z - p1->y*p3->z*p->x
//                    - p1->z*p3->x*p->y + p1->z*p3->y*p->x + p2->x*p3->y*p->z - p2->x*p3->z*p->y - p2->y*p3->x*p->z + p2->y*p3->z*p->x + p2->z*p3->x*p->y - p2->z*p3->y*p->x);
                q = tetra_volume(p1, p2, p3, p);

                if (d*q >= - NEARZERO_STRICT_TOLERANCE)
                {

                    return -1;
                }
            }

        }

    }

    return 0;
}

void interpolate_tetra(point3d *p, point3d *points, tetra *tetr, doublecomplex *field, field3d *f)
{
    // Interpolates the field on a tetrahedron using the edge-elements
    // basis functions.
    // INPUT:
    // point2d* p: a pointer to the 3D point in which the interpolation
    //      must be carried on.
    // point2d* p1, point2d* p2, point2d* p3, point2d* p4: pointers to the vertices of
    //      the tetrahedron.
    // doublecomplex* field: an array of size 6 containing the value of the field
    //      on each edge.
    // field3d* f: pointer to a field3d.
    // OUTPUT:
    // field3d* f: on exit, the vector value of the field at the interpolation point.
    int i, j;

    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    double L[4], gradL[4][3], N[3];
    double be[4], ce[4], de[4], lengths[6], volume;
    point3d *pt[4];

    int edge_i[6][2];

    edge_i[0][0] = 0;
    edge_i[0][1] = 1;
    edge_i[1][0] = 0;
    edge_i[1][1] = 2;
    edge_i[2][0] = 0;
    edge_i[2][1] = 3;
    edge_i[3][0] = 1;
    edge_i[3][1] = 2;
    edge_i[4][0] = 3;
    edge_i[4][1] = 1;
    edge_i[5][0] = 2;
    edge_i[5][1] = 3;

    // Computes the edge lengths
    for(i = 0; i < 6; i++)
    {
        pt[0] = &points[tetr->p[edge_i[i][0]]];
        pt[1] = &points[tetr->p[edge_i[i][1]]];

        lengths[i] = sqrt(square(pt[0]->x - pt[1]->x) + square(pt[0]->y - pt[1]->y) + square(pt[0]->z - pt[1]->z));

    }

    // Computes the volume
    volume = tetra_volume(&points[tetr->p[0]],&points[tetr->p[1]],&points[tetr->p[2]],&points[tetr->p[3]]);


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

    // Computes coefficients
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

    // Computes the Linear interpolation functions and their
    // gradients.
    for(i = 0; i < 4; i++)
    {
        pt[0] = &points[tetr->p[0]];
        pt[1] = &points[tetr->p[1]];
        pt[2] = &points[tetr->p[2]];
        pt[3] = &points[tetr->p[3]];

        pt[i] = p;

        L[i] = tetra_volume(pt[0], pt[1], pt[2], pt[3])/volume;


        gradL[i][0] = -be[i]/(6*volume);
        gradL[i][1] =  ce[i]/(6*volume);
        gradL[i][2] = -de[i]/(6*volume);
    }

    f->x.re = f->x.i = f->y.re = f->y.i = f->z.re = f->z.i = 0;
    for(i = 0; i < 6; i++)
    {
        for(j = 0; j < 3; j++)
        {
            // Computes the j-th component of the i-th basis function.
            N[j] = lengths[i]*(L[edge_i[i][0]]*gradL[edge_i[i][1]][j] - L[edge_i[i][1]]*gradL[edge_i[i][0]][j]);

        }

        // Adds the value of the field relative to this edge to the total.
        f->x.re += N[0]*field[i].re;
        f->x.i  += N[0]*field[i].i;
        f->y.re += N[1]*field[i].re;
        f->y.i  += N[1]*field[i].i;
        f->z.re += N[2]*field[i].re;
        f->z.i  += N[2]*field[i].i;
    }
}

int interpolate_3d(point3d *grid, int gridn, tetra *mtetras, int ntetras, point3d *mpoints, doublecomplex *mfield, field3d *efield)
{
    // Interpolates the electric field, given as edge-elements coefficients,
    // on a provided three-dimensional grid.
    // INPUT:
    // point3d* grid: a vector of size gridn containing the interpolation points.
    // int gridn: the number of points in the interpolation grid.
    // tetra *mtetras: the array containing the tetrahedra data for the mesh.
    // int ntetras: the number of elements in the mtetras[] array.
    // point3d *mpoints: the array containing the mesh nodes data.
    // doublecomplex *mfield: the field expressed as edge-element coefficients.
    // field3d *efield: pointer to an array of size gridn.
    // OUTPUT:
    // The function return 0 on success, and 1 if the grid is empty.
    // field3d *efield: On exit, the array of the interpolated field on the points
    //      contained in grid[].

    int i, j, k;
    bool found;
    point3d *p, *p1, *p2, *p3, *p4;
    doublecomplex *field;

    // Check if the grid is not empty.
    if (gridn <= 0)
        return 1;

    // Allocates an auxiliary vector to store the field values
    // for the current tetrahedron
    field = malloc(6*sizeof(doublecomplex));

    for (i = 0; i < gridn; i++)
    {
        // Interpolates on each grid point.
        p = &grid[i];

        // Looks for the tetrahedron containing the current point.
        found = 0;
        for (j = 0; j < ntetras && !found; j++)
        {
            p1 = &mpoints[mtetras[j].p[0]];
            p2 = &mpoints[mtetras[j].p[1]];
            p3 = &mpoints[mtetras[j].p[2]];
            p4 = &mpoints[mtetras[j].p[3]];

            if (is_in_tetra(p, p1, p2, p3, p4))
            {
                // If the point is in the current tetrahedron, saves the values of
                // the field on the edges (with the relative orientations)
                // in the auxiliary array field[], that will be then passed to interpolate_tetra.
                for (k = 0; k < 6; k++)
                {
                    field[k].re = mfield[mtetras[j].e[k]].re*mtetras[j].o[k];
                    field[k].i = mfield[mtetras[j].e[k]].i*mtetras[j].o[k];
                }
                // Interpolates the field in the tetrahedron.
                interpolate_tetra(p, mpoints, &mtetras[j], field, &efield[i]);

                found = -1;
            }
        }

        if (!found)
        {
            // No tetrhedron containing the point has been found. The field is set to zero arbitrarily.
            printf("WARNING: Interpolation point (%g, %g, %g) outside of mesh.\n", p->x, p->y, p->z);
            efield[i].x.re = efield[i].x.i = efield[i].y.re = efield[i].y.i = efield[i].z.re = efield[i].z.i =0;
        }

    }

    return 0;

}


