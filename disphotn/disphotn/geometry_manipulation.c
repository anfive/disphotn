// This file contains the function used to manipulate the geometric information of the problem.

int identify_periodic_boundaries(point2d *points, edge *edges, int *boundary[], int *nbound, int b1, int b2, point2d *a, byte *boundary_orientations)
{
    // Identifies the corresponding edges on two boundaries. Used to apply periodic boundary conditions.
    // INPUT:
    // point2d* points: the array containing the information on the nodes.
    // edge* edges: the array containing the information on the edges.
    // int* boundary[]: an array of pointers. The i-th element of the array points to another array containing the
    //                  number of the edges on the i-th boundary.
    // int* nbound: an array containing the number of elements on each boundary.
    // int b1: the index (within boundary[]) of the source boundary.
    // int b2: the index (within boundary[]) of the destination boundary.
    // point2d* a: the Bravais lattice vector corresponding to the pair of boundaries, expressed as a 2 element vector.
    // OUTPUT:
    // The function returns 0 on success. A positive return value i indicates that the i-th edge in the Source boundary has
    // no corresponding edge.
    // int* boundary[]: the array corresponding to the Destination boundary (specified by b2) is ordered such that
    //                  each element correspond to the related edge in the Source array.
    // byte* boundary_orientations: an array of the orientations (+1 or -1) of the corresponding edges.

    int *ordered, i, j;
    point2d p1, p2;
    point2d exp1, exp2, dist1, dist2;

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
        exp1 = vector_subtract(p1, *a);
        exp2 = vector_subtract(p2, *a);

        for (j = 0; j < nbound[b2] && ordered[i] == -1; j++)
        {

            p1 = points[edges[bound2[j]].p[0]];
            p2 = points[edges[bound2[j]].p[1]];


            // Computes the difference between the expected position of the corresponding edge
            // and the position of the edge under examination.
            dist1 = vector_subtract(exp1, p1);
            dist2 = vector_subtract(exp2, p2);

            // If this is the corresponding edge, the difference is zero (within rounding errors)
            if (near_zero(dist1.x) && near_zero(dist1.y) && near_zero(dist2.x) && near_zero(dist2.y))
            {
                // The edge is the corresponding one. Its index is saved in the auxiliary array.
                ordered[i] = bound2[j];
                boundary_orientations[i] = 1;
            }

            // Computes the difference again. but this time inverting the edge.
            dist1 = vector_subtract(exp1, p2);
            dist2 = vector_subtract(exp2, p1);

            if (near_zero_tol(dist1.x) && near_zero_tol(dist1.y) && near_zero_tol(dist2.x) && near_zero_tol(dist2.y))
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

int identify_edges(edge **edges, int *nedges, triangle *tries, int *ntries)
{
    // Inspects the mesh triangles and assigns an unique edge index to each edge (group of 2 connected nodes).
    // INPUT:
    // triangle *tries: an array containing the information on the mesh triangle. Only the node fields need to be
    //                  initialized (the first 3 int fields).
    // int *ntries: the number of elements in the tries array.
    // OUTPUT:
    // The function returns 0 in case of success. It returns an error code > 0 in case of error.
    // edge **edges: a pointer to an array containing the newly created edge information.
    // int *nedges: the number of edge created
    // triangle *tries: the edge (int e) and orientation (int o) fields of each triangle contain now the edge
    //                  numbers and their orientations.

    *nedges = 0;
    edge *edg;
    edg = (edge*)malloc(3*(*ntries)*sizeof(edge));

    int edge_order[4], i, j, k;
    int *p, *e, check = 0;
    bool *o;

    // The convention for numbering edges is the following:
    // Edge 0: from node 0 to node 1.
    // Edge 1: from node 1 to node 2.
    // Edge 2: from node 2 to node 0.
    edge_order[0] = 0;
    edge_order[1] = 1;
    edge_order[2] = 2;
    edge_order[3] = 0;

    for (i = 0; i < *ntries; i++)
    {
        // check is an int variable, initialized to zero. When an edge is found,
        // check is incremented by one. At the end of the cycle, if check < 3 then
        // some edges have not been found and must be created.
        check = 0;

        // Retrieves the pointer to the triangle data.
        p = tries[i].p;

        e = tries[i].e;
        o = tries[i].o;

        // Initializes the edge indices.
        e[0] = e[1] = e[2] = -1;

        // Scans the array of the existing edges, to
        // see if the edges for this triangle have already been created.
        for (j = 0; j < *nedges && check < 3; j++)
        {
            for (k = 0; k < 3; k++)
            {
                if (e[k] == -1)
                {
                    // The k-th edge has not been found yet.

                    if (p[edge_order[k]] == edg[j].p[0] && p[edge_order[k + 1]] == edg[j].p[1])
                    {
                        // The j-th edge corresponds to the edge we are looking for. Its index
                        // is saved in the triangle data, with the orientation.

                        e[k] = j;
                        o[k] = 1;
                        check++;

                    }
                    else if (p[edge_order[k]] == edg[j].p[1] && p[edge_order[k + 1]] == edg[j].p[0])
                    {
                        // We found a corresponding edge, but with opposite orientation.

                        e[k] = j;
                        o[k] = -1;
                        check++;

                    }

                }
            }

        }
        if(check < 3)
        {
            // Some edges are missing
            for (k = 0; k < 3; k++)
            {
                if (e[k] == -1)
                {
                    // If the edge is missing, a new entry in the edg[] array is created
                    // with the missing edge information.
                    edg[*nedges].p[0] = p[edge_order[k]];
                    edg[*nedges].p[1] = p[edge_order[k+1]];
                    e[k] = *nedges;
                    o[k] = 1;
                    (*nedges)++;
                }
            }
        }
    }

    // Resizes the edges array to the actual size.
    *edges = (edge*)realloc(edg, *nedges*sizeof(edge));
    if (*edges == NULL)
    {
        free(edg);
        return 1;
    }

    return 0;


}


bool is_in_triangle(point2d *p, point2d *p1, point2d *p2, point2d *p3)
{
    // Determines if a certain point is inside a triangle.
    // INPUT:
    // point2d* p: pointer to a 2D point.
    // point2d* p1, point2d* p2, point2d* p3: pointers to the vertices of the
    //      triangle.
    // OUTPUT:
    // The function returns true (-1) if the point p is inside the triangle defined
    // by *p1, *p2 and *p3, 0 otherwise.

    double d, q;

    // To determine if the point is inside the triangle, it computes the signed areas
    // of the triangles obtained substituting one of the vertices with p. If all the areas
    // and the signed area of the triangle have the same sign, then the point is inside the triangle.
    // However, due to rounding errors, if the point lies on an edge, the signed area may be a
    // small but non-zero number, of either sign. So we must allow a certain tolerance on the area
    // signs.
    d = (p2->x*p3->y + p3->x*p1->y + p1->x*p2->y - p2->x*p1->y - p3->x*p2->y - p1->x*p3->y);
    q = (p2->x*p3->y + p3->x*p->y + p->x*p2->y - p2->x*p->y - p3->x*p2->y - p->x*p3->y);

    if (d*q >= - NEARZERO_STRICT_TOLERANCE)
    {
        q = (p->x*p3->y + p3->x*p1->y + p1->x*p->y - p->x*p1->y - p3->x*p->y - p1->x*p3->y);

        if (d*q >= - NEARZERO_STRICT_TOLERANCE)
        {
            q = (p2->x*p->y + p->x*p1->y + p1->x*p2->y - p2->x*p1->y - p->x*p2->y - p1->x*p->y);

            if (d*q >= - NEARZERO_STRICT_TOLERANCE)
            {
                return -1;
            }

        }

    }

    return 0;
}


void interpolate_triangle(point2d *p, point2d *p1, point2d *p2, point2d *p3, doublecomplex *field, field2d *f)
{
    // Interpolates the field on a triangle using the edge-elements
    // basis functions.
    // INPUT:
    // point2d* p: a pointer to the 2D point in which the interpolation
    //      must be carried on.
    // point2d* p1, point2d* p2, point2d* p3: pointers to the vertices of
    //      the triangle.
    // doublecomplex* field: an array of size 3 containing the value of the field
    //      on each edge.
    // field2d* f: pointer to a field2d.
    // OUTPUT:
    // field2d* f: on exit, the vector value of the field at the interpolation point.

    double x = p->x;
    double y = p->y;

    double inv4sqarea = 1/(4*square(0.5*(p2->x*p3->y + p3->x*p1->y + p1->x*p2->y - p2->x*p1->y - p3->x*p2->y - p1->x*p3->y)));

    double ln1, ln2, ln3;

    double x1, x2, x3, y1, y2, y3;
    double N[6];

    // Edge lengths.
    ln1 = sqrt(square(p2->x - p1->x) + square(p2->y - p1->y));
    ln2 = sqrt(square(p3->x - p2->x) + square(p3->y - p2->y));
    ln3 = sqrt(square(p1->x - p3->x) + square(p1->y - p3->y));

    x1 = p1->x;
    x2 = p2->x;
    x3 = p3->x;
    y1 = p1->y;
    y2 = p2->y;
    y3 = p3->y;

    // Computes the x components of the 3 basis functions.
    N[0] = (ln1*(y - y3)*(x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3)))*(inv4sqarea);
    N[1] = (ln2*(y - y1)*(x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3)))*(inv4sqarea);
    N[2] = (ln3*(y - y2)*(x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3)))*(inv4sqarea);

    // Computes the y component of the 3 basis functions.
    N[3] = (ln1*(x - x3)*(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3)))*(inv4sqarea);
    N[4] = (ln2*(x - x1)*(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3)))*(inv4sqarea);
    N[5] = -(ln3*(x - x2)*(x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3)))*(inv4sqarea);


    // Computes the interpolated field.
    (f->x).re = N[0]*field[0].re + N[1]*field[1].re + N[2]*field[2].re;
    (f->x).i = N[0]*field[0].i + N[1]*field[1].i + N[2]*field[2].i;

    (f->y).re = N[3]*field[0].re + N[4]*field[1].re + N[5]*field[2].re;
    (f->y).i = N[3]*field[0].i + N[4]*field[1].i + N[5]*field[2].i;

}

int interpolate_2d(point2d *grid, int gridn, triangle *mtries, int ntries, point2d *mpoints, doublecomplex *mfield, field2d *efield)
{
    // Interpolates the electric field, given as edge-elements coefficients,
    // on a provided two-dimensional grid.
    // INPUT:
    // point2d* grid: a vector of size gridn containing the interpolation points.
    // int gridn: the number of points in the interpolation grid.
    // triangle *mtries: the array containing the triangle data for the mesh.
    // int ntries: the number of elements in the mtries[] array.
    // point2d *mpoints: the array containing the mesh nodes data.
    // doublecomplex *mfield: the field expressed as edge-element coefficients.
    // field2d *efield: pointer to an array of size gridn.
    // OUTPUT:
    // The function return 0 on success, and 1 if the grid is empty.
    // field2d *efield: On exit, the array of the interpolated field on the points
    //      contained in grid[].

    int i, j, k;
    bool found;
    point2d *p, *p1, *p2, *p3;
    doublecomplex *field;

    // Check if the grid is not empty.
    if (gridn <= 0)
        return 1;

    // Allocates an auxiliary vector to store the field values
    // for the current triangle
    field = malloc(3*sizeof(doublecomplex));

    for (i = 0; i < gridn; i++)
    {
        // Interpolates on each grid point.
        p = &grid[i];

        // Looks for the triangle containing the current point.
        found = 0;
        for (j = 0; j < ntries && !found; j++)
        {
            p1 = &mpoints[mtries[j].p[0]];
            p2 = &mpoints[mtries[j].p[1]];
            p3 = &mpoints[mtries[j].p[2]];

            if (is_in_triangle(p, p1, p2, p3))
            {
                // If the point is in the current triangle, saves the values of
                // the field on the edges (with the relative orientations)
                // in the auxiliary array field[].
                for (k = 0; k < 3; k++)
                {
                    field[k].re = mfield[mtries[j].e[k]].re*mtries[j].o[k];
                    field[k].i = mfield[mtries[j].e[k]].i*mtries[j].o[k];
                }
                // Interpolates the field in the triangle.
                interpolate_triangle(p, p1, p2, p3, field, &efield[i]);


                found = -1;
            }
        }

        if (!found)
        {
            // No triangle containing the point has been found. The field is set to zero arbitrarily.
            printf("WARNING: Interpolation point (%g, %g) outside of mesh.\n", p->x, p->y);
            efield[i].x.re = efield[i].x.i = efield[i].y.re = efield[i].y.i = 0;
        }

    }

    return 0;

}


