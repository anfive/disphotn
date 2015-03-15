
int scan_line(FILE* fid, char* line)
{
    // Reads the next valid line from a file.
    // The function keeps reading lines from the stream pointed by fid
    // removing comments (identified by the comment character #) until
    // a non-empty string is obtained or the stream is ended.
    // Input: FILE *fid: a valid pointer to a file stream;
    // char* line: a pointer to a char array of minimum 256 elements.
    // Output: the function returns -1 if a valid line has been read.
    // char* line: the first valid line read.

    int i;

    if (fid == NULL)
    {

        return 0;
    }

    if (feof(fid))
    {
        return 0;
    }


    do
    {
        // Reads a line from the stream
        fgets(line, 256, fid);

        // Removes comments and newline characters
        for (i = 0; line[i] != '\0'; i++)
        {
            if (line[i] == '#' || line[i] == '\n' || line[i] == '\r')
            {
                // The delimiter character is replaced with an end-of-string character.
                // The program will ignore the subsequent characters in all the other functions.
                line[i] = '\0';
                break;
            }
        }

    }
    while (line[0] == '\0' && feof(fid) == 0);

    if (feof(fid))
        return 0;
    else
        return -1;

}

int read_configuration(char *filename, strlist **var, strlist **val)
{
    // Reads the configuration file. The generic format for a configuration
    // file line is:
    // parameter-name [value1 [value2 [...]]]
    // The function creates two lists of strings, one containing the parameter names
    // (from the beginning of the line to the first blank), and the other containing the
    // rest of the line (to the newline or the comment character '#').
    // INPUT:
    // char *filename: a string containing the configuration file name.
    // strlist** var, strlist** val: two pointers to pointers.
    // OUTPUT:
    // The function returns -1 on success and 0 if the file could not be open.
    // strlist** var, strlist** val: on exit, the pointers point to the address of the
    //      beginning of each list. The var list contains the parameter names, the
    //      val list contains the parameter values (in string form).

    strlist *pvar, *pval;
    char *tmp1, *tmp2;

    // Opens file
    FILE *fid = fopen(filename, "r");
    if (fid == NULL)
    {
        return 0;
    }
    *var = *val = NULL;

    // Allocates auxiliary strings
    tmp1 = malloc(sizeof(char)*256);
    tmp2 = malloc(sizeof(char)*256);

    // Scans the first line
    scan_line(fid, tmp1);

    // Splits the line in two tokens at the first blank space.
    tmp1 = strtok(tmp1, " ;#\r\n");
    tmp2 = strtok(NULL, ";#\r\n");

    while (feof(fid) == 0)
    {
        if (*var != NULL)
        {
            // Adds an element to the list
            pvar = pvar->n = malloc(sizeof(strlist));
            pval = pval->n = malloc(sizeof(strlist));
        }
        else
        {
            // First list element
            pvar = *var = malloc(sizeof(strlist));
            pval = *val = malloc(sizeof(strlist));
        }

        // Saves the current values in the list element.
        pvar->v = tmp1;
        pval->v = tmp2;
        pvar->n = NULL;
        pval->n = NULL;

        // Allocates memory for the next list values
        tmp1 = malloc(sizeof(char)*256);
        tmp2 = malloc(sizeof(char)*256);

        // Reads the next line.
        scan_line(fid, tmp1);

        // Splits the line in two tokens at the first blank space.
        tmp1 = strtok(tmp1, " ;#\r\n");
        tmp2 = strtok(NULL, ";#\r\n");



    }

    fclose(fid);
    return -1;

}

int read_grid_file_2d(char *filename, point2d **grid, int *gridn)
{
    // Reads a file containing the interpolation points.
    // The first line of the file must be an integer equal to the
    // number of grid points.
    // The following lines must be couples of double values that are
    // interpreted as the (x,y) coordinates of the interpolation points.
    // INPUT:
    // char* filename: a string containing the grid file name.
    // point2d** grid: a pointer to a point2d pointer.
    // int* gridn: the address of an int.
    // OUTPUT:
    // The function returns -1 on success, or 0 if the file could not be open.
    // point2d** grid: on exit, it points to the array of grid points read.
    // int* gridn: on exit, it point to the number of grid points read.

    FILE *fid;
    char line[256];
    int i;

    fid = fopen(filename, "r");
    if(fid == NULL)
        return 0;

    // Reads the number of points.
    *gridn = atoi(line);

    // Allocates the memory.
    *grid = malloc(*gridn*sizeof(point2d));
    for(i = 0; i < *gridn; i++)
    {
        // Reads the following lines an retrieves the coordinates of the points.
        if(!scan_line(fid, line))
            return 0;
        sscanf(line, "%lg %lg", &((*grid)[i].x), &((*grid)[i].y));
    }

    return -1;

}

// ********* RESULTS WRITING FUNCTIONS **********
// The following 3 functions are used to write the results (eigenfrequencies,
// e-field and h-field) to file. They are quite self-explanatory.

int write_bands(char *filename, doublecomplex *bands, int nev, int kvals)
{
    FILE *fid;
    int k, i;

    fid = fopen(filename, "w");
    if(fid == NULL)
    {
        return 1;
    }

    for(i = 0; i < nev; i++)
    {
        for(k = 0; k < kvals; k++)
        {
            if(near_zero_tol(bands[k*nev + i].i))
                fprintf(fid, "%+.50g ", bands[k*nev + i].re);
            else
                fprintf(fid, "%+.50g%+.50gi ", bands[k*nev + i].re, bands[k*nev + i].i);
        }
        fprintf(fid, "\r\n");
    }

    fclose(fid);
    return 0;
}

int write_modes(char *filename, field2d *modes, point2d *grid, int gridn, int nev, int kvals)
{
    FILE *fid;
    int k, i,j;

    fid = fopen(filename, "w");
    if(fid == NULL)
    {
        return 1;
    }

    for(j = 0; j < gridn; j++)
    {
        fprintf(fid, "%g ", grid[j].x);
    }
    fprintf(fid, "\r\n");
    for(j = 0; j < gridn; j++)
    {
        fprintf(fid, "%g ", grid[j].y);
    }
    fprintf(fid, "\r\n");

    // x component
    for(k = 0; k < kvals; k++)
        for(i = 0; i < nev; i++)
        {
            for(j = 0; j < gridn; j++)
            {
                fprintf(fid, "%+g%+gi ", modes[(k*nev + i)*gridn + j].x.re, modes[(k*nev + i)*gridn + j].x.i);
            }
            fprintf(fid, "\r\n");
        }

    // y component
    for(k = 0; k < kvals; k++)
        for(i = 0; i < nev; i++)
        {
            for(j = 0; j < gridn; j++)
            {
                fprintf(fid, "%+g%+gi ", modes[(k*nev + i)*gridn + j].y.re, modes[(k*nev + i)*gridn + j].y.i);
            }
            fprintf(fid, "\r\n");
        }

    fclose(fid);
    return 0;
}

int write_hfield(char *filename, doublecomplex *hfield, int ntries, int nev, int kvals)
{
    FILE *fid;
    int k, i,j;

    fid = fopen(filename, "w");
    if(fid == NULL)
    {
        return 1;
    }

    for(k = 0; k < kvals; k++)
        for(i = 0; i < nev; i++)
        {
            for(j = 0; j < ntries; j++)
            {
                fprintf(fid, "%+g%+gi ", hfield[(k*nev + i)*ntries + j].re, hfield[(k*nev + i)*ntries + j].i);
            }
            fprintf(fid, "\r\n");
        }

    fclose(fid);
    return 0;
}

// ******* END OF RESULTS WRITING FUNCTIONS *********

void read_sequence(char *line, void *v, int *nvals, char const *format)
{
    // Reads a series of values from a string.
    // INPUT:
    // char* line: the string containing the values.
    // OUTPUT:
    // double* kv: the values read from the string.
    // int* kvals: the number of values read.

    char *t;
    *nvals = 0;
    t = strtok(line, " ;#\r\n");
    while(t != NULL && *nvals < MAX_SEQ_LEN)
    {
        sscanf(t, format, v + (*nvals)*sizeof(v));
        (*nvals)++;
        t = strtok(NULL, " ;#\r\n");
    }
}

int read_mesh_file_2d(char* filename, point2d** points, int* npoints, edge** edges, int* nedges, triangle** tries, int* ntries)
{
    // Reads a COMSOL mesh file containing a 2-dimensional mesh composed by triangles.
    // Input:
    // char* filename: the path and name of the mesh file
    // Output:
    // The function returns 0 if the reading succeeded. The function returns an error code > 0 in case of error.
    // point2d** points: a pointer to an array containing the information on the nodes.
    // int* npoints: a pointer to an integer containing the number of nodes.
    // edge** edges: a pointer to an array containing the information on the edges (2-node elements).
    // int* nedges: a pointer to an integer containing the number of edges.
    // triangle** tries: a pointer to an array containing the information on the triangles (3-node elements).
    // int* ntries: a pointer to an integer containing the number of triangle.
    // ERROR CODES:
    // 0: No error, reading succeeded.
    // 1: Error in opening file.
    // 2: Error in reading an element: not enough nodes.
    // 3: Inconsistency on the number of elements and the number of domains for a certain element type.

    FILE* fid;

    point2d *pnts;
    edge *edgs;
    triangle *tris;
    void *current;
    char *tp;


    char line[256], a[256];
    int dummy, i, j, k;
    int num_el_types, num_nodes, num_el, size_el;

    // Opens file
    fid = fopen(filename, "r");
    if (fid == NULL)
    {
        return 1;
    }

    // Reads the file header, containing information unused by the program.
    scan_line(fid, line); // Version
    scan_line(fid, line); // Number of tags
    dummy = atoi(line);
    for (i = 0; i < dummy; i++)
    {
        scan_line(fid, line); // Tags
    }

    scan_line(fid, line); // Number of types
    dummy = atoi(line);
    for (i = 0; i < dummy; i++)
    {
        scan_line(fid, line); // Types
    }

    scan_line(fid, line); // ???
    scan_line(fid, line); // Class
    scan_line(fid, line); // Version
    scan_line(fid, line); // SDIM

    // *********** POINTS ***********
    // Reads the information on points (nodes)
    scan_line(fid, line); // Number of points
    *npoints = atoi(line);
    scan_line(fid, line); // Lowest mesh point index


    pnts = (point2d*)malloc(*npoints * sizeof(point2d));
    for (i = 0; i < *npoints; i++)
    {
        scan_line(fid, line);
        sscanf(line, "%lf %lf", &pnts[i].x, &pnts[i].y);

        //printf("x = %lf, y = %lf\n", pnts[i].x, pnts[i].y);
    }


    // *********** ELEMENTS **********
    // Reads the information on elements (groups of n nodes)
    scan_line(fid, line); // Number of element types
    num_el_types = atoi(line);

    //printf("Reading %d element types\n", num_el_types);

    for (i = 0; i < num_el_types; i++)
    {
        scan_line(fid, line); // Type name
        scan_line(fid, line); // Number of Nodes per Element
        num_nodes = atoi(line);
        scan_line(fid,line); // Number of Elements
        num_el = atoi(line);
        size_el = 0;

        //printf("Element type %d, %d nodes per element, %d elements\n", i+1, num_nodes, num_el);

        // Identifies the element type by the number of nodes. The COMSOL element type name is not used.
        // The reading process is the same for each element type. Here are specified the arrays in which
        // store the information read from the stream. The "current" void pointer points to the destination array,
        // and the "size_el" int specifies the size of each element in the array.
        // This method allows the reading of arrays of different types (edge, triangle, etc.) using the
        // same shared code.
        // if current == NULL, the information is read from the stream but is not saved anywhere. This ensures
        // robustness against unrecognized element types.
        switch (num_nodes)
        {
        case 1:
        {
            // Vertex (point)
            current = NULL; // unused
        }
        break;
        case 2:
        {
            // Edge
            edgs = (edge*)malloc(num_el*sizeof(edge));
            current = (void*)edgs;
            *nedges = num_el;
            size_el = sizeof(edge);
        }
        break;
        case 3:
        {
            // Triangle
            tris = (triangle*)malloc(num_el*sizeof(triangle));
            current = (void*)tris;
            *ntries = num_el;
            size_el = sizeof(triangle);
        } break;
        default:
        {
            // Unrecognized
            current = NULL;
        }
        break;


        }

        //printf("reading nodes\n");


        // Reads the element information
        for (j = 0; j < num_el; j++)
        {
            scan_line(fid, line);
            if (current != NULL) // If the current n-nodes element must be stored,...
            {
                strcpy(a, line);

                // ...the line is split into n tokens, each containing a node number...
                tp = strtok(a, " ");
                for (k = 0; k < num_nodes; k++)
                {
                    if(tp == NULL)
                    {
                        return 2;
                    }


                    //sscanf(a, "%s", b);
                    // ...that is saved into the desired array. The following line is
                    // a bit cryptic, but the only thing it does is to save the k-th node
                    // number in the current (j-th) element of the array.
                    // To work properly, all the element data types (structs) must have
                    // the n node numbers (ints) as their first fields.
                    *((int*)(current + j*size_el + k*sizeof(int))) = atoi(tp);

                    tp = strtok(NULL, " ");
                    //strcpy(a, c);

                }
            }
        }

        //printf("Reading parameters for elements %d\n", i+1);

        // The information on parameter values is unused by the program.
        scan_line(fid, line); // Number of parameter values per element
        scan_line(fid, line); // Number of parameters
        dummy = atoi(line);
        for (j = 0; j < dummy; j++)
        {
            scan_line(fid, line); // Parameters
        }

        // ********** DOMAINS ************
        // Reads information on the domains. Domains are used to define freespace/dielectric/metallic
        // regions in the simulation.
        scan_line(fid, line); // Number of domains
        dummy = atoi(line);
        //printf("Reading %d domains for elements %d\n",dummy, i+1);

        // Checks that the number of domains specified is equal to the number of elements.
        if(dummy != num_el)
        {
            return 3;
        }

        // The following piece of code works similarly to the one that reads the nodes.
        // To work properly, the element data struct for the n-noded element type must have the
        // domain number (an int) as the (n+1)th field.
        for (j = 0; j < dummy; j++)
        {
            scan_line(fid, line); // Domains
            if(current != NULL)
            {
                *((int*)(current + j*size_el + num_nodes*sizeof(int)))  = atoi(line);
            }
        }

        // The information on up/down pairs is unused by the program.
        scan_line(fid, line); // Number of up/down pairs
        dummy = atoi(line);
        for (j = 0; j < dummy; j++)
        {
            scan_line(fid, line); // Pairs
        }
    }


    // Copies the pointers to the newly created arrays to the output variables.
    *points = pnts;
    *edges = edgs;
    *tries = tris;

    // Closes the file.
    fclose(fid);
    return 0;
}
