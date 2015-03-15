// Utility Functions. Used troughout the program.

double det2(point2d v1, point2d v2)
{
    // Determinant of a 2x2 matrix obtained by juxtaposition of two 2-elements vectors.
    return (v1.x*v2.y)-(v1.y*v2.x);
}


int near_zero(double a)
{
    // Checks if a certain value is smaller than a certain tolerance. Used in geometrical
    // comparisons. The tolerance is specified in the NEARZERO_TOLERANCE constant.
    return (a < NEARZERO_TOLERANCE && a > - NEARZERO_TOLERANCE);
}

int near_zero_tol(double a)
{
    // Checks if a certain value is smaller than a certain tolerance. Used to correct
    // rounding errors. The tolerance is specified in the NEARZERO_STRICT_TOLERANCE constant.
    return (a < NEARZERO_STRICT_TOLERANCE && a > - NEARZERO_STRICT_TOLERANCE);

}

point2d vector_subtract(point2d v1, point2d v2)
{
    // Subtracts two 2-dimensional vectors (v1 - v2)
    point2d p;
    p.x = v1.x - v2.x;
    p.y = v1.y - v2.y;
    return p;
}

double square(double a)
{
    // Square of a double-precision number
    return a*a;
}

/*
void read_k_values(char *line, double *kv, int *kvals)
{
    // Reads a series of double values from a string.
    // INPUT:
    // char* line: the string containing the values.
    // OUTPUT:
    // double* kv: the values read from the string.
    // int* kvals: the number of values read.

    char *t;
    *kvals = 0;
    t = strtok(line, " ;#\r\n");
    while(t != NULL && *kvals < MAX_SEQ_LEN)
    {
        sscanf(t, "%lg", &(kv[*kvals]));
        (*kvals)++;
        t = strtok(NULL, " ;#\r\n");
    }
}
*/

/*
int check_sorting(int *ia, int *ja, int nnza)
{
    int i;
        //sort check
    for(i = 1;i < nnza; i++)
    {
        if(ia[i] < ia[i-1] || (ia[i] == ia[i-1] && ja[i] < ja[i-1]))
        {
            return 0;
        }
    }
    return -1;



}
*/
