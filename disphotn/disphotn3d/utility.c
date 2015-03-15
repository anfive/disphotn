// Utility Functions. Used troughout the program.

double det3(point3d v1, point3d v2, point3d v3)
{
    // Determinant of a 3x3 matrix obtained by juxtaposition of three 3-elements vectors.
    return v1.x*v2.y*v3.z - v1.x*v2.z*v3.y - v1.y*v2.x*v3.z + v1.y*v2.z*v3.x + v1.z*v2.x*v3.y - v1.z*v2.y*v3.x;
}

double dot_product_3d(point3d *v1, point3d *v2)
{
    // Vector dot product.
    return (v1->x*v2->x + v1->y*v2->y + v1->z*v2->z);

}

/*
int valid_number(double a)
{
    return (a == a && (a + 1) != a && (a - 1) != a);


}
*/

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

point3d vector_subtract_3d(point3d v1, point3d v2)
{
    // Subtracts two 2-dimensional vectors (v1 - v2)
    point3d p;
    p.x = v1.x - v2.x;
    p.y = v1.y - v2.y;
    p.z = v1.z - v2.z;
    return p;
}

double square(double a)
{
    // Square of a double-precision number
    return a*a;
}

doublecomplex compute_bloch_phase_3d(double kx, double ky, double kz, point3d *a)
{
    // Computes the Bloch phase coefficient for given values of the vectors k and a.
    // The formula is phase = exp(i*(k.a))
    doublecomplex phase;
    phase.re = cos(kx*a->x + ky*a->y + kz*a->z);
    phase.i = sin(kx*a->x + ky*a->y + kz*a->z);
    return phase;
}


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
