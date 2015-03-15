

void sparse_mul(int rows, int* ia, int* ja, doublecomplex* a, doublecomplex* y, doublecomplex* x)
{
    // Sparse Matrix-Vector multiplication

    int i, j, colpt;
    for (i = 0; i < rows; i++)
    {

        y[i].re = 0;
        y[i].i = 0;
        for (j = ia[i] - 1; j < ia[i+1] - 1; j++)
        {
            colpt = ja[j] - 1 ; // because of the zero-based C arrays
            y[i].re += (a[j].re * x[colpt].re - a[j].i*x[colpt].i);
            y[i].i += (a[j].i * x[colpt].re + a[j].re*x[colpt].i);
        }
    }
}


void sparse_aminussigmab(int n, doublecomplex* sigma, int *nnza, int** im, int** jm, doublecomplex** m,
                         int nnzb, int* ib, int* jb, doublecomplex* b)
{
    /*
    Sum of square sparse matrices in compressed row format.
    OUTPUT
    ic, jb, c must be arrays of length nnza+nnzb.
    */

    int i, cpa, cpb, cpam, cpbm, cpc=0;
    int *ia, *ja, *ic, *jc, nnzc;
    doublecomplex *a, *c;
    double minsigmar = -(sigma->re), minsigmai = -(sigma->i);

    nnzc = *nnza+nnzb;

    ia = *im;
    ja = *jm;
    a = *m;

    ic = (int *)malloc(sizeof(int)*nnzc);
    jc = (int *)malloc(sizeof(int)*nnzc);
    c = (doublecomplex *)malloc(sizeof(doublecomplex) * nnzc);



    for (i=0;i<n;i++)
    {
        cpa = ia[i];
        cpb = ib[i];

        cpam = ia[i+1];
        cpbm = ib[i+1];
        ic[i] = cpc;

        while (cpa < cpam || cpb < cpbm) // Cycle if there are still elements in the row
        {

            if (cpa >= cpam) // Row in a used, read from row in b.
            {
                jc[cpc] = jb[cpb];
                c[cpc].re = minsigmar * b[cpb].re - minsigmai * b[cpb].i; // THE MINUS IS CORRECT!
                c[cpc].i = minsigmar * b[cpb].i + minsigmai * b[cpb].re;

                cpb++;
                cpc++;
            }
            else if (cpb >= cpbm) // Row in b used, read from row in a.
            {
                jc[cpc] = ja[cpa];
                c[cpc] = a[cpa];

                cpa++;
                cpc++;
            }
            else if (ja[cpa] < jb[cpb]) // Both rows still active
            {
                jc[cpc] = ja[cpa];
                c[cpc] = a[cpa];

                cpa++;
                cpc++;
            }
            else if (jb[cpb] < ja[cpa])
            {
                jc[cpc] = jb[cpb];
                c[cpc].re = minsigmar * b[cpb].re - minsigmai * b[cpb].i; // THE MINUS IS CORRECT!
                c[cpc].i = minsigmar * b[cpb].i + minsigmai * b[cpb].re;

                cpb++;
                cpc++;
            }
            else // Same position, sum
            {
                jc[cpc] = jb[cpb];
                c[cpc].re = a[cpa].re + minsigmar * b[cpb].re - minsigmai * b[cpb].i; // THE MINUS IS CORRECT!
                c[cpc].i = a[cpa].i + minsigmar * b[cpb].i + minsigmai * b[cpb].re;

                cpa++;
                cpb++;
                cpc++;
            }
        }
    }

    nnzc = cpc;

    for(;i<nnzc;i++)
        ic[i] = cpc;

    free(ia);
    free(ja);
    free(a);

    *im = realloc(ic, nnzc*sizeof(int));
    *jm = realloc(jc, nnzc*sizeof(int));
    *m = realloc(c, nnzc*sizeof(doublecomplex));

    // TO DO: pointers check

}

int sparse_orthogonalize(int *n, int *id, int *jd, doublecomplex *d, doublecomplex *v)
{
    return 0;


}
