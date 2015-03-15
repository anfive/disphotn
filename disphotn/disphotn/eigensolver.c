
// PARDISO routines. These are contained in the libpardiso* library.
extern  int pardisoinit_
    (void *, int *, int *, int *, double *, int *);

extern  int pardiso_
    (void *, int *, int *, int *, int *, int *,
     doublecomplex *, int *, int *, int *, int *, int *,
     int *, doublecomplex *, doublecomplex *, int *, double *);

// ARPACK routines.
extern void znaupd_(int *, char *, int *, char *, int *, double *,
                        doublecomplex *, int *, doublecomplex *, int *, int *, int *, doublecomplex *, doublecomplex *,
                        int *, double *, int *);

extern void zneupd_(int *, char *, int *, doublecomplex *, doublecomplex *,
                        int *, doublecomplex *, doublecomplex *, char *, int *, char *, int *,
                        double *, doublecomplex *, int *, doublecomplex *, int *, int *, int *,
                        doublecomplex *, doublecomplex *, int *, double *, int *);


int eigensolver(int **ia, int **ja, doublecomplex **a, int nnza, int **ib, int **jb, doublecomplex **b, int nnzb, int n,
                doublecomplex sigma, int nev, int ncv, doublecomplex *eigs, doublecomplex **eigv, int solver, int *iparam, int processors)
{
    // Eigensolving routine. Used to find nev eigenvalues for the problem
    // (A - sigma*B)^-1 *B*x = lambda*Bx
    // INPUT:
    // int *ia, int *ja, doublecomplex *a: the arrays containing the CRS data
    //      of the matrix A.
    // int nnza: the size of the arrays.
    // int *ib, int *jb, doublecomplex *b: doublecomplex *a: the arrays containing the CRS data
    //      of the matrix B.
    // int nnzb: the size of the arrays.
    // int n: the size of the problem (number of rows).
    // doublecomplex sigma: the eigenvalue shift.
    // int nev: the number of eigenvalues requestes.
    // int ncv: the size of Arnoldi factorization (used by ARPACK). If <= 0, it will be set to 8*nev
    //      by default.
    // doublecomplex* eigs: an array of size nev.
    // doublecomplex** eigv: a pointer to a pointer.
    // int solver: the solver used by PARDISO. 0 = iterative, 1 = direct.
    // int *iparam: an array of size 14 with the parameters for ARPACL.
    // OUTPUT:
    // The function returns 0 on success and a nonzero error code in case of error. It will print to
    // the console information on the eigensolving process, as well as any error that may occur.
    // doublecomplex *eigs: an array of size nev containing the computed eigenvalues.
    // doublecomplex **eigv: a pointer to an array of size nev*n containing the nev eigenvectors.


    // VARIABLE DECLARATION

    int i;
    doublecomplex *aux;

    // PARDISO
    long int pt[64]; // This should be int for 32bit systems
    int mtype, iparm[64], error;
    double dparm[64];
    int maxfct = 1, mnum = 1, phase = 13, nrhs = 1, *perm, msglvl = 0;

    // ARPACK
    char which[2], bmat;
    int ido, lworkl, info, ipntr[14], ldv;
    double tol;
    doublecomplex *workl, *workd, *resid, *v;
    double *rwork;

    // ARPACK Postprocessing
    int rvec, *select;
    int ldz = n;
    char howmny;
    doublecomplex *workenv;



    printf("    Generalized Eigenvalue problem with shift sigma = %g + %gi\n", sigma.re, sigma.i);


    if (sigma.re != 0 || sigma.i != 0)
    {
        // If a nonzero eigenvalue shift is required, computes A - sigma*B.
        // Note that due to the problem properties, an eigenvalue shift should be
        // always provided.
        printf("    Computing A - sigma*B...\n");
        sparse_aminussigmab(n, &sigma, &nnza, ia, ja, a, nnzb, *ib, *jb, *b);
        printf("    Done.\n");
    }

    // Since both PARDISO and ARPACK are written in FORTRAN,
    // we must update the index arrays to take into account the
    // fact that FORTRAN arrays are 1-based.
    for(i=0;i<nnza;i++)
    {
        (*ia)[i]++;
        (*ja)[i]++;
    }
    for(i=0;i<nnzb;i++)
    {
        (*ib)[i]++;
        (*jb)[i]++;
    }

    // Shortens the rows arrays (remove this?)
    *ia = realloc(*ia, (n+1)*sizeof(int));
    *ib = realloc(*ib, (n+1)*sizeof(int));


    // Initializes PARDISO parameters.

    error = 0;
    mtype = 13; // Complex nonsymmetric matrix
    maxfct = 1;
    mnum = 1;
    phase = 13;
    nrhs = 1;
    msglvl = 0;

    iparm[0] = 0;   // Use default parameters.
    iparm[2] = processors; // Number of processors.

    perm = malloc(n*sizeof(int));


    // Initializes ARPACK parameters.
    which[0] = 'L'; which[1] = 'M'; // "LM" must be used for the shift-and-invert method
    bmat = 'G'; // Generalized Eigenvalue problem


    if (ncv == 0)   // If the number of arnoldi vectors is not set, a standard size will be used.
        ncv = 8*nev;
    tol = 0;    // Tolerance

    ido = 0;
    info = 0; // Random starting vector

    // ARPACK Workspace
    lworkl = 3*(ncv*ncv) + 6*ncv;

    workl = malloc(lworkl*sizeof(doublecomplex));
    workd = malloc(3*n*sizeof(doublecomplex));
    resid = malloc(n*sizeof(doublecomplex));
    v = malloc(n*ncv*sizeof(doublecomplex));

    rwork = malloc(3*ncv*sizeof(double));

    ldv = n;

    iparam[0] = 1; // Shift strategy
    iparam[2] = ARPACK_MAX_ITERATIONS; // Max Iterations
    iparam[6] = 3; // Shift-invert mode


    // PARDISO Initialization
    printf("    Initializing PARDISO...\n");

    // Calls the PARDISO initialization routine (check license, etc.)
    pardisoinit_(pt, &mtype, &solver, iparm, dparm, &error);

    if (error)
    {
        // An error in pardisoinit has occurred.
        printf("ERROR: Pardisoinit Error %d\n", error);
        return 100;
    }
    printf("    Done.\n");


    printf("    Computing...\n");

    // Allocates auxiliary vector.
    aux = malloc(sizeof(doublecomplex) * (n + 1));

    // Main loop
    do
    {
        // Calls ARPACK update routine.
        znaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv,
                iparam, ipntr, workd, workl, &lworkl, rwork, &info );

        // According to ARPACK response (value of ido), decides what to do next.
        if (ido == -1)
        {
            // First iteration, or restart.
            // Perform inv(A - sigma*B)*B*x

            // aux = B*x...
            sparse_mul(n,*ib,*jb,*b,aux,&workd[ipntr[0]-1]);

            // ... and y = inv(A - sigma*B)*aux, using PARDISO.
            pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &n, *a, *ia, *ja,
                     perm, &nrhs, iparm, &msglvl, aux, &workd[ipntr[1]-1], &error, dparm);

            if (error != 0)
            {
                printf("ERROR: Pardiso Error %d\n", error);
                return 100;
            }
        }
        else if (ido == 1)
        {
            // Performs y = inv(A - sigma*B)*x, using PARDISO.
            pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &n, *a, *ia, *ja,
                     perm, &nrhs, iparm, &msglvl, &workd[ipntr[2]-1], &workd[ipntr[1]-1], &error, dparm);

            if (error != 0)
            {
                printf("ERROR: Pardiso Error %d\n", error);
                return 100;
            }



        }
        else if (ido == 2)
        {

            // Matrix-Vector Multiplication: y = B*x
            sparse_mul(n,*ib,*jb,*b,&workd[ipntr[1]-1],&workd[ipntr[0]-1]);
        }

    }
    while (ido == -1 || ido == 1 || ido == 2);

    // Frees the auxiliary vector.
    free(aux);

    if (info != 0)
    {
        // An ARPACK error has occurred.
        printf("ERROR: ARPACK Error %d\n", info);
        return 100;
    }

    // Releasing PARDISO memory
    phase = -1;
    pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &n, *a, *ia, *ja,
             perm, &nrhs, iparm, &msglvl, &workd[ipntr[2]-1], &workd[ipntr[1]-1], &error, dparm);

    printf("    Computation completed.\n");


    // ARPACK Postprocessing

    // Allocate space for the eigenvectors.
    *eigv = malloc((nev*n)*sizeof(doublecomplex));


    rvec = -1;
    howmny = 'A';
    select = malloc(ncv*sizeof(int));
    ldz = n;
    workenv = malloc(3*ncv*sizeof(doublecomplex));

    // Calls the ARPACK postprocessing routine to retrieve the eigenvectors.
    printf("    Postprocessing..\n");
    zneupd_(&rvec, &howmny, select, eigs, *eigv, &ldz, &sigma, workenv, &bmat, &n, which,
            &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
    printf("    Done.\n");

    if (info != 0)
    {
        // An ARPACK error has occurred.
        printf("ERROR: ARPACK Postprocessing Error %d\n", info);
        return 100;
    }

    // Frees workspace
    free(workl);
    free(workd);
    free(resid);
    free(v);
    free(rwork);

    free(*ia);
    free(*ja);
    free(*a);
    free(*ib);
    free(*jb);
    free(*b);

    *ia = *ja = *ib = *jb = NULL;
    *a = *b = NULL;

    // Returns.
    return 0;
}
