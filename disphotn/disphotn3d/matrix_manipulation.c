
// Macros to swap inline two integers or two double complex numbers.
#define SWAP(a,b,t) (t) = (a); (a) = (b); (b) = (t)
#define SWAP_COMPLEX(a,b,t) (t) = ((a).re); ((a).re) = ((b).re); ((b).re) = (t); (t) = ((a).i); ((a).i) = ((b).i); ((b).i) = (t)


int compare_matrix_entries(int *ia, int *ja, int i, int j)
{
    // Compares the position of two matrix entries in the same array.
    // INPUT:
    // int* ia, int* ja: the row and column indices arrays
    // int i, int j: the position of the two entries to compare within ia and ja.
    // OUTPUT:
    // Returns 1 if the entry at i follows the entry at j; -1 if precedes;
    // 0 if they are at the same position in the matrix.

    if(ia[i] < ia[j] || (ia[i] == ia[j] && ja[i] < ja[j]))
        return -1;
    else
        if(ia[i] == ia[j] && ja[i] == ja[j])
            return 0;
    return 1;
}

void swap_matrix_entries(int *ia, int *ja, doublecomplex *a, int i, int j)
{
    // Swaps two matrix entries in the same matrix.
    // INPUT:
    // int *ia, int *ja, doublecomplex *a: the arrays containing the information on the matrix
    // int i, int j: the indices (within the arrays) of the entries to swap.
    // OUTPUT:
    // the entries in the matrix arrays are swapped.

    int t;
    double d;

    SWAP(ia[i], ia[j], t);
    SWAP(ja[i], ja[j], t);
    SWAP_COMPLEX(a[i], a[j], d);

}

void sort_matrix_entries(int *ia, int *ja, doublecomplex *a, int first, int last)
{
    // Implementation of a quicksort algorithm that sorts the matrix entries
    // according to the rows index first and the column index second.
    // This function is recursive.
    // INPUT:
    // int *ia, int *ja, doublecomplex *a: the arrays containing the information on
    //      the matrix to sort.
    // int first, int last: the first and last VALID index within the arrays that
    //      define the section of the array to sort.

    int i, p;

    if(last - first == 0)
        return;

    if(last - first == 1)
    {
        if(compare_matrix_entries(ia, ja, last, first) < 0)
        {
            swap_matrix_entries(ia, ja, a, first, last);

        }
        return;

    }

    p = (int)(first + (last - first)/2);

    swap_matrix_entries(ia, ja, a, p, last);

    p = last-1;
    for(i = first; i <= p;)
    {
        if(compare_matrix_entries(ia, ja, i, last) > 0)
        {
            swap_matrix_entries(ia, ja, a, i, p);
            p--;
        }
        else
        {
            i++;
        }
    }

    if(compare_matrix_entries(ia, ja, p, last) < 0)
    {
        p++;
    }
    swap_matrix_entries(ia, ja, a, last, p);

    p--;
    if(p > first)
        sort_matrix_entries(ia, ja, a, first, p);

    p+=2;
    if(p < last)
        sort_matrix_entries(ia, ja, a, p, last);

}



int assembly_matrix(int **iap, int **jap, doublecomplex **ap, int n)
{
    // Sums all the entries at the same i,j position. The matrix entries
    // must be sorted first by row and then by column (i.e. this function
    // must be called after sort_matrix_entries). The function resizes the
    // arrays shortening them to the new size.
    // INPUT:
    // int **iap, int **jap, doublecomplex **ap: the pointers to the arrays
    //      containing the matrix information.
    // int n: the size of the arrays.
    // OUTPUT:
    // The function returns the new size of the arrays, or -1 if the reallocation
    //      of the memory failed.
    // The input pointers now point to new arrays containing the matrix information.

    int i = 0, j = 1;
    int *ia, *ja;
    doublecomplex *a;
    ia = *iap;
    ja = *jap;
    a = *ap;

    while(j < n)
    {
        if(ia[i] == ia[j] && ja[i] == ja[j])
        {
            a[i].re += a[j].re;
            a[i].i += a[j].i;
            j++;
        }
        else
        {
            i++;
            ia[i] = ia[j];
            ja[i] = ja[j];
            a[i].re = a[j].re;
            a[i].i = a[j].i;

            j++;
        }
    }
    i++;

    *iap = (int*)realloc(ia, i*sizeof(int));
    *jap = (int*)realloc(ja, i*sizeof(int));
    //printf("i = %d \n\n\n", i*sizeof(doublecomplex));
    *ap = (doublecomplex*)realloc(a, i*sizeof(doublecomplex));

    if(*iap == NULL || *jap == NULL || *ap == NULL)
    {
        free(ia);
        free(ja);
        free(a);
        return -1;


    }



    return i;

}


void set_PEC_conditions(int *ia, int *ja, doublecomplex *a, int nnz, bool *on_PEC, int nedges)
{
    // Imposes PEC conditions on a matrix (i.e. sets to zero the rows
    // and columns corresponding to edges on a PEC). The function does NOT
    // resize the matrix; the final matrix will therefore have empty rows and
    // colums. The task of removing them and resizing the matrix is carried on
    // by the purge_matrix function.
    // INPUT:
    // int *ia, int *ja, doublecomplex *a: the arrays containing the information on the matrix.
    // int nnz: the size of the arrays.
    // bool *on_PEC: an array of length nedges containing -1 if the edge is on a PEC boundary.
    // int nedges: the number of the edges in the problem (i.e. the size of boundary_info[]).
    // OUTPUT:
    // The values in the rows and columns of the matrix contained in the three arrays ia, ja and a
    // corresponding to edges on a non-periodic boundary are set to zero.

    int i, row, col;
    for(i = 0; i < nnz; i++)
    {
        // The matrix might be extended because of the presence of
        // metallic domains. Retrieve the edge number in this case
        // (the size of the matrix is always a multiple of nedges).
        row = ia[i];
        while(row - nedges >= 0)
            row -= nedges;
        col = ja[i];
        while(col - nedges >= 0)
            col -= nedges;

        if(on_PEC[row] || on_PEC[col])
        {
            // If the entry is in a "PEC" row or column, the value is set to 0.
            a[i].re = 0;
            a[i].i = 0;
        }

    }
}

void set_periodic_conditions(int *ia, int *ja, doublecomplex *a, int nnz, int nedges, int nbound,
                             int *source, int *destination, byte *orientations, doublecomplex phase)
{
    // Imposes Bloch Periodic Condition on a matrix.
    // The entries on the rows and columns corresponding to the source boundary edges
    // are moved on the rows and columns corresponding to the destination boundary edges,
    // and multiplied by a phase factor.
    // INPUT:
    // int *ia, int *ja, doublecomplex *a: the arrays containing the information on the matrix.
    // int nnz: the size of the arrays.
    // int nedges: the number of the edges of the problem.
    // int nbound: the number of the edges on the source and destination boundaries (should be equal).
    // int *source, int *destination: arrays containing the indices of the edges on the source and
    //      destination boundaries. The arrays must be ordered by correspondence, i.e. the i-th edge
    //      of the source array corresponds to the i-th element of the destination array.
    // byte *orientations: an array of length nbound containing the relative orientations of the
    //      corresponding edges (1 or -1).
    // doublecomplex phase: the (complex) Bloch phase factor.
    // OUTPUT:
    // on exit, the matrix information arrays contain the matrix with the boundary conditions imposed.
    // Note that now the matrix contains empty rows and columns corresponding to the indices in the
    // source array.

    int i, j, edgen;
    doublecomplex temp;

    for(i = 0; i < nnz; i++)
    {
        for(j = 0; j < nbound; j++)
        {
            // Retrieves the actual edge index even for expanded matrices.
            edgen = ia[i];
            while(edgen - nedges >= 0)
                edgen -= nedges;

            if(edgen == source[j])
            {
                // If the entry is on a row corresponding to a source edge, it is moved on the
                // corresponding destination row.

                ia[i] = destination[j];

                temp = a[i];

                a[i].re = orientations[j]*(phase.re*temp.re - phase.i*temp.i);
                a[i].i = orientations[j]*(phase.i*temp.re + phase.re*temp.i);

            }

            edgen = ja[i];
            while(edgen - nedges >= 0)
                edgen -= nedges;

            if(edgen == source[j])
            {
                // Same as above, but with columns.

                ja[i] = destination[j];

                temp = a[i];
                // For columns, the complex conjugate of the phase must be used
                a[i].re = orientations[j]*(phase.re*temp.re + phase.i*temp.i);
                a[i].i = orientations[j]*(-phase.i*temp.re + phase.re*temp.i);

            }
        }
    }
}


int remove_zeros(int **im, int **jm, doublecomplex **m, int *nnza)
{
    // Removes the entries with value zero from the matrix. The
    // matrix size is not changed, but the number of entries is.
    // INPUT:
    // int **im, int **jm, doublecomplex **m: pointers to the
    //      arrays containing matrix information.
    // int *nnza: pointer to an int containing the size of the arrays.
    // OUTPUT:
    // the function returns 0 on success, 1 if the resulting matrix is empty
    //      (i.e. the starting matrix was empty or all zeros), and 2 if the
    //      reallocation failed.
    // int **im, int **jm, doublecomplex **m: pointers to new arrays containing
    //      the information on the matrix with the zeros removed.
    // int *nnza: on exit, the address of an int containing the new number of
    //      nonzero elements (i.e. the length of the arrays).

    int i, j;
    int *ia, *ja, *ib, *jb;
    doublecomplex *a, *b;

    // Stores the addresses for convenience.
    ia = *im;
    ja = *jm;
    a = *m;

    // Some auxiliary arrays are used.
    ib = (int*)malloc(*nnza*sizeof(int));
    jb = (int*)malloc(*nnza*sizeof(int));
    b = (doublecomplex*)malloc(*nnza*sizeof(doublecomplex));

    // The counter j will be incremented with every valid (non-zero) entry that is preserved.
    j = 0;

    for(i = 0; i < *nnza; i++)
    {
        // Check if the value is zero within a certain tolerance
        // (to correct rounding errors).
        if(!near_zero_tol(a[i].re) || !near_zero_tol(a[i].i))
        {
            ib[j] = ia[i];
            jb[j] = ja[i];
            b[j] = a[i];
            j++;
        }
    }

    if(j == 0)
    {
        // The matrix is empty
        return 1;


    }

    *nnza = j;

    // Reallocates the arrays to the new size.
    free(ia);
    *im = (int*)realloc(ib, j*sizeof(int));
    if(*im == NULL)
    {
        free(ib);
        return 2;
    }

    free(ja);
    *jm = (int*)realloc(jb, j*sizeof(int));
    if(*jm == NULL)
    {
        free(jb);
        return 2;
    }

    free(a);
    *m = (doublecomplex*)realloc(b, j*sizeof(doublecomplex));
    if(*m == NULL)
    {
        free(b);
        return 2;
    }

    return 0;

}


int purge_matrix(int **im, int **jm, doublecomplex **m, int *nnza, int *nKeep, int **keepIndices )
{
    // This functions removes the rows and columns consisting of only zeros.
    // If the i-th row AND the i-th columns are empty, they are removed and the matrix shrunk.
    // However, only the emptyness of the rows is checked, since a correct imposition
    // of boundary conditions always produces pairs of empty rows and columns.
    // nKeep will contain the number of retained indices (the new dimension of the matrix),
    // while keepIndices will be an array of size nKeep containing the indices retained.
    // INPUT:
    // int **im, int **jm, doublecomplex **m: pointers to the
    //      arrays containing matrix information.
    // int *nnza: pointer to an int containing the size of the arrays.
    // OUTPUT:
    // the function returns 0 on success or 1 if the reallocation failed.
    //      If the function remove_zeros returns an error code, this function will return
    //      100 + the error code. Please see the code for remove_zeros for its error codes.
    // int **im, int **jm, doublecomplex **m: pointers to new arrays containing the information
    //      on the shrunk matrix.
    // int *nKeep: pointer to an int containin the new size of the matrix (i.e. the number of
    //      non-empty rows and columns found.
    // int **keepIndices: pointer to an array containing the indices of the non-empty rows and columns
    //      found.

    int i, n, j, k;
    int *ia, *ja, *keep;

    // Calls remove_zeros to remove the entries with value zero.
    i = remove_zeros(im, jm, m, nnza);

    if(i > 0)
    {
        // If remove_zeros returned an error, the function returns an error as well.
        return 100 + i;
    }

    ia = *im;
    ja = *jm;

    // Since the row array is ordered, its last element is the index of the last row,
    // i.e. the size of the matrix.
    n = ia[*nnza - 1];

    // Allocates space for the non-zeros rows and columns indices.
    keep = malloc(n*sizeof(int));
    *nKeep = 0;
    int currentRow = -1;

    // Scans all the matrix entries keeping track of the non-empty rows.
    // The indices of the non-empty rows are added to the keep array.
    for(i = 0; i < *nnza; i++)
    {
        if(ia[i] > currentRow)
        {
            currentRow = ia[i];
            keep[*nKeep] = ia[i];
            (*nKeep)++;
        }
    }

    // Resizes the array of the valid row indices to the actual size.
    *keepIndices = realloc(keep, *nKeep*sizeof(int));
    if(*keepIndices == NULL)
    {
        free(keep);
        return 1;
    }

    // Modifies the matrix entries updating the row and column indices
    // to reflect the fact that some rows and columns have been removed.

    // j is the index used to scan the keep[] array when comparing with the rows.
    j = 0;

    // k is the index used to scan the keep[] array when comparing with the columns.
    k = 0;

    // Scans all the matrix entries.
    for(i = 0; i < *nnza; i++)
    {

        if(ia[i] > keep[j])
        {
            // We moved to a new row, so the position in keep[] is increased...
            j++;
            // ...and the column-comparison counter is reset
            k = 0;
        }


        if(ia[i] == keep[j])
        {
            // The row is a valid one, so we will update the row index
            // to take into account the removed rows.
            ia[i] = j;

            // We scan the keep array and find the new index corresponding
            // to the column index. We don't need the "k=0" initialization
            // in the for cycle, since the column entries are sorted as well.
            for(;k < *nKeep; k++)
            {
                if(ja[i] == keep[k])
                {
                    // Update the column index with the new index that takes into
                    // account the removed columns.
                    ja[i] = k;
                    break;
                }
            }

        }
    }

    return 0;

}

void compress_rows(int *ia, int *ja, doublecomplex *a, int nnza)
{
    // Converts the matrix into Compressed Row Storage. On exit, the n-th
    // element in ia[] is the index in ja[] and a[] of the first element
    // of the n-th row.
    // INPUT:
    // int *ia, int *ja, doublecomplex *a: the arrays containing
    //      the information on the matrix
    // int nnza: the size of the arrays.
    // OUTPUT:
    // the array ia[] is changed to CRS. The other arrays are untouched.

    // currentRow is the
    int i, j = 1, currentRow = 0;
    int *in;

    // Allocates an auxiliary array
    in = malloc(nnza*sizeof(int));
    in[0] = 0;
    for (i = 1; i < nnza; i++)
    {
        if (ia[i] != currentRow)
        {
            // ia[] moved to a new row. Will add the current index (i)
            // to the in[] array until the current row number is moved
            // up to ia[].
            for (;currentRow < ia[i]; currentRow++)
            {
                in[j] = i;
                j++;

            }

            // This should be useless, but you may never know...
            currentRow = ia[i];
        }

    }

    // Copies the values from the auxiliary vector.
    for(i=0; i<j; i++)
    {
        ia[i] = in[i];
    }

    // Fills the ia vector with the last index.
    for (i = j; i < nnza;i++)
    {
        ia[i] = nnza;
    }

    // Frees the auxiliary vector.
    free(in);
}

void sparse_mul(int rows, int* ia, int* ja, doublecomplex* a, doublecomplex* y, doublecomplex* x)
{
    // Sparse Matrix-Vector multiplication. Performs y = A * x, with
    // A in Compressed Row Storage.
    // INPUT:
    // int rows: the size of the matrix (number of valid elements in the ia[] array)
    // int* ia, int* ja, doublecomplex* a: 1-based arrays containing the CRS information on the matrix.
    // doublecomplex* y: an array of size rows. It will be overwritten by this function.
    // doublecomplex* x: an array of size rows containing the x vector.
    // OUTPUT
    // doublecomplex* y: an array of size rows containing the result of the multiplication A * x

    // Note: since ia[] and ja[] are 1-based, all the values retrieved from these arrays mus be decremented by 1.

    int i, j, colpt;
    // Scans the rows
    for (i = 0; i < rows; i++)
    {

        // Computes the row-column product.
        y[i].re = 0;
        y[i].i = 0;
        // Scans the entries on the i-th row, retrieving the
        // column index from ja[].
        for (j = ia[i] - 1; j < ia[i+1] - 1; j++)
        {

            colpt = ja[j] - 1 ; // Position of the current entry
            y[i].re += (a[j].re * x[colpt].re - a[j].i*x[colpt].i);
            y[i].i += (a[j].i * x[colpt].re + a[j].re*x[colpt].i);
        }
    }
}


void sparse_aminussigmab(int n, doublecomplex* sigma, int *nnza, int** ia, int** ja, doublecomplex** a,
                         int nnzb, int* ib, int* jb, doublecomplex* b)
{
    // Computes A - sigma*B, where A and B are square matrices of the same size
    // stored in Compressed Row Storage and sigma is a complex number. The resulting
    // matrix is stored in A.
    // INPUT:
    // int n: the size of the matrices.
    // doublecomplex* sigma: pointer to the parameter sigma.
    // int nnza: the size of the arrays containing the A matrix information.
    // int** ia, int** ja, doublecomplex** a: pointers to the arrays containing
    //      information on the matrix A.
    // int nnzb: the size of the arrays containing the B matrix information.
    // int* ib, int* jb, doublecomplex* b: arrays containing
    //      information on the matrix B.
    // OUTPUT:
    // On exit, ia, ja, and a point to new arrays containing the information on
    // the A - sigma * B matrix.


    int i, cpa, cpb, cpam, cpbm, cpc;
    int *im, *jm, *ic, *jc, nnzc;
    doublecomplex *m, *c;

    // Computes the coefficients for the B matrix.
    double minsigmar = -(sigma->re), minsigmai = -(sigma->i);

    // Worst case prediction for the size of the new matrix.
    nnzc = *nnza+nnzb;

    // Retrieves the addresses for convenience.
    im = *ia;
    jm = *ja;
    m = *a;

    // Allocates memory in auxiliary arrays. Something like C = A - sigma * B is computed,
    // and then C will be copied into A.
    ic = malloc(sizeof(int)*nnzc);
    jc = malloc(sizeof(int)*nnzc);
    c =  malloc(sizeof(doublecomplex) * nnzc);

    // cpc is the current entry index for the C matrix.
    cpc = 0;


    // Scans the rows.
    for (i = 0; i < n; i++)
    {
        // Initializes the counters.
        // cpa and cpb are the indices of the beginning of the current row in ja[]
        // for A and B respectively. cpam and cpbm are the indices for the beginning
        // of the next row.
        cpa = im[i];
        cpb = ib[i];

        cpam = im[i+1];
        cpbm = ib[i+1];

        // This row in the C matrix begins at the current position.
        ic[i] = cpc;


        while (cpa < cpam || cpb < cpbm) // Cycle if there are still elements in the row.
        {

            if (cpa >= cpam)
            {
                // Row in A used (the cpa index has reached the end of the currend row,
                // saved in cpam). Will proceed copying the entries in B on this row and multiplying
                // them by -sigma.

                jc[cpc] = jb[cpb];
                c[cpc].re = minsigmar * b[cpb].re - minsigmai * b[cpb].i; // THE MINUS IS CORRECT!
                c[cpc].i = minsigmar * b[cpb].i + minsigmai * b[cpb].re;

                // It moved on one position in the B matrix, so the B counter is incremented.
                cpb++;
                // It also wrote a new entry in C, so the C counter must be incremented.
                cpc++;
            }
            else if (cpb >= cpbm)
            {
                // Row in B used. Simply copies the entries from A to C.
                jc[cpc] = jm[cpa];
                c[cpc] = m[cpa];

                // Increments the A counter,
                cpa++;
                // increments the C counter.
                cpc++;
            }
            else // There are still elements in both rows. Three situations can occur:
                if (jm[cpa] < jb[cpb])
            {
                // 1. There is an element in A without an element in B at the same position.
                // The value is copied.

                jc[cpc] = jm[cpa];
                c[cpc] = m[cpa];

                cpa++;
                cpc++;
            }
            else if (jb[cpb] < jm[cpa])
            {
                // 2. There is an element in B without an element in A at the same position.
                // The value is multiplied by -sigma and copied.

                jc[cpc] = jb[cpb];
                c[cpc].re = minsigmar * b[cpb].re - minsigmai * b[cpb].i; // THE MINUS IS CORRECT!
                c[cpc].i = minsigmar * b[cpb].i + minsigmai * b[cpb].re;

                cpb++;
                cpc++;
            }
            else
            {
                // 3. In this position, there is both an element in A and one in B.
                // Computes A(i,j)-sigma*B(i,j) and copies.
                jc[cpc] = jb[cpb];
                c[cpc].re = m[cpa].re + minsigmar * b[cpb].re - minsigmai * b[cpb].i; // THE MINUS IS CORRECT!
                c[cpc].i = m[cpa].i + minsigmar * b[cpb].i + minsigmai * b[cpb].re;

                cpa++;
                cpb++;
                cpc++;
            }
        }
    }

    // At the end of the cycle, cpc is the index of the last element (i.e. the size of the arrays)
    nnzc = cpc;

    // Fill the rest of the ic[] array with the last position.
    for(;i<nnzc;i++)
        ic[i] = cpc;

    // Frees the original A matrix.
    free(im);
    free(jm);
    free(m);

    // Resizes the resulting matrix arrays and put their addresses in the A pointers.
    *ia = realloc(ic, nnzc*sizeof(int));
    *ja = realloc(jc, nnzc*sizeof(int));
    *a = realloc(c, nnzc*sizeof(doublecomplex));

    // TO DO: pointers check

}
