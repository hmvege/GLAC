int index(int i, int j, int k, int l, int N)
{
    /*
     * Function for contigious memory allocation.
     */
    return (N*(N*(N*i + j) + k) + l);
}

int index2(int i, int j, int N)
{
    /*
     * Function for contigious memory allocation for the SU(3) matrices.
     */
    return N*i + j;
}
