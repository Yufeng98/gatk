package org.broadinstitute.hellbender.utils.pairhmm;

/**
 * Superclass for PairHMM that want to use a full read x haplotype matrix for their match, insertion, and deletion matrix
 */
abstract class N2MemoryPairHMM extends PairHMM {
    protected float[][] transition = null; // The transition probabilities cache
    protected float[][] prior = null;      // The prior probabilities cache
    protected int[][] app_transition = null; // The transition probabilities cache
    protected int[][] app_prior = null;      // The prior probabilities cache
    protected float[][] matchMatrix = null;
    protected float[][] insertionMatrix = null;
    protected float[][] deletionMatrix = null;
    protected int[][] app_matchMatrix=null;
    protected int[][] app_insertionMatrix=null;
    protected int[][] app_deletionMatrix=null;
    protected int[][]   significance_i = null;
    protected int[][]   significance_j = null;
    //static final int MIN_INTEGER = Integer.MIN_VALUE/128;
    //static final int MIN_INTEGER = -1*(int)Math.pow(2,LoglessPairHMM.NUM_INT_BITS+LoglessPairHMM.NUM_FRACTION_BITS);//20'h80000
    static final int MIN_INTEGER = -524288;
    @Override
    public void doNotUseTristateCorrection() {
        doNotUseTristateCorrection = true;
    }

    /**
     * Initialize this PairHMM, making it suitable to run against a read and haplotype with given lengths
     *
     * Note: Do not worry about padding, just provide the true max length of the read and haplotype. The HMM will take care of the padding.
     *
     * @param haplotypeMaxLength the max length of haplotypes we want to use with this PairHMM
     * @param readMaxLength the max length of reads we want to use with this PairHMM
     */
    @Override
    public void initialize( final int readMaxLength, final int haplotypeMaxLength ) {
        super.initialize(readMaxLength, haplotypeMaxLength);

        matchMatrix = new float[paddedMaxReadLength][paddedMaxHaplotypeLength];
        insertionMatrix = new float[paddedMaxReadLength][paddedMaxHaplotypeLength];
        deletionMatrix = new float[paddedMaxReadLength][paddedMaxHaplotypeLength];
    	app_matchMatrix=new int[paddedMaxReadLength][paddedMaxHaplotypeLength];
    	app_insertionMatrix=new int[paddedMaxReadLength][paddedMaxHaplotypeLength];
    	app_deletionMatrix=new int[paddedMaxReadLength][paddedMaxHaplotypeLength];
	    //skip_match=new boolean[paddedMaxReadLength][paddedMaxHaplotypeLength];
	    //skip_delete=new boolean[paddedMaxReadLength][paddedMaxHaplotypeLength];
	    //skip_insert=new boolean[paddedMaxReadLength][paddedMaxHaplotypeLength];
        significance_i = new int[paddedMaxReadLength][paddedMaxHaplotypeLength];
        significance_j = new int[paddedMaxReadLength][paddedMaxHaplotypeLength];
        transition = PairHMMModel.createTransitionMatrix(maxReadLength);
        prior = new float[paddedMaxReadLength][paddedMaxHaplotypeLength];
        app_transition = new int[maxReadLength+1][6];
        app_prior = new int[paddedMaxReadLength][paddedMaxHaplotypeLength];
        for(int j=0;j<paddedMaxHaplotypeLength;j++){
            app_matchMatrix[0][j]=MIN_INTEGER;
            app_insertionMatrix[0][j]=MIN_INTEGER;
            //System.err.printf("%d\n",j);
            //System.err.println("Ileft:"+Integer.toHexString(app_matchMatrix[0][j]));
        }
        for(int i=0;i<paddedMaxReadLength;i++){
            app_matchMatrix[i][0] = MIN_INTEGER;
            app_insertionMatrix[i][0] = MIN_INTEGER;
            app_deletionMatrix[i][0] = MIN_INTEGER;
        }
    }

    /**
     * Print out the core hmm matrices for debugging
     */
    protected void dumpMatrices() {
        dumpMatrix("matchMetricArray", matchMatrix);
        dumpMatrix("insertionMatrix", insertionMatrix);
        dumpMatrix("deletionMatrix", deletionMatrix);
    }

    /**
     * Print out in a human readable form the matrix for debugging
     * @param name the name of this matrix
     * @param matrix the matrix of values
     */
    private void dumpMatrix(final String name, final float[][] matrix) {
        System.out.printf("%s%n", name);
        for ( int i = 0; i < matrix.length; i++) {
            System.out.printf("\t%s[%d]", name, i);
            for ( int j = 0; j < matrix[i].length; j++ ) {
                if ( Double.isInfinite(matrix[i][j]) )
                    System.out.printf(" %15s", String.format("%f", matrix[i][j]));
                else
                    System.out.printf(" % 15.5e", matrix[i][j]);
            }
            System.out.println();
        }
    }
}
