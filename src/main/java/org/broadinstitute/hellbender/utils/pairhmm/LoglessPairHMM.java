package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.utils.QualityUtils;

import static org.broadinstitute.hellbender.utils.pairhmm.PairHMMModel.*;
import java.io.*;
import java.util.Random;
import java.util.Arrays;

public class LoglessPairHMM extends N2MemoryPairHMM {
	//static float CURRENT_ERROR=(float)0.2;
	//static float LOG2CURRENT_ERROR=(float)(Math.log(CURRENT_ERROR)/Math.log(2));
    static int LOG2CURRENT_ERROR = Float2Fix((float)(Math.log(0.2)/Math.log(2)));
    static final float INITIAL_CONDITION = (float)Math.pow(2, 127);//originally 1020
    static final float INITIAL_CONDITION_LOG10 = (float)Math.log10(INITIAL_CONDITION);
    static int NUM_FRACTION_BITS=14;////was 12
    static int NUM_INT_BITS=5; //was 10
    static int MAX_INTEGER=Find_Max_Integer(NUM_INT_BITS)<<NUM_FRACTION_BITS;
    static final float INITIAL_CONDITION_UP = (float)Math.pow(2, Find_Max_Integer(NUM_INT_BITS));//originally 1020
    static int MAX_RANGE=NUM_FRACTION_BITS;//was NUM_FRACTION_BITS
    static int MAX_RANGE_BIN=8<<NUM_FRACTION_BITS;//was MAX_RANGE
    static int[] LOGSUM_TABLE_ENTRY=createLogSumTableEntry(MAX_RANGE);
    static int [] LOGSUM_VALUE_ENTRY=createLogSumValueEntry();
    static int MIN_FRACTION = LOGSUM_VALUE_ENTRY[LOGSUM_VALUE_ENTRY.length-1];
	static int[] sum_bin_match= new int[2*(int)Math.ceil((float)(MAX_INTEGER+1)/MAX_RANGE_BIN)];
	static int[] sum_bin_insertion= new int[2*(int)Math.ceil((float)(MAX_INTEGER+1)/MAX_RANGE_BIN)];
        //use this only to match hardware
	static int[] sum_bin_deletion= new int[2*(int)Math.ceil((float)(MAX_INTEGER+1)/MAX_RANGE_BIN)];
        //use this only to match hardware
    private static int[] createLogSumTableEntry(int range){
        int[] result = new int[range+1];
        for(int i=0;i<result.length;i++){
            result[i]=-1*i;
        }
        return result;
    }
    private static int[] createLogSumValueEntry(){
        int[] result = new int[LOGSUM_TABLE_ENTRY.length];
        //System.err.println("logsumvalue");
        for(int i=0;i<LOGSUM_TABLE_ENTRY.length;i++){
            result[i]=Float2Fix((float)(Math.log(1+Math.pow(2,LOGSUM_TABLE_ENTRY[i]))/Math.log(2)));
            //System.err.printf("%d %f\n",i,Fix2Float(result[i]));
            //System.err.printf("i:%d "+Integer.toHexString(result[i])+"\n",i);
        }
        return result;
    }
    static final float NOCOMPUTETH = -10;

    static final int[] INITIALVALUE_TABLE=createInitialValueTable();
    static final float[] OFFSET_TABLE = createOffsetTable();
    static final int TRISTATE_CORRECTION = 3;
    static final byte INSQ = 45;
    static final byte DELQ = 45;
    static final byte GCPQ = 10;
    //long upperlogsumTime = 0;
    static final int[] TRANSITION_TABLE = {
        Upper_LOG2_accurate(PairHMMModel.matchToMatchProb(INSQ,DELQ)),//M2M
        Upper_LOG2_accurate((float)QualityUtils.qualToProb(GCPQ)),//indel2M
        Upper_LOG2_accurate((float)QualityUtils.qualToErrorProb(INSQ)),//M2I
        Upper_LOG2_accurate((float)QualityUtils.qualToErrorProb(GCPQ)),//I2I
        Upper_LOG2_accurate((float)QualityUtils.qualToErrorProb(DELQ)),//M2D
        Upper_LOG2_accurate((float)QualityUtils.qualToErrorProb(GCPQ))//D2D
    };

    static final int[] QUAL2PROB_TABLE= createQual2ProbTable();
    static final int[] QUAL2ERROR_DIV3_TABLE = createQual2ErrorDiv3Table();
    
    private static int[] createQual2ProbTable(){
        int[] result = new int[256];//256 because quality score is a byte
        for(int i=0;i<64;i++){
            result[i] = Upper_LOG2_accurate((float)QualityUtils.qualToProb(i));
        }
        for(int i=64;i<result.length;i++){
            result[i]=result[63];
        }
        //print for debug
        //System.err.println("q2p");
        //for(int i=0; i<result.length;i++){
        //    System.err.println(Integer.toBinaryString(result[i]));
        //}   
        //print for debug
        return result; 
    }
    
    private static int[] createQual2ErrorDiv3Table(){
        int[] result = new int[256];
        for(int i=0;i<64;i++){
            result[i] = Upper_LOG2_accurate((float)QualityUtils.qualToErrorProb(i) / TRISTATE_CORRECTION);
        }
        for(int i=64;i<result.length;i++){
            result[i]=result[63];       
        }
        //print for debug
        //System.err.println("q2e/3");
        //for(int i=0; i<result.length;i++){
        //    System.err.println(Integer.toBinaryString(result[i]));
        //}   
        //print for debug
        return result; 
    }
    
    private static int[] createInitialValueTable(){
        int[] initialTable = new int[700];
        initialTable[0]=0;//meaningless
        for(int i=1;i<initialTable.length;i++){
            initialTable[i] = Upper_LOG2_accurate(INITIAL_CONDITION_UP / i);
            //System.err.printf("i=%d log2=%f result=%f\n",i,Math.log(INITIAL_CONDITION_UP / i)/Math.log(2),Fix2Float(initialTable[i]));
        }
        return initialTable;
    }
    private static float[] createOffsetTable(){
        float [] offsetTable = new float[700];
        offsetTable[0]=0;//meaningless
        for(int i=1; i<offsetTable.length;i++){
            offsetTable[i] = (float)Math.log10(Math.pow(2,Fix2Float(INITIALVALUE_TABLE[i]))* i);
        }
        return offsetTable;
    }
    private static float Fix2Float(int integer){
        float result = (float)(integer/Math.pow(2,NUM_FRACTION_BITS));
        //System.err.printf("fix=%d fix>>=%d float=%f\n",integer,integer>>NUM_FRACTION_BITS,result);
        return result;
    }
    
    // we divide e by 3 because the observed base could have come from any of the non-observed alleles
    public static int Find_Max_Integer(int bits){
	return (int)Math.pow(2,bits)-1; 
    }
    
    //Convert a 32bit floating point number to a fix point number represented with integer
    public static int Float2Fix(float exact_float){
	    if(exact_float==Float.NEGATIVE_INFINITY){
	    	return MIN_INTEGER;
	    }
        int result = (int)Math.ceil(exact_float*Math.pow(2,NUM_FRACTION_BITS));
        //System.err.printf("exact_float=%f result=%f\n",exact_float,Fix2Float(result)); 
        return result;
    }

    //More accurate log2
    private static int Upper_LOG2_accurate(float num){
        //long start=System.nanoTime();
        float numLog2 = (float)(Math.log(num)/Math.log(2));
        int result = Float2Fix(numLog2);
        //System.err.printf("exact_fix=%f fix=%f\n",numLog2,Fix2Float(result));
        //long elapsedTime = System.nanoTime()-start;
        //System.err.printf("time to do upper_log2_accurate=%d\n",elapsedTime);
        return result;
    }

    //fixed point implementation with exact diff (different from current hardware)
    private static int Upper_LOGSUM(int num1,int num2){
        int diff = Math.min(MAX_RANGE,(Math.abs(num1 - num2))>>NUM_FRACTION_BITS);
        //if(diff>MAX_RANGE){return Math.max(num1,num2)+MIN_FRACTION;}
        //System.err.printf("diff=%d ",diff);
        int result=Math.max(num1,num2)+LOGSUM_VALUE_ENTRY[diff];
        //System.err.printf("num1=%f num2=%f diff=%d result=%f\n",Fix2Float(num1),Fix2Float(num2),diff,Fix2Float(result));
        return result;
    }
    
    




    
    public double subComputeReadLikelihoodGivenHaplotypeLog10_exact( final byte[] haplotypeBases,
                                                               final byte[] readBases,
                                                               final byte[] readQuals,
                                                               final byte[] insertionGOP,
                                                               final byte[] deletionGOP,
                                                               final byte[] overallGCP,
                                                               final int hapStartIndex,
                                                               final boolean recacheReadValues,
                                                               final int nextHapStartIndex) {

	
	
	    if (previousHaplotypeBases == null || previousHaplotypeBases.length != haplotypeBases.length) {
                final float initialValue = INITIAL_CONDITION / haplotypeBases.length;
                // set the initial value (free deletions in the beginning) for the first row in the deletion matrix
                for( int j = 0; j < paddedHaplotypeLength; j++ ) {
                    deletionMatrix[0][j] = initialValue;
                }
            }
	    	

            if ( ! constantsAreInitialized || recacheReadValues ) {
                initializeProbabilities(transition, insertionGOP, deletionGOP, overallGCP);

                // note that we initialized the constants
                constantsAreInitialized = true;
            }

            initializePriors(haplotypeBases, readBases, readQuals, hapStartIndex);
            
	    int endI=paddedReadLength-1;
	    for (int i = 1; i < paddedReadLength; i++) {
	        // +1 here is because hapStartIndex is 0-based, but our matrices are 1 based
                //for (int j = hapStartIndex+1; j < paddedHaplotypeLength; j++) {
                for (int j = 1; j < paddedHaplotypeLength; j++) {
                    //Inlined the code from updateCell - helps JIT to detect hotspots and produce good native code
                    //Note that I modified this part to be bit compatible with verilog
                    //matchMatrix[i][j] = prior[i][j] * ( matchMatrix[i - 1][j - 1] * transition[i][matchToMatch] +
                    //        insertionMatrix[i - 1][j - 1] * transition[i][indelToMatch] +
                    //        deletionMatrix[i - 1][j - 1] * transition[i][indelToMatch] );
                    //insertionMatrix[i][j] = matchMatrix[i - 1][j] * transition[i][matchToInsertion] + insertionMatrix[i - 1][j] * transition[i][insertionToInsertion];
                    //deletionMatrix[i][j] = matchMatrix[i][j - 1] * transition[i][matchToDeletion] + deletionMatrix[i][j - 1] * transition[i][deletionToDeletion];
	        		
                    matchMatrix[i][j] =  prior[i][j] * transition[i][matchToMatch]* matchMatrix[i - 1][j - 1]  +
                        		 prior[i][j] *transition[i][indelToMatch]* (insertionMatrix[i - 1][j - 1] +
                        		deletionMatrix[i - 1][j - 1] );
                   	insertionMatrix[i][j] =  matchMatrix[i - 1][j] * transition[i][matchToInsertion] + insertionMatrix[i - 1][j] * transition[i][insertionToInsertion];
                    deletionMatrix[i][j] =  matchMatrix[i][j - 1] * transition[i][matchToDeletion] + deletionMatrix[i][j - 1] * transition[i][deletionToDeletion];
	            }
         }
            // final log probability is the log10 sum of the last element in the Match and Insertion state arrays
            // this way we ignore all paths that ended in deletions! (huge)
            // but we have to sum all the paths ending in the M and I matrices, because they're no longer extended.
            float finalSumProbabilities = (float)0.0;
        for (int j = 1; j < paddedHaplotypeLength; j++) {
                finalSumProbabilities += matchMatrix[endI][j] + insertionMatrix[endI][j];
	    }
	    double result = (double)((float)Math.log10(finalSumProbabilities) - INITIAL_CONDITION_LOG10);
	    //System.err.printf("exact=%f\n",result);
        ////////////////////////////////////Input to Verilog EXACT ///////////////////////////////////////
        //int endJ = paddedHaplotypeLength -1;
        ////if(endI==10 && endJ==142){
        //    System.err.print("read:\n");
        //    for (int i=readBases.length-1; i>=0;i--){
        //        String base = "xx";
        //        //System.err.println((char)readBases[i]);
        //        if(readBases[i]=='A'){
        //            base = "00";
        //        }
        //        if(readBases[i]=='T'){
        //            base = "01";
        //        }
        //        if(readBases[i]=='G'){
        //            base = "10";
        //        }
        //        if(readBases[i]=='C'){
        //            base = "11";
        //        }
        //        //System.err.printf("%d:",i);
        //        System.err.print(base);
        //    }
        //    System.err.print("\n");
	    //    System.err.printf("readlength:%d\n",readBases.length);
        //    System.err.println(Integer.toBinaryString(readBases.length));
	    //    System.err.print("haplotype:\n");
        //    for (int i=haplotypeBases.length-1; i>=0;i--){
        //        String base = "xx";
        //        if(haplotypeBases[i]=='A'){
        //            base = "00";
        //        }
        //        if(haplotypeBases[i]=='T'){
        //            base = "01";
        //        }
        //        if(haplotypeBases[i]=='G'){
        //            base = "10";
        //        }
        //        if(haplotypeBases[i]=='C'){
        //            base = "11";
        //        }
        //        //System.err.printf("%d:",i);
        //        System.err.print(base);
        //    }
        //    System.err.print("\n");
	    //    System.err.printf("haplength:%d\n",haplotypeBases.length);
        //    System.err.println(Integer.toBinaryString(haplotypeBases.length));
	    //    System.err.println("qual");
        //    for (int i=readQuals.length-1; i>=0;i--){
        //        String qs=Integer.toBinaryString(readQuals[i]);
        //        //System.err.println(qs);
        //        while(qs.length()<8){
        //            qs = '0'+qs ;
        //        }
        //        if((int)readQuals[i]<64){
        //            //System.err.print("enter");
        //            System.err.print(qs.substring(2,8));
        //        }else{
        //            System.err.print("111111");
        //        }
        //    }
        //    System.err.print("\n");
        //    System.err.println("initialvalue");
        //    System.err.println(Integer.toBinaryString(Float.floatToIntBits(deletionMatrix[0][0])));
        //    System.err.printf("result:%f"+Integer.toHexString(Float.floatToIntBits(finalSumProbabilities))+"\n",result);
        //    System.err.println(Integer.toBinaryString(Float.floatToIntBits(finalSumProbabilities)));
        //    System.err.println("skipmarker");
        //    System.err.printf("%d %d\n",endI,endJ);
        //    //for (int j = 1; j <=endJ; j++) {
        //    //    for (int i = 1; i <=endI; i++) {
        //    //        if(!(matchMatrix[i][j]==0 && insertionMatrix[i][j]==0 && deletionMatrix[i][j]==0) ){
	    //    //		    System.err.println("Cell:"); 
        //    //            System.err.printf("i %d j %d M "+Integer.toHexString(Float.floatToIntBits(matchMatrix[i][j]))+" I "+Integer.toHexString(Float.floatToIntBits(insertionMatrix[i][j]))+" D "+Integer.toHexString(Float.floatToIntBits(deletionMatrix[i][j]))+"\n" ,i,j);
        //    //            String hapb = "x";
        //    //            String readb = "x";
        //    //            if(haplotypeBases[j-1]=='A'){
        //    //                hapb = "00";
        //    //            }
        //    //            if(haplotypeBases[j-1]=='T'){
        //    //                hapb = "01";
        //    //            }
        //    //            if(haplotypeBases[j-1]=='G'){
        //    //                hapb = "10";
        //    //            }
        //    //            if(haplotypeBases[j-1]=='C'){
        //    //                hapb = "11";
        //    //            }
        //    //            if(readBases[i-1]=='A'){
        //    //                readb = "00";
        //    //            }
        //    //            if(readBases[i-1]=='T'){
        //    //                readb = "01";
        //    //            }
        //    //            if(readBases[i-1]=='G'){
        //    //                readb = "10";
        //    //            }
        //    //            if(readBases[i-1]=='C'){
        //    //                readb = "11";
        //    //            }
        //    //            System.err.print("hapbp:"+hapb+" readbp:"+readb + " q="+readQuals[i-1]+"\n");
        //    //            System.err.print("M: "+Integer.toBinaryString(Float.floatToIntBits(matchMatrix[i][j])));
        //    //            System.err.print(" I: "+Integer.toBinaryString(Float.floatToIntBits(insertionMatrix[i][j])));
        //    //            System.err.print(" D: "+Integer.toBinaryString(Float.floatToIntBits(deletionMatrix[i][j])) + "\n");
        //    //            System.err.print("prior="+Integer.toHexString(Float.floatToIntBits(prior[i][j]))+"\n");
        //    //            System.err.print("probmm:"+Integer.toHexString(Float.floatToIntBits(prior[i][j] * transition[i][matchToMatch])) + " matchim1jm1:"+Integer.toHexString(Float.floatToIntBits(matchMatrix[i - 1][j - 1]))+" match0:"+Integer.toHexString(Float.floatToIntBits(prior[i][j] * transition[i][matchToMatch]* matchMatrix[i - 1][j - 1])));
        //    //            System.err.print("probim:"+Integer.toHexString(Float.floatToIntBits(prior[i][j] * transition[i][indelToMatch])) + " insertionim1jm1:"+Integer.toHexString(Float.floatToIntBits(insertionMatrix[i - 1][j - 1]))+" deletionim1jm1:"+Integer.toHexString(Float.floatToIntBits(deletionMatrix[i - 1][j - 1]))+" match1:"+Integer.toHexString(Float.floatToIntBits(prior[i][j] *transition[i][indelToMatch]* (insertionMatrix[i - 1][j - 1] +deletionMatrix[i - 1][j - 1])))+"\n");
        //    //            System.err.print("matchim1j:"+Integer.toHexString(Float.floatToIntBits(matchMatrix[i - 1][j]))+" m2i:"+Integer.toHexString(Float.floatToIntBits(transition[i][matchToInsertion])));
        //    //            System.err.print("insertionim1j:"+Integer.toHexString(Float.floatToIntBits(insertionMatrix[i - 1][j]))+" i2i:"+Integer.toHexString(Float.floatToIntBits(transition[i][insertionToInsertion])));
        //    //            System.err.print("matchijm1:"+Integer.toHexString(Float.floatToIntBits(matchMatrix[i][j-1]))+" m2d:"+Integer.toHexString(Float.floatToIntBits(transition[i][matchToDeletion])));
        //    //            System.err.print("deletionijm1:"+Integer.toHexString(Float.floatToIntBits(deletionMatrix[i - 1][j]))+" d2d:"+Integer.toHexString(Float.floatToIntBits(transition[i][deletionToDeletion]))+"\n");
        //    //        }
        //    //    
        //    //    }
        //    //}
        //
	    ////}    
        //////////////////////////////////END Input to Verilog EXACT///////////////////////////////////////
	    return result;
   }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] subComputeReadLikelihoodGivenHaplotypeLog10_approximate( final byte[] haplotypeBases,
                                                               final byte[] readBases,
                                                               final byte[] readQuals,
                                                               final byte[] insertionGOP,
                                                               final byte[] deletionGOP,
                                                               final byte[] overallGCP,
                                                               final int hapStartIndex,
                                                               final boolean recacheReadValues,
                                                               final int nextHapStartIndex) {
        //System.err.println("transition");
        //for(int i=0; i<TRANSITION_TABLE.length;i++){
        //    System.err.println(Integer.toBinaryString(TRANSITION_TABLE[i]));
        //}
        //upperlogsumTime=0;	
        //long start = System.nanoTime();
        if ( ! constantsAreInitialized || recacheReadValues ) {
            //initializeProbabilities(transition, insertionGOP, deletionGOP, overallGCP);
                // note that we initialized the constants
            for (int i=0; i<transition.length;i++){
	        	for(int j=0; j<transition[i].length;j++){
	        		app_transition[i][j] = TRANSITION_TABLE[j];
	        	}
	        }
            //constantsAreInitialized = true;
        }

        initializePriors(haplotypeBases, readBases, readQuals, hapStartIndex,app_prior);
        //initializePriors(haplotypeBases, readBases, readQuals, hapStartIndex);
	    int endI=paddedReadLength-1;
	    
	    ////print for debug
        //for(int i=0; i<LOGSUM_VALUE_ENTRY.length;i++){
        //    System.err.printf("logsum %d=%f\n",i,LOGSUM_VALUE_ENTRY[i]);
        //}
        ////end print for debug

        

        //for (int i=0; i<haplotypeBases.length;i++){
        //    System.err.print(haplotypeBases[i]);
        //}
        //System.err.print("\n");
        //System.err.printf("hap len=%d read len=%d\n",haplotypeBases.length, readBases.length);
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    //Redo from the begining of the function using table lookup in log2 domain////////////////////////////////////	
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            //Check prior matrix
	    	
	    int endJ=paddedHaplotypeLength-1;
	    //Copy match,deletion, insertion matrix 
	    //for (int i=0;i<paddedReadLength;i++){
        //        app_deletionMatrix[i][0] = Upper_LOG2(0);
        //        app_insertionMatrix[i][0] = Upper_LOG2(0);
        //        app_matchMatrix[i][0] = Upper_LOG2(0);
	    //	
	    //}
        //    for( int j = 0; j < paddedHaplotypeLength; j++ ) {
        //        app_insertionMatrix[0][j] = Upper_LOG2(0);
        //        app_matchMatrix[0][j] = Upper_LOG2(0);
	    //}

	    if (previousHaplotypeBases == null || previousHaplotypeBases.length != haplotypeBases.length) {
                //final float initialValue_upper_log2 = Upper_LOG2(INITIAL_CONDITION / haplotypeBases.length);
                final int initialValue_upper_log2 = INITIALVALUE_TABLE[haplotypeBases.length];
                //System.err.printf("app_start initialvalue=%d float=%f hap length=%d \n",initialValue_upper_log2,Fix2Float(initialValue_upper_log2),haplotypeBases.length);
                //System.err.print("app_end initialvalue\n");
                
            	for( int j = 0; j < paddedHaplotypeLength; j++ ) {
                		app_deletionMatrix[0][j] = initialValue_upper_log2;
	    	    }
	    }
	    
	    //END Copy match,deletion, insertion matrix 
	    
	    //Construct error matrix
	    //END Construct error matrix
	    
	    
	    //Copy transition matrix and prior matrix
	    //float[][] app_transition=new float[transition.length][transition[0].length];
	    //float[][] app_prior=new float[prior.length][prior[0].length];
        
        //System.err.print("app_end transition\n");
        //System.err.print("app_end priors\n");
        //for (int i=0; i<prior.length;i++){
	    //	for(int j=0; j<prior[i].length;j++){
	    //		app_prior[i][j]=Upper_LOG2_accurate(prior[i][j]);
	    //	}
	    //} 
	    //END Copy transition matrix and prior matrix
	    
	    
	    
	    
	    //Start log2 pairhmm matrix computation
	    //int[][] significance_i=new int [paddedReadLength][paddedHaplotypeLength];
	    //int[][] significance_j=new int [paddedReadLength][paddedHaplotypeLength];
	    //for (int i = 1; i < paddedReadLength; i++) {
        //    	for (int j = hapStartIndex+1; j < paddedHaplotypeLength; j++) {
	    //    		significance_i[i][j]=0;//0 means no consecutive skippable diagonal line
	    //    		//significance_j[i][j]=0;//0 means no consecutive skippable diagonal line
	    //	}
	    //}
	    
	    //////////////////////////////////////Input to verilog////////////////////////////////////
            //PrintFix("log2initial",app_deletionMatrix[0][1],false);
	    //PrintFix("errorthresh",log2current_error,false);

        //System.err.printf("sum_bin lenggth=%d MAX_INTEGER=%d MAX_RANGE_BIN=%d \n",sum_bin.length, MAX_INTEGER,MAX_RANGE_BIN);
	    //for(int index=0; index<sum_bin.length;index++){
	    //	PrintFix("bin",MAX_RANGE_BIN*index,false);
	    //}
	    //////////////////////////////////////END Input to verilog/////////////////////////////////

        //long elapsedTime = System.nanoTime()-start;
        //System.err.printf("time spent before matrix=%d\n",elapsedTime);
	     
        //long start = System.nanoTime();
        int app_mm,app_im,app_dm;
        //int earlyStopRowIndex=(int)Math.floor(0.5*endI);
	    int earlyStopRowIndex=endI;
        //float log10Offset = OFFSET_TABLE[haplotypeBases.length];
        //final int EARLYSTOP_THRESHOLD = Float2Fix((float)(Math.log(Math.pow(10,NOCOMPUTETH+2+log10Offset))/Math.log(2)));
        //final int EARLYSTOP_THRESHOLD = Float2Fix((float)(Math.log(Math.pow(10,NOCOMPUTETH+log10Offset))/Math.log(2)));
        
        final int EARLYSTOP_THRESHOLD = Float2Fix((float)(Math.log(Math.pow(10,NOCOMPUTETH))/Math.log(2)))+Upper_LOG2_accurate(INITIAL_CONDITION_UP);
        
        //int earlyStopIndexJ=16;
        //final boolean doEarlyStop = endJ<=128 ? true : false;
        //System.err.printf("deletion0=%f\n",Fix2Float(app_deletionMatrix[0][1]));
        //int numBigReduce = 0;
        //int preMaxNum = Integer.MIN_VALUE;
        
        for (int i = 1; i <=endI ; i++) {
            //System.err.printf("endI=%d endIndex=%d i=%d\n",endI,endIndex,i);
            int maxNum= Integer.MIN_VALUE;
            for (int j = hapStartIndex+1; j < paddedHaplotypeLength; j++) {
                int tad = app_deletionMatrix[i-1][j-1]+app_transition[i][indelToMatch]+app_prior[i][j];
                int tai = app_insertionMatrix[i-1][j-1]+app_transition[i][indelToMatch]+app_prior[i][j];
                int tb = app_matchMatrix[i-1][j-1]+app_transition[i][matchToMatch]+app_prior[i][j];
                tad = tad < MIN_INTEGER ? MIN_INTEGER : tad;
                tai = tai < MIN_INTEGER ? MIN_INTEGER : tai;
                tb = tb < MIN_INTEGER ? MIN_INTEGER : tb;
                //System.err.printf("i=%d j=%d\n",i,j);
                //System.err.print(" tai:"+Integer.toHexString(tai));
                //System.err.print(" tad:"+Integer.toHexString(tad));
                //System.err.print(" tb:"+Integer.toHexString(tb));
                //int ta = Upper_LOGSUM(tad,tai);
                //System.err.print(" ta:"+Integer.toHexString(ta));
                app_matchMatrix[i][j]=Upper_LOGSUM(Upper_LOGSUM(tad,
                                            tai),
                                            tb);
                 
                //System.err.print(" match:"+Integer.toHexString(app_matchMatrix[i][j]));
                //System.err.print(" anothermatch:"+Integer.toHexString(Upper_LOGSUM(Upper_LOGSUM(tad,
                //                            tai),
                //                            tb)));
                //System.err.print("\n");
                int i0 = app_matchMatrix[i-1][j]+app_transition[i][matchToInsertion];
                int i1 = app_insertionMatrix[i - 1][j] +app_transition[i][insertionToInsertion];
                i0 = i0 < MIN_INTEGER ? MIN_INTEGER : i0;
                i1 = i1 < MIN_INTEGER ? MIN_INTEGER : i1;
                app_insertionMatrix[i][j]=Upper_LOGSUM(i0,i1);

                int d0 = app_matchMatrix[i][j - 1] +app_transition[i][matchToDeletion];
                int d1 = app_deletionMatrix[i][j - 1] +app_transition[i][deletionToDeletion];
                d0 = d0 < MIN_INTEGER ? MIN_INTEGER : d0;
                d1 = d1 < MIN_INTEGER ? MIN_INTEGER : d1;
                app_deletionMatrix[i][j]=Upper_LOGSUM(d0,d1);
	        }
            if(i%16==0){
                for(int j=1;j<=endJ;j++){
	                maxNum=Math.max(app_matchMatrix[i][j],maxNum);	
	                //maxNum=Math.max(app_insertionMatrix[i][j],maxNum);	
	                //maxNum=Math.max(app_deletionMatrix[i][j],maxNum);	
                }
                if(maxNum<EARLYSTOP_THRESHOLD){
                    earlyStopRowIndex = i;
                    //System.err.printf("earlystoprow=%d maxNum=%f "+Integer.toHexString(maxNum)+"\n",i,Fix2Float(maxNum));
                    break;
                }
            }
	    }
        
        
        //Try to stop early and binsum
        //Seperating match,insertion and deletion is really to match to hardware implementation
        int app_finalSumProbabilities_match = MIN_INTEGER;
        int app_finalSumProbabilities_insertion = MIN_INTEGER;
        int app_finalSumProbabilities_deletion = MIN_INTEGER;
        int app_finalSumProbabilities = MIN_INTEGER;
        double finalSumProbabilities_upperbound;
        for (int index=0;index<sum_bin_match.length;index++){
	    	sum_bin_match[index]=MIN_INTEGER;
	    	sum_bin_insertion[index]=MIN_INTEGER;
	    	sum_bin_deletion[index]=MIN_INTEGER;
	    }
	    //For match
        for (int j = 1; j < paddedHaplotypeLength; j++) {
	        int toadd=app_matchMatrix[earlyStopRowIndex][j];
            int index= Math.max((int)Math.floor((float)toadd/MAX_RANGE_BIN)+sum_bin_match.length/2,0);
            sum_bin_match[index]=Upper_LOGSUM(sum_bin_match[index],toadd);
            //index = index<4? index+4 : index-4;
            //System.err.printf("j:%d M:%d "+Integer.toHexString(toadd)+" %f\n",j,index,Fix2Float(toadd));
        }
	    //For insertion
        for (int j = 1; j < paddedHaplotypeLength; j++) {
            int toadd=app_insertionMatrix[earlyStopRowIndex][j];
            int index = Math.max((int)Math.floor((float)toadd/MAX_RANGE_BIN)+sum_bin_insertion.length/2,0);
            sum_bin_insertion[index]=Upper_LOGSUM(sum_bin_insertion[index],toadd);
            //index = index<4? index+4 : index-4;
            //System.err.printf("j:%d I:%d "+Integer.toHexString(toadd)+" %f\n",j,index,Fix2Float(toadd));
        }
	    //For deletion
        if(earlyStopRowIndex!=endI){
            for (int j = 1; j < paddedHaplotypeLength; j++) {
                int toadd=app_deletionMatrix[earlyStopRowIndex][j];
                int index = Math.max((int)Math.floor((float)toadd/MAX_RANGE_BIN)+sum_bin_deletion.length/2,0);
                sum_bin_deletion[index]=Upper_LOGSUM(sum_bin_deletion[index],toadd);
                //System.err.printf("j:%d D:%d "+Integer.toHexString(toadd)+" %f\n",j,index,Fix2Float(toadd));
            }
        }
	    //For match & insertion
	    for (int index=0;index<sum_bin_match.length;index++){
	    	app_finalSumProbabilities_match=Upper_LOGSUM(sum_bin_match[index],app_finalSumProbabilities_match);
	    	app_finalSumProbabilities_insertion=Upper_LOGSUM(sum_bin_insertion[index],app_finalSumProbabilities_insertion);
            //System.err.printf("bins:%d M:"+Integer.toHexString(sum_bin_match[index])+" I:"+Integer.toHexString(sum_bin_insertion[index])+"\n",index);
        }
        app_finalSumProbabilities = Upper_LOGSUM(app_finalSumProbabilities_match,app_finalSumProbabilities_insertion);
        //System.err.print("matchbin:"+Integer.toHexString(app_finalSumProbabilities_match) + " insertionbin:"+Integer.toHexString(app_finalSumProbabilities_insertion)+"deletiontin"+Integer.toHexString(app_finalSumProbabilities_deletion)+"\n");
        //System.err.print("misum:"+Integer.toHexString(app_finalSumProbabilities)+"\n");
	    //For deletion
        if(earlyStopRowIndex!=endI){
	        for (int index=0;index<sum_bin_deletion.length;index++){
	    	    app_finalSumProbabilities_deletion=Upper_LOGSUM(sum_bin_deletion[index],app_finalSumProbabilities_deletion);
                //System.err.printf("bins:%d D:"+Integer.toHexString(sum_bin_deletion[index])+"\n",index);
            }
            app_finalSumProbabilities = Upper_LOGSUM(app_finalSumProbabilities_deletion,app_finalSumProbabilities);
            //System.err.print("midsum:"+Integer.toHexString(app_finalSumProbabilities)+"\n");
        }
        finalSumProbabilities_upperbound=(double)((float)Math.log10((float)Math.pow(2,Fix2Float(app_finalSumProbabilities-Upper_LOG2_accurate(INITIAL_CONDITION_UP)))));
        //finalSumProbabilities_upperbound=(double)((float)Math.log10((float)Math.pow(2,Fix2Float(app_finalSumProbabilities)))-log10Offset);
        
        //System.err.printf("temp=%f upperbound=%f\n",temp_upperbound, finalSumProbabilities_upperbound); 
        
        System.err.printf("interResult: unskipedratio %f unskip %d tot %d result=%f \n",(float)(earlyStopRowIndex)/(endI),earlyStopRowIndex*endJ,endI*endJ,finalSumProbabilities_upperbound);
	    
       int nocompute_th = Float2Fix((float)(Math.log(Math.pow(10,NOCOMPUTETH))/Math.log(2)))+Upper_LOG2_accurate(INITIAL_CONDITION_UP);
       //float reverseth = ((float)Math.log10((float)Math.pow(2,Fix2Float(nocompute_th-Upper_LOG2_accurate(INITIAL_CONDITION_UP)))));
       //System.err.printf("NOCOMPUTE=%f reverse=%f\n",NOCOMPUTETH,(double)reverseth);
       //System.err.printf("EARLYSTOP_THRESHOLD=%f nocompute_th=%f\n",Fix2Float(EARLYSTOP_THRESHOLD),Fix2Float(nocompute_th));
        

	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    //END Redo from the begining of the function using table lookup in log2 domain////////////////////////////////////	
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ////Accurate calcultaion of active bands from the begining of the function based on skip matrix////////////////////	
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double finalSumProbabilities_lowerbound=Double.NEGATIVE_INFINITY;
        int skip_marker_i = 0;
        int skip_marker_j = 0;
        float finalSumProbabilities_approximate = (float)0.0;
        if(finalSumProbabilities_upperbound>NOCOMPUTETH && earlyStopRowIndex==endI){
	        //start= System.nanoTime();
            int skip_index_lastrow=endJ;
	        //Only leave the biggest number method	
	        //Find the biggest number in match in the earlyStopRowIndex. Insertion is not considered.Maybe risky!
	        int max_m=app_matchMatrix[earlyStopRowIndex][endJ];
	        for (int j=endJ; j>=1; j--){
                //System.err.printf("match=%f deletion=%f insertion=%f j=%d\n",Fix2Float(app_matchMatrix[earlyStopRowIndex][j]),Fix2Float(app_deletionMatrix[earlyStopRowIndex][j]),Fix2Float(app_insertionMatrix[earlyStopRowIndex][j]),j);
	        	//System.err.printf("i=%d j=%d match=%f:"+Integer.toHexString(app_matchMatrix[earlyStopRowIndex][j]) + "max_m=%f:"+Integer.toHexString(max_m)+"\n",earlyStopRowIndex,j,Fix2Float(app_matchMatrix[earlyStopRowIndex][j]),Fix2Float(max_m));
                if(app_matchMatrix[earlyStopRowIndex][j] > max_m){
	        		max_m=app_matchMatrix[earlyStopRowIndex][j];
	        		skip_index_lastrow=j;
	        	}
	            //System.err.printf("max=%f match[%d]=%f\n",Fix2Float(max_m),j,Fix2Float(app_matchMatrix[earlyStopRowIndex][j]));
                //System.err.print(Integer.toHexString(max_m));
                //System.err.print(" "+Integer.toHexString(app_matchMatrix[earlyStopRowIndex][j]));
                //System.err.print("\n");
	        }
            //System.err.printf("skip_index_lastrow=%d\n",skip_index_lastrow);
            if(skip_index_lastrow+(endI-earlyStopRowIndex)>endJ){
                skip_index_lastrow = endJ - (endI - earlyStopRowIndex);
                //System.err.printf("adjusted skip_index_lastrow=%d\n",skip_index_lastrow);
            }
	        //check if significance_min_i correspond to the skip_index_lastrow
            skip_marker_i = earlyStopRowIndex;
            skip_marker_j = skip_index_lastrow;
            if(skip_index_lastrow>0){
                int j=skip_index_lastrow;
                //System.err.printf("skip_index_lastrow=%d earlyStopRowIndex=%d\n",skip_index_lastrow,earlyStopRowIndex);
                for(int i=earlyStopRowIndex;i>1;i--){
                    int tad = app_deletionMatrix[i-1][j-1]+app_transition[i][indelToMatch]+app_prior[i][j];
                    int tai = app_insertionMatrix[i-1][j-1]+app_transition[i][indelToMatch]+app_prior[i][j];
                    int tb = app_matchMatrix[i-1][j-1]+app_transition[i][matchToMatch]+app_prior[i][j];
                    tad = tad < MIN_INTEGER ? MIN_INTEGER : tad;
                    tai = tai < MIN_INTEGER ? MIN_INTEGER : tai;
                    tb = tb < MIN_INTEGER ? MIN_INTEGER : tb;
                    int ta = Upper_LOGSUM(tad,tai);

                   if(ta-tb<LOG2CURRENT_ERROR){
	                	skip_marker_i = i-1;
                        skip_marker_j = j-1;
                        //System.err.printf("i=%d j=%d skipI=%d skipJ=%d\n",i-1,j-1,skip_marker_i,skip_marker_j);
                        //System.err.printf("ta=%f tb=%f error=%f\n",Fix2Float(ta),Fix2Float(tb),Fix2Float(LOG2CURRENT_ERROR));
                        //System.err.print("ta:"+Integer.toHexString(ta)+" tb:"+Integer.toHexString(tb)+"\n");
                        j--;
	               }else{
                        break;
                   }
                }
            }else{//handle wierd situation where read length>haplen && stop early && diagnal cannot reach final row
                skip_marker_i = 2;
                skip_marker_j = 2;
            }

	        //int skip_marker_i=significance_i[endI-1][skip_index_lastrow-1]==0 ? endI : significance_i[endI-1][skip_index_lastrow-1];
	        //int skip_marker_j=significance_j[endI-1][skip_index_lastrow-1]==0 ? skip_index_lastrow : significance_j[endI-1][skip_index_lastrow-1];
	        

	        //System.err.printf("earlystoprow=%d skip_marker_i=%d j=%d\n",earlyStopRowIndex,skip_marker_i,skip_marker_j);
	        if (previousHaplotypeBases == null || previousHaplotypeBases.length != haplotypeBases.length) {
                final float initialValue = INITIAL_CONDITION / haplotypeBases.length;
                // set the initial value (free deletions in the beginning) for the first row in the deletion matrix
                for( int j = 0; j < paddedHaplotypeLength; j++ ) {
                    deletionMatrix[0][j] = initialValue;
                }
	        }
            if ( ! constantsAreInitialized || recacheReadValues ) {
                initializeProbabilities(transition, insertionGOP, deletionGOP, overallGCP);
                    // note that we initialized the constants
                constantsAreInitialized = true;
            }

            initializePriors(haplotypeBases, readBases, readQuals, hapStartIndex);
	        
	        	
	        float unskip_count=0;
	        float total_count=endI*endJ;
	        int line = 0;
	        //System.err.printf("hapStartIndex=%d\n",hapStartIndex);
	        //Reset last row value
	        for (int j = 1; j < paddedHaplotypeLength; j++) {
                matchMatrix[endI][j]=0;
	            insertionMatrix[endI][j]=0;
	            deletionMatrix[endI][j]=0;
            }
	        ///////////////////////////////////////Currently only implement combined dm im skip marker/////////////////////////////
	        int combined_marker_i = skip_marker_i ;	
	        int combined_marker_j = skip_marker_j ;
            //Only do this for verilog debug. No need to do this in real sw
            for(int i=1;i<=endI;i++){
                for(int j=1;j<=endJ; j++){
	        	   matchMatrix[i][j]=0;
	        	   insertionMatrix[i][j]=0;
	        	   deletionMatrix[i][j]=0;
                }
            }
            //END Only do this for verilog debug. No need to do this in real sw
            
            //elapsedTime = System.nanoTime()-start;
            //System.err.printf("time spent before lowerbound=%d\n",elapsedTime);
	        //start = System.nanoTime(); 
            if(combined_marker_i == endI && combined_marker_j==endJ){
                for (int i = 1; i <=endI; i++) {
                    for (int j = 1; j <= endJ; j++) {
	            	    matchMatrix[i][j] =  prior[i][j] * transition[i][matchToMatch]* matchMatrix[i - 1][j - 1]  +
                            		 prior[i][j] *transition[i][indelToMatch]* (insertionMatrix[i - 1][j - 1] +
                            		deletionMatrix[i - 1][j - 1] );
                        insertionMatrix[i][j] =  matchMatrix[i - 1][j] * transition[i][matchToInsertion] + insertionMatrix[i - 1][j] * transition[i][insertionToInsertion];
                        deletionMatrix[i][j] =  matchMatrix[i][j - 1] * transition[i][matchToDeletion] + deletionMatrix[i][j - 1] * transition[i][deletionToDeletion];
                    
                    }
                }
            
            }else{
                for (int i = 1; i < paddedReadLength; i++) {
	            	// +1 here is because hapStartIndex is 0-based, but our matrices are 1 based
                    for (int j = 1; j <= endI-combined_marker_i+combined_marker_j; j++) {
	            	    if(i<combined_marker_i && j<combined_marker_j || (i==combined_marker_i+line && j==combined_marker_j+line)) {
	            		    matchMatrix[i][j] =  prior[i][j] * transition[i][matchToMatch]* matchMatrix[i - 1][j - 1]  +
                                		 prior[i][j] *transition[i][indelToMatch]* (insertionMatrix[i - 1][j - 1] +
                                		deletionMatrix[i - 1][j - 1] );
                       		insertionMatrix[i][j] =  matchMatrix[i - 1][j] * transition[i][matchToInsertion] + insertionMatrix[i - 1][j] * transition[i][insertionToInsertion];
                        	deletionMatrix[i][j] =  matchMatrix[i][j - 1] * transition[i][matchToDeletion] + deletionMatrix[i][j - 1] * transition[i][deletionToDeletion];
                            //System.err.printf("i=%d j=%d marker_i=%d marker_%d line=%d match=%.10f insert=%.10f delete=%.10f\n",i,j,combined_marker_i,combined_marker_j,line,matchMatrix[i][j],insertionMatrix[i][j],deletionMatrix[i][j]);
	            		    unskip_count++;
	            		    if(i>=combined_marker_i && j>=combined_marker_j){
	            			    line++;
	            		        
                            }				

	            	    }
	            	    else {
	            	    	matchMatrix[i][j]=0;
	            	    	insertionMatrix[i][j]=0;
	            	    	deletionMatrix[i][j]=0;
	            	    }
	            	}
	            }
            
            }

            // final log probability is the log10 sum of the last element in the Match and Insertion state arrays
            // this way we ignore all paths that ended in deletions! (huge)
            // but we have to sum all the paths ending in the M and I matrices, because they're no longer extended.
	        for (int j = 1; j < paddedHaplotypeLength; j++) {
                finalSumProbabilities_approximate += matchMatrix[endI][j] + insertionMatrix[endI][j];
                //System.err.printf("match=%f insert=%f sum=%f\n",matchMatrix[endI][j],insertionMatrix[endI][j],finalSumProbabilities_approximate);
            }
	        finalSumProbabilities_lowerbound=(double)((float)Math.log10(finalSumProbabilities_approximate) - INITIAL_CONDITION_LOG10);
            //elapsedTime = System.nanoTime()-start;
            //System.err.printf("time spent to reach lowerbound=%d\n",elapsedTime);
	        
            System.err.printf("inaccuracy=%f skip_ratio=%.1f read_length=%d hap_length=%d approximate_result=%f total_count=%f unskip_count=%.1f upperbound=%f diff=%f\n",1.0,1.0-unskip_count/total_count,endI,endJ,finalSumProbabilities_lowerbound,total_count,unskip_count,finalSumProbabilities_upperbound,finalSumProbabilities_upperbound-finalSumProbabilities_lowerbound);
            //System.err.printf("time spent on matrix=%d\n",elapsedTime);
	        
        }else{
            finalSumProbabilities_lowerbound = Double.NEGATIVE_INFINITY;
            int unskip_count=0; 
            int total_count=endI*endJ;
	        System.err.printf("inaccuracy=%f skip_ratio=%.1f read_length=%d hap_length=%d approximate_result=%f total_count=%d unskip_count=%d upperbound=%f diff=%f\n",1.0,1.0-unskip_count/total_count,endI,endJ,finalSumProbabilities_lowerbound,total_count,unskip_count,finalSumProbabilities_upperbound,finalSumProbabilities_upperbound-finalSumProbabilities_lowerbound);
            //System.err.printf("time spent on matrix=%d\n",elapsedTime);
	        
	        //return new double[]{finalSumProbabilities_lowerbound, finalSumProbabilities_upperbound};

        }     
        
        ////////////////////////////////Input to Verilog UPPERBOUND///////////////////////////////////////
        //if(endI==41 && endJ==41){
            System.err.print("read:\n");
            for (int i=readBases.length-1; i>=0;i--){
                String base = "xx";
                //System.err.println((char)readBases[i]);
                if(readBases[i]=='A'){
                    base = "00";
                }
                if(readBases[i]=='T'){
                    base = "01";
                }
                if(readBases[i]=='G'){
                    base = "10";
                }
                if(readBases[i]=='C'){
                    base = "11";
                }
                //System.err.printf("%d:",i);
                System.err.print(base);
            }
            System.err.print("\n");
	        System.err.printf("readlength:%d\n",readBases.length);
            System.err.println(Integer.toBinaryString(readBases.length));
	        System.err.print("haplotype:\n");
            for (int i=haplotypeBases.length-1; i>=0;i--){
                String base = "xx";
                if(haplotypeBases[i]=='A'){
                    base = "00";
                }
                if(haplotypeBases[i]=='T'){
                    base = "01";
                }
                if(haplotypeBases[i]=='G'){
                    base = "10";
                }
                if(haplotypeBases[i]=='C'){
                    base = "11";
                }
                //System.err.printf("%d:",i);
                System.err.print(base);
            }
            System.err.print("\n");
	        System.err.printf("haplength:%d\n",haplotypeBases.length);
            System.err.println(Integer.toBinaryString(haplotypeBases.length));
	        System.err.println("qual");
            for (int i=readQuals.length-1; i>=0;i--){
                String qs=Integer.toBinaryString(readQuals[i]);
                //System.err.println(qs);
                while(qs.length()<8){
                    qs = '0'+qs ;
                }
                if((int)readQuals[i]<64){
                    //System.err.print("enter");
                    System.err.print(qs.substring(2,8));
                }else{
                    System.err.print("111111");
                }
            }
            System.err.print("\n");
            //System.err.println("log2initialvalue");
            //System.err.println(Integer.toBinaryString(app_deletionMatrix[0][0]));
            //System.err.println("earlystop_th");
            //System.err.println(Integer.toBinaryString(EARLYSTOP_THRESHOLD));
            //System.err.println("error_th");
            //System.err.println(Integer.toBinaryString(LOG2CURRENT_ERROR));
            System.err.printf("result_upperbound:%f:"+Integer.toHexString(app_finalSumProbabilities)+"\n",finalSumProbabilities_upperbound);
            System.err.println(Integer.toBinaryString(app_finalSumProbabilities));
            System.err.printf("result_lowerbound:%f:"+Integer.toHexString(Float.floatToIntBits(finalSumProbabilities_approximate))+"\n",finalSumProbabilities_lowerbound);
            System.err.println(Integer.toBinaryString(Float.floatToIntBits(finalSumProbabilities_approximate)));
            System.err.println("skipmarker");
            System.err.printf("%d %d\n",skip_marker_i,skip_marker_j);
            //for (int i = 1; i <=earlyStopRowIndex ; i++) {
            //    for (int j = 1; j < paddedHaplotypeLength; j++) {
            //        System.err.printf("i=%d j=%d\n", i,j);
            //        System.err.print("M:"+Integer.toHexString(app_matchMatrix[i][j]));
            //        System.err.print(" I:"+Integer.toHexString(app_insertionMatrix[i][j]));
            //        System.err.print(" D:"+Integer.toHexString(app_deletionMatrix[i][j])+"\n");
            //        
            //        String hapb = "x";
            //        String readb = "x";
            //        if(haplotypeBases[j-1]=='A'){
            //            hapb = "00";
            //        }
            //        if(haplotypeBases[j-1]=='T'){
            //            hapb = "01";
            //        }
            //        if(haplotypeBases[j-1]=='G'){
            //            hapb = "10";
            //        }
            //        if(haplotypeBases[j-1]=='C'){
            //            hapb = "11";
            //        }
            //        if(readBases[i-1]=='A'){
            //            readb = "00";
            //        }
            //        if(readBases[i-1]=='T'){
            //            readb = "01";
            //        }
            //        if(readBases[i-1]=='G'){
            //            readb = "10";
            //        }
            //        if(readBases[i-1]=='C'){
            //            readb = "11";
            //        }
            //        System.err.print("hap:"+hapb+" read:"+readb+"\n");
            //        //System.err.print("M "+Integer.toHexString(app_matchMatrix[i][j]));
            //        //System.err.print(" I "+Integer.toHexString(app_insertionMatrix[i][j]));
            //        //System.err.print(" D "+Integer.toHexString(app_deletionMatrix[i][j]));
            //        //System.err.print("\n");
            //        
            //        
            //        System.err.printf("M:%f I:%f D:%f\n",Fix2Float(app_matchMatrix[i][j]),Fix2Float(app_insertionMatrix[i][j]),Fix2Float(app_deletionMatrix[i][j]));
            //        int tad = app_deletionMatrix[i-1][j-1]+app_transition[i][indelToMatch]+app_prior[i][j];
            //        int tai = app_insertionMatrix[i-1][j-1]+app_transition[i][indelToMatch]+app_prior[i][j];
            //        int tb = app_matchMatrix[i-1][j-1]+app_transition[i][matchToMatch]+app_prior[i][j];
            //        System.err.print("i_im1jm1:"+Integer.toHexString(app_insertionMatrix[i-1][j-1]));
            //        System.err.print(" d_im1jm1:"+Integer.toHexString(app_deletionMatrix[i-1][j-1]));
            //        System.err.print(" m_im1jm1:"+Integer.toHexString(app_matchMatrix[i-1][j-1]));
            //        System.err.print("tai:"+Integer.toHexString(tai));
            //        System.err.print(" tad:"+Integer.toHexString(tad));
            //        System.err.print(" tb:"+Integer.toHexString(tb));
            //        System.err.print("\n");
            //        System.err.print("i_im1j:"+Integer.toHexString(app_insertionMatrix[i-1][j]));
            //        System.err.print(" m_im1j:"+Integer.toHexString(app_matchMatrix[i-1][j]));
            //        System.err.print("\n");
            //        System.err.print("d_ijm1:"+Integer.toHexString(app_deletionMatrix[i][j-1]));
            //        System.err.print(" m_ijm1:"+Integer.toHexString(app_matchMatrix[i][j-1]));
            //        System.err.print("\n");
            //        tad = tad < MIN_INTEGER ? MIN_INTEGER : tad;
            //        tai = tai < MIN_INTEGER ? MIN_INTEGER : tai;
            //        tb = tb < MIN_INTEGER ? MIN_INTEGER : tb;
            //        System.err.print(" tai:"+Integer.toHexString(tai));
            //        System.err.print(" tad:"+Integer.toHexString(tad));
            //        System.err.print(" tb:"+Integer.toHexString(tb));
            //        System.err.print("\n");
            //        int ta = Upper_LOGSUM(tad,tai);
            //        System.err.print("TA:"+Integer.toHexString(ta));
            //        System.err.print(" TB:"+Integer.toHexString(tb));
            //        System.err.print("\n");
            //        System.err.printf("TA:%f TB:%f\n",Fix2Float(ta),Fix2Float(tb));
            //        System.err.print("logsumtatb="+Integer.toHexString(Upper_LOGSUM(ta,tb))+"\n");
            //        
            //        
            //        //System.err.println("Ileft:"+Integer.toHexString(app_matchMatrix[i-1][j])+" "+Integer.toHexString(app_matchMatrix[i-1][j]+app_transition[i][matchToInsertion])+" Iright:"+Integer.toHexString(app_insertionMatrix[i - 1][j])+" "+Integer.toHexString(app_insertionMatrix[i - 1][j] +app_transition[i][insertionToInsertion]));
            //        //System.err.println("Dleft:"+Integer.toHexString(app_matchMatrix[i][j - 1])+" "+Integer.toHexString(app_matchMatrix[i][j - 1] +app_transition[i][matchToDeletion])+" Dright:"+Integer.toHexString(app_deletionMatrix[i][j - 1] )+" "+Integer.toHexString(app_deletionMatrix[i][j - 1] +app_transition[i][deletionToDeletion]));
            //        System.err.print("MM:"+Integer.toHexString(app_transition[i][matchToMatch]+app_prior[i][j]));
            //        System.err.print(" IM:"+Integer.toHexString(app_transition[i][indelToMatch]+app_prior[i][j])+"\n");
            //        //System.err.printf(" MM:%f IM:%f\n",Fix2Float(app_transition[i][matchToMatch]+app_prior[i][j]),Fix2Float(app_transition[i][indelToMatch]+app_prior[i][j]));
            //        //System.err.print(" II:"+Integer.toHexString(app_transition[i][insertionToInsertion]));
            //        //System.err.print(" DD:"+Integer.toHexString(app_transition[i][deletionToDeletion]));
            //        //System.err.print(" MD:"+Integer.toHexString(app_transition[i][matchToDeletion]));
            //        //System.err.print(" MI:"+Integer.toHexString(app_transition[i][matchToInsertion]));
            //        //System.err.print("\n");
            //   }
            //} 
        
        //}
	      
        //////////////////////////////END Input to Verilog UPPERBOUND///////////////////////////////////////
	    
        
        ////////////////////////////////Input to Verilog LOWERBOUND ///////////////////////////////////////
        //   //System.err.print("transitionm2m: "+Integer.toBinaryString(Float.floatToIntBits(transition[1][matchToMatch]))+"\n");
        //   //System.err.print("transitioni2m: "+Integer.toBinaryString(Float.floatToIntBits(transition[1][indelToMatch]))+"\n");
        //   //System.err.print("transitionm2i: "+Integer.toBinaryString(Float.floatToIntBits(transition[1][matchToInsertion]))+"\n");
        //   //System.err.print("transitioni2i: "+Integer.toBinaryString(Float.floatToIntBits(transition[1][insertionToInsertion]))+"\n");
        //   //System.err.print("transitionm2d: "+Integer.toBinaryString(Float.floatToIntBits(transition[1][matchToDeletion]))+"\n");
        //   //System.err.print("transitiond2d: "+Integer.toBinaryString(Float.floatToIntBits(transition[1][deletionToDeletion]))+"\n");
        //   //for(int i=0;i<64;i++){
        //   //      System.err.println(Integer.toBinaryString(Float.floatToIntBits((float)QualityUtils.qualToErrorProb(i))));
        //   //}
        //if(skip_marker_i>8){
        ////if(endI==10 && endJ==142){
        //    System.err.print("approximate:\n");
        //    System.err.print("read:\n");
        //    for (int i=readBases.length-1; i>=0;i--){
        //        String base = "xx";
        //        //System.err.println((char)readBases[i]);
        //        if(readBases[i]=='A'){
        //            base = "00";
        //        }
        //        if(readBases[i]=='T'){
        //            base = "01";
        //        }
        //        if(readBases[i]=='G'){
        //            base = "10";
        //        }
        //        if(readBases[i]=='C'){
        //            base = "11";
        //        }
        //        //System.err.printf("%d:",i);
        //        System.err.print(base);
        //    }
        //    System.err.print("\n");
	    //    System.err.printf("readlength:%d\n",readBases.length);
        //    System.err.println(Integer.toBinaryString(readBases.length));
	    //    System.err.print("haplotype:\n");
        //    for (int i=haplotypeBases.length-1; i>=0;i--){
        //        String base = "xx";
        //        if(haplotypeBases[i]=='A'){
        //            base = "00";
        //        }
        //        if(haplotypeBases[i]=='T'){
        //            base = "01";
        //        }
        //        if(haplotypeBases[i]=='G'){
        //            base = "10";
        //        }
        //        if(haplotypeBases[i]=='C'){
        //            base = "11";
        //        }
        //        //System.err.printf("%d:",i);
        //        System.err.print(base);
        //    }
        //    System.err.print("\n");
	    //    System.err.printf("haplength:%d\n",haplotypeBases.length);
        //    System.err.println(Integer.toBinaryString(haplotypeBases.length));
	    //    System.err.println("qual");
        //    for (int i=readQuals.length-1; i>=0;i--){
        //        String qs=Integer.toBinaryString(readQuals[i]);
        //        //System.err.println(qs);
        //        while(qs.length()<8){
        //            qs = '0'+qs ;
        //        }
        //        if((int)readQuals[i]<64){
        //            //System.err.print("enter");
        //            System.err.print(qs.substring(2,8));
        //        }else{
        //            System.err.print("111111");
        //        }
        //    }
        //    System.err.print("\n");
        //    System.err.println("initialvalue");
        //    System.err.println(Integer.toBinaryString(Float.floatToIntBits(deletionMatrix[0][0])));
        //    System.err.printf("result:%f"+Integer.toHexString(Float.floatToIntBits(finalSumProbabilities_approximate))+"\n",finalSumProbabilities_lowerbound);
        //    System.err.println(Integer.toBinaryString(Float.floatToIntBits(finalSumProbabilities_approximate)));
        //    System.err.println("skipmarker");
        //    System.err.printf("%d %d\n",skip_marker_i,skip_marker_j);
        //for (int j = 1; j <=endJ; j++) {
        //    for (int i = 1; i <=endI; i++) {
        //        if(!(matchMatrix[i][j]==0 && insertionMatrix[i][j]==0 && deletionMatrix[i][j]==0) ){
	    //		    System.err.println("Cell:"); 
        //            System.err.printf("i %d j %d M "+Integer.toHexString(Float.floatToIntBits(matchMatrix[i][j]))+" I "+Integer.toHexString(Float.floatToIntBits(insertionMatrix[i][j]))+" D "+Integer.toHexString(Float.floatToIntBits(deletionMatrix[i][j]))+"\n" ,i,j);
        //            String hapb = "x";
        //            String readb = "x";
        //            if(haplotypeBases[j-1]=='A'){
        //                hapb = "00";
        //            }
        //            if(haplotypeBases[j-1]=='T'){
        //                hapb = "01";
        //            }
        //            if(haplotypeBases[j-1]=='G'){
        //                hapb = "10";
        //            }
        //            if(haplotypeBases[j-1]=='C'){
        //                hapb = "11";
        //            }
        //            if(readBases[i-1]=='A'){
        //                readb = "00";
        //            }
        //            if(readBases[i-1]=='T'){
        //                readb = "01";
        //            }
        //            if(readBases[i-1]=='G'){
        //                readb = "10";
        //            }
        //            if(readBases[i-1]=='C'){
        //                readb = "11";
        //            }
        //            System.err.print("hapbp:"+hapb+" readbp:"+readb + " q="+readQuals[i-1]+"\n");
        //            System.err.print("M: "+Integer.toBinaryString(Float.floatToIntBits(matchMatrix[i][j])));
        //            System.err.print(" I: "+Integer.toBinaryString(Float.floatToIntBits(insertionMatrix[i][j])));
        //            System.err.print(" D: "+Integer.toBinaryString(Float.floatToIntBits(deletionMatrix[i][j])) + "\n");
        //            System.err.print("prior="+Integer.toHexString(Float.floatToIntBits(prior[i][j]))+"\n");
        //            System.err.print("probmm:"+Integer.toHexString(Float.floatToIntBits(prior[i][j] * transition[i][matchToMatch])) + " matchim1jm1:"+Integer.toHexString(Float.floatToIntBits(matchMatrix[i - 1][j - 1]))+" match0:"+Integer.toHexString(Float.floatToIntBits(prior[i][j] * transition[i][matchToMatch]* matchMatrix[i - 1][j - 1])));
        //            System.err.print("probim:"+Integer.toHexString(Float.floatToIntBits(prior[i][j] * transition[i][indelToMatch])) + " insertionim1jm1:"+Integer.toHexString(Float.floatToIntBits(insertionMatrix[i - 1][j - 1]))+" deletionim1jm1:"+Integer.toHexString(Float.floatToIntBits(deletionMatrix[i - 1][j - 1]))+" match1:"+Integer.toHexString(Float.floatToIntBits(prior[i][j] *transition[i][indelToMatch]* (insertionMatrix[i - 1][j - 1] +deletionMatrix[i - 1][j - 1])))+"\n");
        //            //System.err.printf("matchim1j:"+Integer.toHexString(Float.floatToIntBits(matchMatrix[i - 1][j]))+" m2i:"+Integer.toHexString(Float.floatToIntBits(transition[i][matchToInsertion])));
        //            //System.err.print(" insertionim1j:"+Integer.toHexString(Float.floatToIntBits(insertionMatrix[i - 1][j]))+" i2i:"+Integer.toHexString(Float.floatToIntBits(transition[i][insertionToInsertion])));
        //            //System.err.printf(" iproduct0:"+Integer.toHexString(Float.floatToIntBits(matchMatrix[i - 1][j] * transition[i][matchToInsertion]))+" iproduct1:"+Integer.toHexString(Float.floatToIntBits(insertionMatrix[i - 1][j] * transition[i][insertionToInsertion])));
        //            //System.err.printf("matchijm1:"+Integer.toHexString(Float.floatToIntBits(matchMatrix[i][j-1]))+" m2d:"+Integer.toHexString(Float.floatToIntBits(transition[i][matchToDeletion])));
        //            //System.err.print(" deletionijm1:"+Integer.toHexString(Float.floatToIntBits(deletionMatrix[i - 1][j]))+" d2d:"+Integer.toHexString(Float.floatToIntBits(transition[i][deletionToDeletion])));
        //            //System.err.print(" dproduct0:"+Integer.toHexString(Float.floatToIntBits(matchMatrix[i][j - 1] * transition[i][matchToDeletion]))+" dproduct1:"+Integer.toHexString(Float.floatToIntBits(deletionMatrix[i][j - 1] * transition[i][deletionToDeletion]))+"\n");
        //        }
        //    
        //    }
        //}
        //
        //} 
	    //    
        ////////////////////////////////END Input to Verilog LOWERBOUND///////////////////////////////////////
        
        return new double[]{finalSumProbabilities_lowerbound, finalSumProbabilities_upperbound};
	    //System.err.printf("\n");	
    }


    public double subComputeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                               final byte[] readBases,
                                                               final byte[] readQuals,
                                                               final byte[] insertionGOP,
                                                               final byte[] deletionGOP,
                                                               final byte[] overallGCP,
                                                               final int hapStartIndex,
                                                               final boolean recacheReadValues,
                                                               final int nextHapStartIndex) {

	return 0.0;
    }












    /**
     * Initializes the matrix that holds all the constants related to the editing
     * distance between the read and the haplotype.
     *
     * @param haplotypeBases the bases of the haplotype
     * @param readBases      the bases of the read
     * @param readQuals      the base quality scores of the read
     * @param startIndex     where to start updating the distanceMatrix (in case this read is similar to the previous read)
     */
    void initializePriors(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final int startIndex) {

        // initialize the prior matrix for all combinations of read x haplotype bases
        // Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.
        //System.err.printf("doNotUseTristateCorrection=%b\n",doNotUseTristateCorrection);
        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = 0; j < haplotypeBases.length; j++) {
                final byte y = haplotypeBases[j];
                prior[i+1][j+1] = x == y || x == (byte) 'N' || y == (byte) 'N' ?
                        (float)QualityUtils.qualToProb(qual) : ((float)QualityUtils.qualToErrorProb(qual) / (float)(doNotUseTristateCorrection ? 1.0 : TRISTATE_CORRECTION));
                //System.err.printf("imatrix=%d jmatrix=%d x==y?=%b qualtoprob="+Integer.toHexString(Float.floatToIntBits((float)QualityUtils.qualToProb(qual)))+"prior=%f "+Integer.toHexString(Float.floatToIntBits(prior[i+1][j+1]))+"\n",i+1,j+1,x==y,prior[i+1][j+1]); 
	    }
        }
    }
    
    void initializePriors(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final int startIndex, int[][] app_prior) {

        // initialize the prior matrix for all combinations of read x haplotype bases
        // Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.

        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = 0; j < haplotypeBases.length; j++) {
                final byte y = haplotypeBases[j];

                app_prior[i+1][j+1] =  x == y || x == (byte) 'N' || y == (byte) 'N' ?
                        QUAL2PROB_TABLE[qual] : QUAL2ERROR_DIV3_TABLE[qual];
	    }
        }
    }

    /**
     * Initializes the matrix that holds all the constants related to quality scores.
     *
     * @param insertionGOP   insertion quality scores of the read
     * @param deletionGOP    deletion quality scores of the read
     * @param overallGCP     overall gap continuation penalty
     */
    static void initializeProbabilities(final float[][] transition, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP) {
        PairHMMModel.qualToTransProbs(transition,insertionGOP,deletionGOP,overallGCP);
    }
}
