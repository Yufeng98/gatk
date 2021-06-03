package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.collections4.ListUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * HaplotypeCaller's genotyping strategy implementation.
 */
public class HaplotypeCallerGenotypingEngine extends AssemblyBasedCallerGenotypingEngine {

    private static final Logger logger = LogManager.getLogger(HaplotypeCallerGenotypingEngine.class);


    private final IndependentSampleGenotypesModel genotypingModel;

    private final PloidyModel ploidyModel;
    private final ReferenceConfidenceMode referenceConfidenceMode;
    protected final double snpHeterozygosity;
    protected final double indelHeterozygosity;

    private final int maxGenotypeCountToEnumerate;
    private final Map<Integer, Integer> practicalAlleleCountForPloidy = new HashMap<>();

    /**
     * {@inheritDoc}
     * @param configuration {@inheritDoc}
     * @param samples {@inheritDoc}
     * @param doPhysicalPhasing whether to try physical phasing.
     */
    public HaplotypeCallerGenotypingEngine(final HaplotypeCallerArgumentCollection configuration, final SampleList samples,
                                           final AFCalculatorProvider afCalculatorProvider, final boolean doPhysicalPhasing) {
        super(configuration, samples, afCalculatorProvider, doPhysicalPhasing);
        ploidyModel = new HomogeneousPloidyModel(samples,configuration.genotypeArgs.samplePloidy);
        genotypingModel = new IndependentSampleGenotypesModel();
        maxGenotypeCountToEnumerate = configuration.genotypeArgs.MAX_GENOTYPE_COUNT;
        referenceConfidenceMode = configuration.emitReferenceConfidence;
        snpHeterozygosity = configuration.genotypeArgs.snpHeterozygosity;
        indelHeterozygosity = configuration.genotypeArgs.indelHeterozygosity;
    }

    @Override
    protected String callSourceString() {
        return "HC_call";
    }

    @Override
    protected boolean forceKeepAllele(final Allele allele) {
        return allele.equals(Allele.NON_REF_ALLELE,false) || referenceConfidenceMode != ReferenceConfidenceMode.NONE
                || configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES;
    }



    /**
     * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
     * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
     *
     * The list of samples we're working with is obtained from the readLikelihoods
     *
     * @param haplotypes                             Haplotypes to assign likelihoods to
     * @param readLikelihoods                        Map from reads->(haplotypes,likelihoods)
     * @param perSampleFilteredReadList              Map from sample to reads that were filtered after assembly and before calculating per-read likelihoods.
     * @param ref                                    Reference bytes at active region
     * @param refLoc                                 Corresponding active region genome location
     * @param activeRegionWindow                     Active window
     * @param activeAllelesToGenotype                Alleles to genotype
     * @param emitReferenceConfidence whether we should add a &lt;NON_REF&gt; alternative allele to the result variation contexts.
     * @param maxMnpDistance Phased substitutions separated by this distance or less are merged into MNPs.  More than
     *                       two substitutions occuring in the same alignment block (ie the same M/X/EQ CIGAR element)
     *                       are merged until a substitution is separated from the previous one by a greater distance.
     *                       That is, if maxMnpDistance = 1, substitutions at 10,11,12,14,15,17 are partitioned into a MNP
     *                       at 10-12, a MNP at 14-15, and a SNP at 17.  May not be negative.
     * @param withBamOut whether to annotate reads in readLikelihoods for future writing to bamout
     *
     * @return                                       A CalledHaplotypes object containing a list of VC's with genotyped events and called haplotypes
     *
     */
    
    //Original Method
    public CalledHaplotypes assignGenotypeLikelihoods(final List<Haplotype> haplotypes,
                                                      final ReadLikelihoods<Haplotype> readLikelihoods,
                                                      final Map<String, List<GATKRead>> perSampleFilteredReadList,
                                                      final byte[] ref,
                                                      final SimpleInterval refLoc,
                                                      final SimpleInterval activeRegionWindow,
                                                      final FeatureContext tracker,
                                                      final List<VariantContext> activeAllelesToGenotype,
                                                      final boolean emitReferenceConfidence,
                                                      final int maxMnpDistance,
                                                      final SAMFileHeader header,
                                                      final boolean withBamOut) {
        // sanity check input arguments
        Utils.nonEmpty(haplotypes, "haplotypes input should be non-empty and non-null");
        Utils.validateArg(readLikelihoods != null && readLikelihoods.numberOfSamples() > 0, "readLikelihoods input should be non-empty and non-null");
        Utils.validateArg(ref != null && ref.length > 0, "ref bytes input should be non-empty and non-null");
        Utils.nonNull(refLoc, "refLoc must be non-null");
        Utils.validateArg(refLoc.size() == ref.length, " refLoc length must match ref bytes");
        Utils.nonNull(activeRegionWindow, "activeRegionWindow must be non-null");
        Utils.nonNull(activeAllelesToGenotype, "activeAllelesToGenotype must be non-null");
        Utils.validateArg(refLoc.contains(activeRegionWindow), "refLoc must contain activeRegionWindow");
        ParamUtils.isPositiveOrZero(maxMnpDistance, "maxMnpDistance may not be negative.");

        // update the haplotypes so we're ready to call, getting the ordered list of positions on the reference
        // that carry events among the haplotypes
        final SortedSet<Integer> startPosKeySet = decomposeHaplotypesIntoVariantContexts(haplotypes, ref, refLoc, activeAllelesToGenotype, maxMnpDistance);

        // Walk along each position in the key set and create each event to be outputted
        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();
        final int ploidy = configuration.genotypeArgs.samplePloidy;
        final List<Allele> noCallAlleles = GATKVariantContextUtils.noCallAlleles(ploidy);

        if (withBamOut) {
            //add annotations to reads for alignment regions and calling regions
            AssemblyBasedCallerUtils.annotateReadLikelihoodsWithRegions(readLikelihoods, activeRegionWindow);
        }

        for( final int loc : startPosKeySet ) {
            if( loc < activeRegionWindow.getStart() || loc > activeRegionWindow.getEnd() ) {
                continue;
            }

            final List<VariantContext> activeEventVariantContexts;
            if( activeAllelesToGenotype.isEmpty() ) {
                activeEventVariantContexts = getVariantContextsFromActiveHaplotypes(loc, haplotypes, true);
            } else { // we are in GGA mode!
                activeEventVariantContexts = getVariantContextsFromGivenAlleles(loc, activeAllelesToGenotype, true);
            }

            final List<VariantContext> eventsAtThisLocWithSpanDelsReplaced =
                    replaceSpanDels(activeEventVariantContexts,
                            Allele.create(ref[loc - refLoc.getStart()], true),
                            loc);

            VariantContext mergedVC = AssemblyBasedCallerUtils.makeMergedVariantContext(eventsAtThisLocWithSpanDelsReplaced);

            if( mergedVC == null ) {
                continue;
            }

            final Map<Allele, List<Haplotype>> alleleMapper = createAlleleMapper(mergedVC, loc, haplotypes, activeAllelesToGenotype);

            if( configuration.debug && logger != null ) {
                logger.info("Genotyping event at " + loc + " with alleles = " + mergedVC.getAlleles());
            }

            mergedVC = removeAltAllelesIfTooManyGenotypes(ploidy, alleleMapper, mergedVC);

            ReadLikelihoods<Allele> readAlleleLikelihoods = readLikelihoods.marginalize(alleleMapper, new SimpleInterval(mergedVC).expandWithinContig(ALLELE_EXTENSION, header.getSequenceDictionary()));
            if (configuration.isSampleContaminationPresent()) {
                readAlleleLikelihoods.contaminationDownsampling(configuration.getSampleContamination());
            }

            if (emitReferenceConfidence) {
                mergedVC = addNonRefSymbolicAllele(mergedVC);
                readAlleleLikelihoods.addNonReferenceAllele(Allele.NON_REF_ALLELE);
            }

            final GenotypesContext genotypes = calculateGLsForThisEvent(readAlleleLikelihoods, mergedVC, noCallAlleles);
            final VariantContext call = calculateGenotypes(new VariantContextBuilder(mergedVC).genotypes(genotypes).make(), getGLModel(mergedVC), header);
            if( call != null ) {

                readAlleleLikelihoods = prepareReadAlleleLikelihoodsForAnnotation(readLikelihoods, perSampleFilteredReadList,
                        emitReferenceConfidence, alleleMapper, readAlleleLikelihoods, call);

                final VariantContext annotatedCall = makeAnnotatedCall(ref, refLoc, tracker, header, mergedVC, readAlleleLikelihoods, call);
                returnCalls.add( annotatedCall );

                if (withBamOut) {
                    AssemblyBasedCallerUtils.annotateReadLikelihoodsWithSupportedAlleles(call, readAlleleLikelihoods);
                }

                // maintain the set of all called haplotypes
                call.getAlleles().stream().map(alleleMapper::get).filter(Objects::nonNull).forEach(calledHaplotypes::addAll);
            }
        }

        final List<VariantContext> phasedCalls = doPhysicalPhasing ? phaseCalls(returnCalls, calledHaplotypes) : returnCalls;
        return new CalledHaplotypes(phasedCalls, calledHaplotypes);
    }
    
    
    
    
    //Prune Method
    public CalledHaplotypes assignGenotypeLikelihoods(final List<Haplotype> haplotypes,
                                                      final List<ReadLikelihoods<Haplotype>> readLikelihoods,
                                                      final Map<String, List<GATKRead>> perSampleFilteredReadList,
                                                      final byte[] ref,
                                                      final SimpleInterval refLoc,
                                                      final SimpleInterval activeRegionWindow,
                                                      final FeatureContext tracker,
                                                      final List<VariantContext> activeAllelesToGenotype,
                                                      final boolean emitReferenceConfidence,
                                                      final int maxMnpDistance,
                                                      final SAMFileHeader header,
                                                      final boolean withBamOut) {
        // sanity check input arguments
        Utils.nonEmpty(haplotypes, "haplotypes input should be non-empty and non-null");
        Utils.validateArg(readLikelihoods.get(0) != null && readLikelihoods.get(0).numberOfSamples() > 0, "readLikelihoods input should be non-empty and non-null");
        Utils.validateArg(ref != null && ref.length > 0, "ref bytes input should be non-empty and non-null");
        Utils.nonNull(refLoc, "refLoc must be non-null");
        Utils.validateArg(refLoc.size() == ref.length, " refLoc length must match ref bytes");
        Utils.nonNull(activeRegionWindow, "activeRegionWindow must be non-null");
        Utils.nonNull(activeAllelesToGenotype, "activeAllelesToGenotype must be non-null");
        Utils.validateArg(refLoc.contains(activeRegionWindow), "refLoc must contain activeRegionWindow");
        ParamUtils.isPositiveOrZero(maxMnpDistance, "maxMnpDistance may not be negative.");

        // update the haplotypes so we're ready to call, getting the ordered list of positions on the reference
        // that carry events among the haplotypes
        final SortedSet<Integer> startPosKeySet = decomposeHaplotypesIntoVariantContexts(haplotypes, ref, refLoc, activeAllelesToGenotype, maxMnpDistance);

        // Walk along each position in the key set and create each event to be outputted
        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();
        final int ploidy = configuration.genotypeArgs.samplePloidy;
        final List<Allele> noCallAlleles = GATKVariantContextUtils.noCallAlleles(ploidy);

        if (withBamOut) {
            //add annotations to reads for alignment regions and calling regions
            AssemblyBasedCallerUtils.annotateReadLikelihoodsWithRegions(readLikelihoods.get(0), activeRegionWindow);
        }

        boolean recompute_done = false;
        for( final int loc : startPosKeySet ) {
            if( loc < activeRegionWindow.getStart() || loc > activeRegionWindow.getEnd() ) {
                continue;
            }

            final List<VariantContext> activeEventVariantContexts;
            if( activeAllelesToGenotype.isEmpty() ) {
                activeEventVariantContexts = getVariantContextsFromActiveHaplotypes(loc, haplotypes, true);
            } else { // we are in GGA mode!
                activeEventVariantContexts = getVariantContextsFromGivenAlleles(loc, activeAllelesToGenotype, true);
            }

            final List<VariantContext> eventsAtThisLocWithSpanDelsReplaced =
                    replaceSpanDels(activeEventVariantContexts,
                            Allele.create(ref[loc - refLoc.getStart()], true),
                            loc);

            VariantContext mergedVC = AssemblyBasedCallerUtils.makeMergedVariantContext(eventsAtThisLocWithSpanDelsReplaced);

            if( mergedVC == null ) {
                continue;
            }

            final Map<Allele, List<Haplotype>> alleleMapper = createAlleleMapper(mergedVC, loc, haplotypes, activeAllelesToGenotype);

            if( configuration.debug && logger != null ) {
                logger.info("Genotyping event at " + loc + " with alleles = " + mergedVC.getAlleles());
            }

            mergedVC = removeAltAllelesIfTooManyGenotypes(ploidy, alleleMapper, mergedVC);
                
            //Partial recompute should start here with a while loop
            boolean[] moreRecompute = new boolean[2];// [0]:crossbound [1]:passEmit
            moreRecompute[0]=false;
            moreRecompute[1]=false;
            int[] readPT=new int[1];
            readPT[0]=0;

            boolean[] exact_only=new boolean[1];
	        VariantContext call;
            List<ReadLikelihoods<Allele>> readAlleleLikelihoods;
            int redolocus=0;
            do{
                int [] redoSquares = new int[1];
                readAlleleLikelihoods = readLikelihoods.get(0).marginalize(haplotypes, readLikelihoods.get(1), readLikelihoods.get(2),alleleMapper, new SimpleInterval(mergedVC).expandWithinContig(ALLELE_EXTENSION, header.getSequenceDictionary()),moreRecompute[0]||moreRecompute[1],readPT,redoSquares);//[0] lowerbound [1] upperbound [2] exact
                //System.err.printf("Xiao: /walkers/haplotypecaller/HaplotypeCallerGenotypingEngine.java/assignGenotypeLikelihoods after marginalization\n");
                
	            if(moreRecompute[0]){
                    System.err.printf("Xiao:redo total squares = %d reason = crossbound\n", redoSquares[0]);
                    redolocus++;
                    System.err.printf("Xiao:redolocus %d\n",redolocus );
                } 
                if(moreRecompute[1]){
                    System.err.printf("Xiao:redo total squares = %d reason = passEmit\n", redoSquares[0]);
                    redolocus++;
                    System.err.printf("Xiao:redolocus %d\n",redolocus );
                }
                //System.err.printf("readPT[0]=%d exact_only=%b \n",readPT[0],exact_only[0] );
	            if (configuration.isSampleContaminationPresent()) {
                    System.err.print("Xiao: /walkers/haplotypecaller/HaplotypeCallerGenotypingEngine.java/assignGenotypeLikelihoods downsample encountered but not implemented\n");
		            readAlleleLikelihoods.get(0).contaminationDownsampling(configuration.getSampleContamination());
                }

                if (emitReferenceConfidence) {
                    mergedVC = addNonRefSymbolicAllele(mergedVC);
                    System.err.print("Xiao: /walkers/haplotypecaller/HaplotypeCallerGenotypingEngine.java/assignGenotypeLikelihoods emitReferenceConfidence encountered but not implemented\n");
                    readAlleleLikelihoods.get(0).addNonReferenceAllele(Allele.NON_REF_ALLELE);
                }
                List<Allele> bestAlleleList = new ArrayList<Allele>();
                
                final List<GenotypesContext> genotypes = calculateGLsForThisEvent(readAlleleLikelihoods, mergedVC, noCallAlleles,exact_only,moreRecompute,bestAlleleList);
                //System.err.printf("Xiao: /walkers/haplotypecaller/HaplotypeCallerGenotypingEngine.java/assignGenotypeLikelihoods aftercalculateGLmoreRecompute[0]=%b exact_only=%b\n",moreRecompute[0],exact_only[0]);
                //About genotypes:.
                    //if exact_only, the entire result is replaced with exact results. So there is only 1 return result;
                    //if bounds are not crossed && biallic case:            
                    //  genotypes[0]=guaranteed pass if pass
                    //  genotypes[1]=guaranteed fail if fail
                    //  genotypes[2]=exact
                    //if bounds are not crossed && multiple alleles with RA best case: 
                    //  genotypes[0]=guarantee to pass with {RR,XR,XX}
                    //  genotypes[1]=guarantee to pass with {XX,XA,AA}
                    //  genotypes[2]=guarantee to fail with {RR,XR,XX}
                    //  genotypes[3]=guarantee to fail with {XX,XA,AA}
                    //  genotypes[4]=exact
                    //if bounds are not crossed && multiple alleles with AB best case:
                    //  genotypes[0]=guarantee to pass with {RR,XR,XX}
                    //  genotypes[1]=guarantee to pass with {XX,XA,AA}
                    //  genotypes[2]=guarantee to pass with {XX,XB,BB}
                    //  genotypes[3]=guarantee to fail with {RR,XR,XX}
                    //  genotypes[4]=guarantee to fail with {XX,XA,AA}
                    //  genotypes[5]=guarantee to fail with {XX,XB,BB}
                    //  genotypes[6]=exact
                if(!moreRecompute[0]){
                    List<VariantContext> vcList = new ArrayList<VariantContext>();
                    for(int i=0; i<genotypes.size();i++){
                        vcList.add(new VariantContextBuilder(mergedVC).genotypes(genotypes.get(i)).make());
                    }
                    //boolean exact_only = genotypes.size() == 1 ? true : false;
                    //recompute_done = recompute_done || exact_only;
                    //boolean[] recompute_done_ref=new boolean[2];
                    call = calculateGenotypes(exact_only[0], vcList, getGLModel(mergedVC), header,readAlleleLikelihoods,moreRecompute,bestAlleleList);
                    //recompute_done = recompute_done || recompute_done_ref[0];
                }else{
                    call=null;
                }
                //System.err.printf("Xiao: /walkers/haplotypecaller/HaplotypeCallerGenotypingEngine.java/assignGenotypeLikelihoods moreRecompute[0]=%b moreRecompute[1]=%b\n",moreRecompute[0],moreRecompute[1]);
                if(exact_only[0]){
                    moreRecompute[0]=false;
                    moreRecompute[1]=false;
                }
            }while(moreRecompute[0]||moreRecompute[1]); 
            if( call != null ) {

                
                readAlleleLikelihoods.set(0,prepareReadAlleleLikelihoodsForAnnotation(readLikelihoods.get(0), perSampleFilteredReadList,
                        emitReferenceConfidence, alleleMapper, readAlleleLikelihoods.get(0), call));

                final VariantContext annotatedCall = makeAnnotatedCall(ref, refLoc, tracker, header, mergedVC, readAlleleLikelihoods.get(0), call);
                returnCalls.add( annotatedCall );
                ////Print out for debug
                //final VariantContext printcall = annotatedCall;
                //System.err.printf("Xiao walkers/genotyper/GnotypingEngine.java/call!==null\n");
                //for(int s=0;s<printcall.getGenotypes().size();s++){
                //    for(int g=0;g<printcall.getGenotypes().get(s).getPL().length;g++){
                //        System.err.printf("Xiao walkers/genotyper/GnotypingEngine.java/calculateGenotypes PL(%d)[%d]=%d genotype=%s\n",s,g,printcall.getGenotypes().get(s).getPL()[g],printcall.getGenotypes().get(s).getGenotypeString(false));
                //    }
                //}
                ////END Print out for debug
                if (withBamOut) {
                    AssemblyBasedCallerUtils.annotateReadLikelihoodsWithSupportedAlleles(call, readAlleleLikelihoods.get(0));
                }

                // maintain the set of all called haplotypes
                call.getAlleles().stream().map(alleleMapper::get).filter(Objects::nonNull).forEach(calledHaplotypes::addAll);
            }
        }

        final List<VariantContext> phasedCalls = doPhysicalPhasing ? phaseCalls(returnCalls, calledHaplotypes) : returnCalls;
        return new CalledHaplotypes(phasedCalls, calledHaplotypes);
    }

    public CalledHaplotypes assignGenotypeLikelihoods(final List<Haplotype> haplotypes,
                                                      final ReadLikelihoods<Haplotype> readLikelihoods,
                                                      final Map<String, List<GATKRead>> perSampleFilteredReadList,
                                                      final byte[] ref,
                                                      final SimpleInterval refLoc,
                                                      final SimpleInterval activeRegionWindow,
                                                      final FeatureContext tracker,
                                                      final List<VariantContext> activeAllelesToGenotype,
                                                      final boolean emitReferenceConfidence,
                                                      final int maxMnpDistance,
                                                      final SAMFileHeader header) {
        return assignGenotypeLikelihoods(haplotypes,readLikelihoods,perSampleFilteredReadList,ref,refLoc,
                activeRegionWindow,tracker,activeAllelesToGenotype,emitReferenceConfidence,maxMnpDistance,header,false);
    }

    @VisibleForTesting
    static List<VariantContext> replaceSpanDels(final List<VariantContext> eventsAtThisLoc, final Allele refAllele, final int loc) {
        return eventsAtThisLoc.stream().map(vc -> replaceWithSpanDelVC(vc, refAllele, loc)).collect(Collectors.toList());
    }

    @VisibleForTesting
    static VariantContext replaceWithSpanDelVC(final VariantContext variantContext, final Allele refAllele, final int loc) {
        if (variantContext.getStart() == loc) {
            return variantContext;
        } else {
            VariantContextBuilder builder = new VariantContextBuilder(variantContext)
                    .start(loc)
                    .stop(loc)
                    .alleles(Arrays.asList(refAllele, Allele.SPAN_DEL))
                    .genotypes(GenotypesContext.NO_GENOTYPES);
            return builder.make();
        }

    }

    /**
     * If the number of alleles is so high that enumerating all possible genotypes is impractical, as determined by
     * {@link #maxGenotypeCountToEnumerate}, remove alt alleles from the input {@code alleleMapper} that are
     * not well supported by good-scored haplotypes.
     * Otherwise do nothing.
     *
     * Alleles kept are guaranteed to have higher precedence than those removed, where precedence is determined by
     * {@link AlleleScoredByHaplotypeScores}.
     *
     * After the remove operation, entries in map are guaranteed to have the same relative order as they were in the input map,
     * that is, entries will be only be removed but not not shifted relative to each other.
     *  @param ploidy        ploidy of the sample
     * @param alleleMapper  original allele to haplotype map
     */
    private VariantContext removeAltAllelesIfTooManyGenotypes(final int ploidy, final Map<Allele, List<Haplotype>> alleleMapper, final VariantContext mergedVC) {

        final int originalAlleleCount = alleleMapper.size();
        practicalAlleleCountForPloidy.putIfAbsent(ploidy, GenotypeLikelihoodCalculators.computeMaxAcceptableAlleleCount(ploidy, maxGenotypeCountToEnumerate));
        final int practicalAlleleCount = practicalAlleleCountForPloidy.get(ploidy);

        if (originalAlleleCount > practicalAlleleCount) {
            final List<Allele> allelesToKeep = whichAllelesToKeepBasedonHapScores(alleleMapper, practicalAlleleCount);
            alleleMapper.keySet().retainAll(allelesToKeep);
            logger.warn(String.format("Removed alt alleles where ploidy is %d and original allele count is %d, whereas after trimming the allele count becomes %d. Alleles kept are:%s",
                    ploidy, originalAlleleCount, practicalAlleleCount, allelesToKeep));
            return removeExcessAltAllelesFromVC(mergedVC, allelesToKeep);
        } else {
            return mergedVC;
        }
    }

    /**
     * Returns a list of alleles that is a subset of the key set of input map {@code alleleMapper}.
     * The size of the returned list is min({@code desiredNumOfAlleles}, alleleMapper.size()).
     *
     * Alleles kept are guaranteed to have higher precedence than those removed, where precedence is determined by
     * {@link AlleleScoredByHaplotypeScores}.
     *
     * Entries in the returned list are guaranteed to have the same relative order as they were in the input map.
     *
     * @param alleleMapper          original allele to haplotype map
     * @param desiredNumOfAlleles   desired allele count, including ref allele
     */
    @VisibleForTesting
    static List<Allele> whichAllelesToKeepBasedonHapScores(final Map<Allele, List<Haplotype>> alleleMapper,
                                                           final int desiredNumOfAlleles) {

        if(alleleMapper.size() <= desiredNumOfAlleles){
            return alleleMapper.keySet().stream().collect(Collectors.toList());
        }

        final PriorityQueue<AlleleScoredByHaplotypeScores> alleleMaxPriorityQ = new PriorityQueue<>();
        for(final Allele allele : alleleMapper.keySet()){
            final List<Double> hapScores = alleleMapper.get(allele).stream().map(Haplotype::getScore).sorted().collect(Collectors.toList());
            final Double highestScore = hapScores.size() > 0 ? hapScores.get(hapScores.size()-1) : Double.NEGATIVE_INFINITY;
            final Double secondHighestScore = hapScores.size()>1 ? hapScores.get(hapScores.size()-2) : Double.NEGATIVE_INFINITY;

            alleleMaxPriorityQ.add(new AlleleScoredByHaplotypeScores(allele, highestScore, secondHighestScore));
        }

        final Set<Allele> allelesToRetain = new LinkedHashSet<>();
        while(allelesToRetain.size()<desiredNumOfAlleles){
            allelesToRetain.add(alleleMaxPriorityQ.poll().getAllele());
        }
        return alleleMapper.keySet().stream().filter(allelesToRetain::contains).collect(Collectors.toList());
    }

    /**
     * A utility class that provides ordering information, given best and second best haplotype scores.
     * If there's a tie between the two alleles when comparing their best haplotype score, the second best haplotype score
     * is used for breaking the tie. In the case that one allele doesn't have a second best allele, i.e. it has only one
     * supportive haplotype, its second best score is set as {@link Double#NEGATIVE_INFINITY}.
     * In the extremely unlikely cases that two alleles, having the same best haplotype score, neither have a second
     * best haplotype score, or the same second best haplotype score, the order is exactly the same as determined by
     * {@link Allele#compareTo(Allele)}.
     */
    private static final class AlleleScoredByHaplotypeScores implements Comparable<AlleleScoredByHaplotypeScores>{
        private final Allele allele;
        private final Double bestHaplotypeScore;
        private final Double secondBestHaplotypeScore;

        public AlleleScoredByHaplotypeScores(final Allele allele, final Double bestHaplotypeScore, final Double secondBestHaplotypeScore){
            this.allele = allele;
            this.bestHaplotypeScore = bestHaplotypeScore;
            this.secondBestHaplotypeScore = secondBestHaplotypeScore;
        }

        @Override
        public int compareTo(final AlleleScoredByHaplotypeScores other) {

            if(allele.isReference() && other.allele.isNonReference()){
                return -1;
            } else if(allele.isNonReference() && other.allele.isReference()){
                return 1;
            } else if(bestHaplotypeScore > other.bestHaplotypeScore) {
                return -1;
            } else if (bestHaplotypeScore < other.bestHaplotypeScore) {
                return 1;
            } else if (!secondBestHaplotypeScore.equals(other.secondBestHaplotypeScore)) {
                return secondBestHaplotypeScore > other.secondBestHaplotypeScore ? -1 : 1;
            } else {
                return allele.compareTo(other.allele);
            }
        }

        public Allele getAllele(){
            return allele;
        }
    }

    /**
     * Returns an VC that is similar to {@code inputVC} in every aspect except that alleles not in {@code allelesToKeep}
     * are removed in the returned VC.
     * @throws IllegalArgumentException if 1) {@code allelesToKeep} is null or contains null elements; or
     *                                     2) {@code allelesToKeep} doesn't contain a reference allele; or
     *                                     3) {@code allelesToKeep} is not a subset of {@code inputVC.getAlleles()}
     */
    @VisibleForTesting
    static VariantContext removeExcessAltAllelesFromVC(final VariantContext inputVC, final Collection<Allele> allelesToKeep){
        Utils.validateArg(allelesToKeep!=null, "alleles to keep is null");
        Utils.validateArg(!allelesToKeep.contains(null), "alleles to keep contains null elements");
        Utils.validateArg(allelesToKeep.stream().anyMatch(Allele::isReference), "alleles to keep doesn't contain reference allele!");
        Utils.validateArg(inputVC.getAlleles().containsAll(allelesToKeep), "alleles to keep is not a subset of input VC alleles");
        if(inputVC.getAlleles().size() == allelesToKeep.size()) return inputVC;

        final VariantContextBuilder vcb = new VariantContextBuilder(inputVC);
        final List<Allele> originalList = inputVC.getAlleles();
        originalList.retainAll(allelesToKeep);
        vcb.alleles(originalList);
        return vcb.make();
    }

    protected VariantContext makeAnnotatedCall(byte[] ref, SimpleInterval refLoc, FeatureContext tracker, SAMFileHeader header, VariantContext mergedVC, ReadLikelihoods<Allele> readAlleleLikelihoods, VariantContext call) {
        final SimpleInterval locus = new SimpleInterval(mergedVC.getContig(), mergedVC.getStart(), mergedVC.getEnd());
        final SimpleInterval refLocInterval= new SimpleInterval(refLoc);
        final ReferenceDataSource refData = new ReferenceMemorySource(new ReferenceBases(ref, refLocInterval), header.getSequenceDictionary());
        final ReferenceContext referenceContext = new ReferenceContext(refData, locus, refLocInterval);

        final VariantContext untrimmedResult =  annotationEngine.annotateContext(call, tracker, referenceContext, readAlleleLikelihoods, a -> true);
        return call.getAlleles().size() == mergedVC.getAlleles().size() ? untrimmedResult
                : GATKVariantContextUtils.reverseTrimAlleles(untrimmedResult);
    }

    private VariantContext addNonRefSymbolicAllele(final VariantContext mergedVC) {
        final List<Allele> alleleList = ListUtils.union(mergedVC.getAlleles(), Arrays.asList(Allele.NON_REF_ALLELE));
        return new VariantContextBuilder(mergedVC).alleles(alleleList).make();
    }

    /**
     * For a particular event described in inputVC, form PL vector for each sample by looking into allele read map and filling likelihood matrix for each allele
     * @param readLikelihoods          Allele map describing mapping from reads to alleles and corresponding likelihoods
     * @param mergedVC               Input VC with event to genotype
     * @return                       GenotypesContext object wrapping genotype objects with PLs
     */
    //Original Method 
    protected GenotypesContext calculateGLsForThisEvent(final ReadLikelihoods<Allele> readLikelihoods, final VariantContext mergedVC, final List<Allele> noCallAlleles ) {
        Utils.nonNull(readLikelihoods, "readLikelihoods");
        Utils.nonNull(mergedVC, "mergedVC");
        final List<Allele> vcAlleles = mergedVC.getAlleles();
        final AlleleList<Allele> alleleList = readLikelihoods.numberOfAlleles() == vcAlleles.size() ? readLikelihoods : new IndexedAlleleList<>(vcAlleles);
        final GenotypingLikelihoods<Allele> likelihoods = genotypingModel.calculateLikelihoods(alleleList,new GenotypingData<>(ploidyModel,readLikelihoods));
        final int sampleCount = samples.numberOfSamples();
        final GenotypesContext result = GenotypesContext.create(sampleCount);
        for (int s = 0; s < sampleCount; s++) {
            result.add(new GenotypeBuilder(samples.getSample(s)).alleles(noCallAlleles).PL(likelihoods.sampleLikelihoods(s).getAsPLs()).make());
        }
        return result;
    }
    
    
    //Prune Method
        //Note that the return list size can be 1 or 3 or 5 or 7. This is assuming only 1 sample
        //if bounds are crossed, the entire result is replaced with exact results. So there is only 1 return result;
        //if bounds are not crossed && biallic case:            
        //  result[0]=guaranteed pass if pass
        //  result[1]=guaranteed fail if fail
        //  result[2]=exact
        //if bounds are not crossed && multiple alleles with RA best case: 
        //  result[0]=guarantee to pass with {RR,XR,XX}
        //  result[1]=guarantee to pass with {XX,XA,AA}
        //  result[2]=guarantee to fail with {RR,XR,XX}
        //  result[3]=guarantee to fail with {XX,XA,AA}
        //  result[4]=exact
        //if bounds are not crossed && multiple alleles with AB best case:
        //  result[0]=guarantee to pass with {RR,XR,XX}
        //  result[1]=guarantee to pass with {XX,XA,AA}
        //  result[2]=guarantee to pass with {XX,XB,BB}
        //  result[3]=guarantee to fail with {RR,XR,XX}
        //  result[4]=guarantee to fail with {XX,XA,AA}
        //  result[5]=guarantee to fail with {XX,XB,BB}
        //  result[6]=exact
    //Prune Method
    protected List<GenotypesContext> calculateGLsForThisEvent(final List<ReadLikelihoods<Allele>> readLikelihoods, final VariantContext mergedVC, final List<Allele> noCallAlleles,final boolean []exact_only, boolean [] moreRecompute, List<Allele> bestAlleleList) {
        Utils.nonNull(readLikelihoods, "readLikelihoods");
        Utils.nonNull(mergedVC, "mergedVC");
        
        final List<GenotypesContext> result = new ArrayList<GenotypesContext>();
        
        //if this region is already recomputed, then no need to combine bounds or examine for recompute
        if(exact_only[0]){
            result.add(calculateGLsForThisEvent(readLikelihoods.get(2),mergedVC,noCallAlleles));
            return result;
        }

        final List<Allele> vcAlleles = mergedVC.getAlleles();
        final AlleleList<Allele> alleleList = readLikelihoods.get(0).numberOfAlleles() == vcAlleles.size() ? readLikelihoods.get(0) : new IndexedAlleleList<>(vcAlleles);
        final List<GenotypingLikelihoods<Allele>> likelihoods = genotypingModel.calculateLikelihoods(alleleList,new GenotypingData<>(ploidyModel,readLikelihoods.get(0)),new GenotypingData<>(ploidyModel,readLikelihoods.get(1)),new GenotypingData<>(ploidyModel,readLikelihoods.get(2)));//[0] lowerbound [1] upperbound [2] exact
         
        
        
        //Begin combining upper and lower bounds
        //Print for debug:
        //    for(int s=0;s<samples.numberOfSamples();s++){
        //        for(int g=0;g<likelihoods.get(0).sampleLikelihoods(s).getAsVector().length;g++){
        //            System.err.printf("likelihoods_lower[%d](%d %d)=%f \n",g,GenotypeLikelihoods.getAllelePair(g).alleleIndex1,GenotypeLikelihoods.getAllelePair(g).alleleIndex2,likelihoods.get(0).sampleLikelihoods(s).getAsVector()[g]);
        //        }
        //        for(int g=0;g<likelihoods.get(1).sampleLikelihoods(s).getAsVector().length;g++){
        //            System.err.printf("likelihoods_upper[%d](%d %d)=%f\n",g,GenotypeLikelihoods.getAllelePair(g).alleleIndex1,GenotypeLikelihoods.getAllelePair(g).alleleIndex2,likelihoods.get(1).sampleLikelihoods(s).getAsVector()[g]);
        //        }
        //        for(int g=0;g<likelihoods.get(2).sampleLikelihoods(s).getAsVector().length;g++){
        //            System.err.printf("likelihoods_exact[%d](%d %d)=%f\n",g,GenotypeLikelihoods.getAllelePair(g).alleleIndex1,GenotypeLikelihoods.getAllelePair(g).alleleIndex2,likelihoods.get(2).sampleLikelihoods(s).getAsVector()[g]);
        //        }
        //    }
        //End print for debug
       
       
        List<GenotypeLikelihoods> genotypelikelihoods_pass = new ArrayList<GenotypeLikelihoods>();
        List<GenotypeLikelihoods> genotypelikelihoods_fail = new ArrayList<GenotypeLikelihoods>();
        boolean cross_bound = false;
        int s=0;//only support single sample for now

        //identify most likely genotype in lowerbound first:
        double [] lowerboundValues = likelihoods.get(0).sampleLikelihoods(s).getAsVector();
        double maxGL = lowerboundValues[0];
        int maxGL_index = 0;
        for(int g=0; g<lowerboundValues.length;g++){
            if(lowerboundValues[g]>maxGL){
                maxGL = lowerboundValues[g];
                maxGL_index = g;
            }
        }

        //make sure upperbound of the nonmax likelihoods are lower than maxGL
        double [] upperboundValues = likelihoods.get(1).sampleLikelihoods(s).getAsVector();
        for(int g=0;g<upperboundValues.length;g++){
            if(upperboundValues[g] > maxGL && g!=maxGL_index){
                cross_bound = true;
                //System.err.printf("Xiao:from tools/walkers/haplotypecaller/HaplotypeCallerGenotypingENgine.java/calculateGLsForThisEvent maxGL=%f maxGL_index=%d upperbound=%f\n",maxGL,maxGL_index,likelihoods.get(1).sampleLikelihoods(s).getAsVector()[g] );
            }
        }
        
        moreRecompute[0]=cross_bound;
        if(cross_bound){
            final GenotypesContext resultN = GenotypesContext.create(samples.numberOfSamples());
            resultN.add(new GenotypeBuilder(samples.getSample(0)).alleles(noCallAlleles).PL(likelihoods.get(0).sampleLikelihoods(0).getAsPLs()).make());
            result.add(resultN); 
            return result;
        }
        
        //find out allele pair corresponding to maxGL_index
        boolean bestAlleleAB = false;
        int [] bestAllelePairIndex=new int[2];
        GenotypeLikelihoods.getAllelePair(maxGL_index);
        bestAllelePairIndex[0] = GenotypeLikelihoods.getAllelePair(maxGL_index).alleleIndex1;
        bestAllelePairIndex[1] = GenotypeLikelihoods.getAllelePair(maxGL_index).alleleIndex2;
        //System.err.printf("exact_only=%b maxGL_index=%d alleleindex1=%d alleleindex2=%d\n",exact_only[0],maxGL_index,bestAllelePairIndex[0],bestAllelePairIndex[1]);
        bestAlleleList.add(vcAlleles.get(bestAllelePairIndex[0]));
        bestAlleleList.add(vcAlleles.get(bestAllelePairIndex[1]));
        //produce guaranteed pass if pass and guaranteed fail if fail likelihoods
        //biallelic case
        if(likelihoods.get(0).numberOfAlleles()==2){
            double [] log10Likelihoods_pass = likelihoods.get(0).sampleLikelihoods(s).getAsVector().clone();//lowerbound
            double [] log10Likelihoods_fail = likelihoods.get(1).sampleLikelihoods(s).getAsVector().clone();//upperbound
            double RR_lower = log10Likelihoods_pass[0];
            double RR_upper = log10Likelihoods_fail[0];
            log10Likelihoods_pass[0] = RR_upper;
            log10Likelihoods_fail[0] = RR_lower;
            genotypelikelihoods_pass.add(GenotypeLikelihoods.fromLog10Likelihoods(log10Likelihoods_pass));
            genotypelikelihoods_fail.add(GenotypeLikelihoods.fromLog10Likelihoods(log10Likelihoods_fail));
            //System.err.printf("Xiao:from tools/walkers/haplotypecaller/HaplotypeCallerGenotypingENgine.java/calculateGLsForThisEvent: biallelic case\n" );
        }
        //multiple allele case
        else{
            //RA case
            int BestaltAlleleIndex = bestAllelePairIndex[0]>0 ? bestAllelePairIndex[0] : bestAllelePairIndex[1];
            //handle passEmit criteria(1) {RR,RX,XX}
            double [] log10Likelihoods_pass = likelihoods.get(0).sampleLikelihoods(s).getAsVector().clone();//lowerbound
            double [] log10Likelihoods_fail = likelihoods.get(1).sampleLikelihoods(s).getAsVector().clone();//upperbound
            double RR_lower = log10Likelihoods_pass[0];
            double RR_upper = log10Likelihoods_fail[0];
            log10Likelihoods_pass[0] = RR_upper;
            log10Likelihoods_fail[0] = RR_lower;
            genotypelikelihoods_pass.add(GenotypeLikelihoods.fromLog10Likelihoods(log10Likelihoods_pass));
            genotypelikelihoods_fail.add(GenotypeLikelihoods.fromLog10Likelihoods(log10Likelihoods_fail));
            
            //handle passEmit criteria(2.1) {XX,XA,AA}
            log10Likelihoods_pass = likelihoods.get(0).sampleLikelihoods(s).getAsVector().clone();//lowerbound
            log10Likelihoods_fail = likelihoods.get(1).sampleLikelihoods(s).getAsVector().clone();//upperbound
            for(int g=0; g<likelihoods.get(0).sampleLikelihoods(s).getAsVector().length;g++ ){
                int [] allelePairIndex = new int[2];
                allelePairIndex[0]= GenotypeLikelihoods.getAllelePair(g).alleleIndex1;
                allelePairIndex[1]= GenotypeLikelihoods.getAllelePair(g).alleleIndex2;
                if(allelePairIndex[0]!=BestaltAlleleIndex && allelePairIndex[1]!=BestaltAlleleIndex){
                    double XX_upper = log10Likelihoods_fail[g];
                    double XX_lower = log10Likelihoods_pass[g];
                    log10Likelihoods_pass[g] = XX_upper;
                    log10Likelihoods_fail[g] = XX_lower;
                }
            }
            genotypelikelihoods_pass.add(GenotypeLikelihoods.fromLog10Likelihoods(log10Likelihoods_pass));
            genotypelikelihoods_fail.add(GenotypeLikelihoods.fromLog10Likelihoods(log10Likelihoods_fail));
            //System.err.printf("Xiao:from tools/walkers/haplotypecaller/HaplotypeCallerGenotypingENgine.java/calculateGLsForThisEvent:multi allele case with RA\n" );
            //AB case
            if( bestAllelePairIndex[0]>0 && bestAllelePairIndex[1]>0 && bestAllelePairIndex[0]!=bestAllelePairIndex[1] ){
                //handle passEmit criteria(2.2) {XX,XB,BB}
                int SecondBestAlleleIndex = bestAllelePairIndex[1];
                log10Likelihoods_pass = likelihoods.get(0).sampleLikelihoods(s).getAsVector().clone();//lowerbound
                log10Likelihoods_fail = likelihoods.get(1).sampleLikelihoods(s).getAsVector().clone();//upperbound
                for(int g=0; g<likelihoods.get(0).sampleLikelihoods(s).getAsVector().length;g++ ){
                    int [] allelePairIndex = new int[2];
                    allelePairIndex[0]= GenotypeLikelihoods.getAllelePair(g).alleleIndex1;
                    allelePairIndex[1]= GenotypeLikelihoods.getAllelePair(g).alleleIndex2;
                    if(allelePairIndex[0]!=SecondBestAlleleIndex && allelePairIndex[1]!=SecondBestAlleleIndex){
                        double XX_upper = log10Likelihoods_fail[g];
                        double XX_lower = log10Likelihoods_pass[g];
                        log10Likelihoods_pass[g] = XX_upper;
                        log10Likelihoods_fail[g] = XX_lower;
                    }
                }
                genotypelikelihoods_pass.add(GenotypeLikelihoods.fromLog10Likelihoods(log10Likelihoods_pass));
                genotypelikelihoods_fail.add(GenotypeLikelihoods.fromLog10Likelihoods(log10Likelihoods_fail));
                //System.err.printf("Xiao:from tools/walkers/haplotypecaller/HaplotypeCallerGenotypingENgine.java/calculateGLsForThisEvent:multi allele case with AB\n" );
            }
        }
        
        
        //Note that the return list size can be 1 or 3 or 5 or 7. This is assuming only 1 sample
        //if bounds are crossed, the entire result is replaced with exact results. So there is only 1 return result;
        //if bounds are not crossed && biallic case:            
        //  result[0]=guaranteed pass if pass
        //  result[1]=guaranteed fail if fail
        //  result[2]=exact
        //if bounds are not crossed && multiple alleles with RA or AA best case: 
        //  result[0]=guarantee to pass with {RR,XR,XX}
        //  result[1]=guarantee to pass with {XX,XA,AA}
        //  result[2]=guarantee to fail with {RR,XR,XX}
        //  result[3]=guarantee to fail with {XX,XA,AA}
        //  result[4]=exact
        //if bounds are not crossed && multiple alleles with AB best case:
        //  result[0]=guarantee to pass with {RR,XR,XX}
        //  result[1]=guarantee to pass with {XX,XA,AA}
        //  result[2]=guarantee to pass with {XX,XB,BB}
        //  result[3]=guarantee to fail with {RR,XR,XX}
        //  result[4]=guarantee to fail with {XX,XA,AA}
        //  result[5]=guarantee to fail with {XX,XB,BB}
        //  result[6]=exact
        //final int sampleCount = samples.numberOfSamples();
        final int sampleCount = 1;
        exact_only[0]=true;
        for(int i=0; i<genotypelikelihoods_pass.size();i++){
            for(int g=0;g<genotypelikelihoods_pass.get(i).getAsVector().length;g++){
                if(genotypelikelihoods_pass.get(i).getAsVector()[g]!=genotypelikelihoods_fail.get(i).getAsVector()[g]){
                    exact_only[0]=false;
                } 
               
            }
        }
        if(!cross_bound){
            for(int i=0; i<genotypelikelihoods_pass.size();i++){
                //for(int g=0;g<genotypelikelihoods_pass.get(i).getAsVector().length;g++){
                    //System.err.printf("Xiao:tools/walkers/haplotypecaller/HaplotypeCallerGenotypingENgine.java/calculateGLsForThisEvent: likelihood set pass %d: likelihood[%d]=%f\n",i,g,genotypelikelihoods_pass.get(i).getAsVector()[g]);
                //} 
                final GenotypesContext result_pass = GenotypesContext.create(sampleCount);
                result_pass.add(new GenotypeBuilder(samples.getSample(s)).alleles(noCallAlleles).PL(genotypelikelihoods_pass.get(i).getAsPLs()).make());
                result.add(result_pass);
            }
            for(int i=0; i<genotypelikelihoods_fail.size();i++){
                //for(int g=0;g<genotypelikelihoods_fail.get(i).getAsVector().length;g++){
                    //System.err.printf("Xiao:tools/walkers/haplotypecaller/HaplotypeCallerGenotypingENgine.java/calculateGLsForThisEvent: likelihood set fail %d: likelihood[%d]=%f\n",i,g,genotypelikelihoods_fail.get(i).getAsVector()[g]);
                    
                //} 
                final GenotypesContext result_fail = GenotypesContext.create(sampleCount);
                result_fail.add(new GenotypeBuilder(samples.getSample(s)).alleles(noCallAlleles).PL(genotypelikelihoods_fail.get(i).getAsPLs()).make());
                result.add(result_fail);
            }
            //check for exact only
            //ArrayList<Genotype> passList=new ArrayList<Genotype>();
            //ArrayList<Genotype> failList=new ArrayList<Genotype>();
            //Iterator<Genotype> it=result.get(0).iterator();
            //while(it.hasNext()){
            //    passList.add(it.next());
            //}
            //it=result.get(1).iterator();
            //while(it.hasNext()){
            //    failList.add(it.next());
            //}
            //exact_only[0]=true;
            //for(int i=0; i<passList.size();i++){
            //    for(int g=0;g<passList.get(i).getPL().length;g++){
            //        System.err.printf("pass[%d][%d]=%d fail=%d\n",i,g,passList.get(i).getPL()[g],failList.get(i).getPL()[g]);
            //        if(passList.get(i).getPL()[g]!=failList.get(i).getPL()[g]){
            //            exact_only[0]=false;
            //        }
            //    }
            //}
            //System.err.printf("exact_only=%b\n",exact_only[0]);
            ////Print for debug
            //for(int g=0;g<passList.get(0).getPL().length;g++){
            //    System.err.printf("result_pass[%d]=%d\n",g,passList.get(0).getPL()[g]);
            //}
            //for(int g=0;g<failList.get(0).getPL().length;g++){
            //    System.err.printf("result_fail[%d]=%d\n",g,failList.get(0).getPL()[g]);
            //}
            ////end print for debug
        
        }
        final GenotypesContext result_exact = GenotypesContext.create(sampleCount);
        result_exact.add(new GenotypeBuilder(samples.getSample(0)).alleles(noCallAlleles).PL(likelihoods.get(2).sampleLikelihoods(s).getAsPLs()).make());
        
        result.add(result_exact);
        return result;
    }

    /**
     * Removes symbolic events from list of haplotypes
     * @param haplotypes       Input/output list of haplotypes, before/after removal
     */
    // TODO - split into input haplotypes and output haplotypes as not to share I/O arguments
    protected static void cleanUpSymbolicUnassembledEvents( final List<Haplotype> haplotypes ) {
        Utils.nonNull(haplotypes);
        final List<Haplotype> haplotypesToRemove = new ArrayList<>();
        for( final Haplotype h : haplotypes ) {
            for( final VariantContext vc : h.getEventMap().getVariantContexts() ) {
                if( vc.isSymbolic() ) {
                    for( final Haplotype h2 : haplotypes ) {
                        for( final VariantContext vc2 : h2.getEventMap().getVariantContexts() ) {
                            if( vc.getStart() == vc2.getStart() && (vc2.isIndel() || vc2.isMNP()) ) { // unfortunately symbolic alleles can't currently be combined with non-point events
                                haplotypesToRemove.add(h);
                                break;
                            }
                        }
                    }
                }
            }
        }
        haplotypes.removeAll(haplotypesToRemove);
    }

    /**
     * Returns the ploidy-model used by this genotyping engine.
     *
     * @return never {@code null}.
     */
    public PloidyModel getPloidyModel() {
        return ploidyModel;
    }

    // check whether all alleles of a vc, including the ref but excluding the NON_REF allele, are one base in length
    protected static GenotypeLikelihoodsCalculationModel getGLModel(final VariantContext vc) {
        final boolean isSNP = vc.getAlleles().stream().filter(a -> !a.isSymbolic()).allMatch(a -> a.length() == 1);
        return isSNP ? GenotypeLikelihoodsCalculationModel.SNP : GenotypeLikelihoodsCalculationModel.INDEL;
    }
}
