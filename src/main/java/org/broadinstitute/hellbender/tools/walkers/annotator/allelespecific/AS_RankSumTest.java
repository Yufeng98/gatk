package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.utils.CompressedDataList;
import org.broadinstitute.hellbender.utils.Histogram;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;

import java.util.*;

/**
 * Allele-specific implementation of rank sum test annotations
 */
public abstract class AS_RankSumTest extends RankSumTest implements ReducibleAnnotation {
    private static final Logger logger = Logger.getLogger(AS_RankSumTest.class);
    public static final String SPLIT_DELIM = "\\|"; //String.split takes a regex, so we need to escape the pipe
    public static final String PRINT_DELIM = "|";
    public static final String RAW_DELIM = ",";
    public static final String REDUCED_DELIM = ",";

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        return annotateRawData(ref, vc, likelihoods);
    }

    /**
     * Generates an annotation by calling the client implementation of getElementForRead(GATKRead read) over each read
     * given its best assigned allele and returns the value of the allele as a double. This data gets condensed into a
     * CompressedDataList (Which operates by runlength encoding) object represented as a string and returns this in a map with its key being the name of the
     * raw annotation.
     *
     * @param ref the reference context for this annotation
     * @param vc the variant context to annotate
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample
     * @return
     */
    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final ReadLikelihoods<Allele> likelihoods ) {
        if ( likelihoods == null) {
            return Collections.emptyMap();
        }

        final AlleleSpecificAnnotationData<CompressedDataList<Integer>> myRawData = initializeNewRawAnnotationData(vc.getAlleles());
        calculateRawData(vc, likelihoods, myRawData);
        Map<Allele, Double> myRankSumStats = calculateRankSum(myRawData.getAttributeMap(), myRawData.getRefAllele());
        final String annotationString = makeRawAnnotationString(vc.getAlleles(),myRankSumStats);
        if (annotationString == null){
            return Collections.emptyMap();
        }
        return Collections.singletonMap(getRawKeyName(), annotationString);
    }

    /**
     * Initializing a AlleleSpecificAnnotationData<CompressedDataList<Integer>> object for annotateRawData() to be used for
     * the per-read data generated by calculateRawData().
     *
     * @param vcAlleles alleles to segment the annotation data on
     * @return A set of CompressedDataLists representing the the values for the reads supporting each allele
     */
    protected AlleleSpecificAnnotationData<CompressedDataList<Integer>> initializeNewRawAnnotationData(final List<Allele> vcAlleles) {
        Map<Allele, CompressedDataList<Integer>> perAlleleValues = new HashMap<>();
        for (Allele a : vcAlleles) {
            perAlleleValues.put(a, new CompressedDataList<>());
        }
        final AlleleSpecificAnnotationData<CompressedDataList<Integer>> ret = new AlleleSpecificAnnotationData<>(vcAlleles, perAlleleValues.toString());
        ret.setAttributeMap(perAlleleValues);
        return ret;
    }

    /**
     * Initializing a AlleleSpecificAnnotationData<CompressedDataList<Integer>> object for annotateRawData().
     * Note: this differs from initializeNewRawAnnotationData() in that it produces a set of histograms to be used for the combined data
     *
     * @param vcAlleles alleles to segment the annotation data on
     * @return A set of Histograms representing the the values for the reads supporting each allele
     */
    private AlleleSpecificAnnotationData<Histogram> initializeNewAnnotationData(final List<Allele> vcAlleles) {
        Map<Allele, Histogram> perAlleleValues = new HashMap<>();
        for (Allele a : vcAlleles) {
            perAlleleValues.put(a, new Histogram());
        }
        final AlleleSpecificAnnotationData<Histogram> ret = new AlleleSpecificAnnotationData<>(vcAlleles, perAlleleValues.toString());
        ret.setAttributeMap(perAlleleValues);
        return ret;
    }

    protected String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, Double> perAlleleValues) {
        String annotationString = "";
        for (int i = 0; i< vcAlleles.size(); i++) {
            if (!vcAlleles.get(i).isReference()) {
                if (i != 0) { //strings will always start with a printDelim because we won't have values for the reference allele, but keep this for consistency with other annotations
                    annotationString += PRINT_DELIM;
                }
                final Double alleleValue = perAlleleValues.get(vcAlleles.get(i));
                //can be null if there are no ref reads
                if (alleleValue != null) {
                    annotationString += outputSingletonValueAsHistogram(alleleValue);
                }
            }
        }
        return annotationString;
    }

    // Generates as CompressedDataList over integer values over each read
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    public void calculateRawData(VariantContext vc, final ReadLikelihoods<Allele> likelihoods, ReducibleAnnotationData myData) {
        if( vc.getGenotypes().getSampleNames().size() != 1) {
            throw new IllegalStateException("Calculating raw data for allele-specific rank sums requires variant context input with exactly one sample, as in a gVCF.");
        }
        if(likelihoods == null) {
            return;
        }

        final int refLoc = vc.getStart();

        final Map<Allele, CompressedDataList<Integer>> perAlleleValues = myData.getAttributeMap();
        for ( final ReadLikelihoods<Allele>.BestAllele bestAllele : likelihoods.bestAlleles() ) {
            if (bestAllele.isInformative() && isUsableRead(bestAllele.read, refLoc)) {
                final OptionalDouble value = getElementForRead(bestAllele.read, refLoc, bestAllele);
                if (value.isPresent() && value.getAsDouble() != INVALID_ELEMENT_FROM_READ && perAlleleValues.containsKey(bestAllele.allele)) {
                    perAlleleValues.get(bestAllele.allele).add((int) value.getAsDouble());
                }
            }
        }
    }

    /**
     * Parses the raw data histograms generated for each allele and outputs the median value for each histogram
     * as a representative double value.
     *
     * @param vc -- contains the final set of alleles, possibly subset by GenotypeGVCFs
     * @param originalVC -- used to get all the alleles for all gVCFs
     * @return
     */
    public  Map<String, Object> finalizeRawData(final VariantContext vc, final VariantContext originalVC) {
        if (!vc.hasAttribute(getRawKeyName())) {
            return new HashMap<>();
        }
        final String rawRankSumData = vc.getAttributeAsString(getRawKeyName(),null);
        if (rawRankSumData == null) {
            return new HashMap<>();
        }
        final Map<String,Object> annotations = new HashMap<>();
        final AlleleSpecificAnnotationData<Histogram> myData = new AlleleSpecificAnnotationData<>(originalVC.getAlleles(), rawRankSumData);
        parseRawDataString(myData);

        final Map<Allele, Double> perAltRankSumResults = calculateReducedData(myData.getAttributeMap(), myData.getRefAllele());
        //shortcut for no ref values
        if (perAltRankSumResults.isEmpty()) {
            return annotations;
        }
        final String annotationString = makeReducedAnnotationString(vc, perAltRankSumResults);
        annotations.put(getKeyNames().get(0), annotationString);
        return annotations;
    }

    // Calculates the median of values per-Allele in the rank-sum tests
    public Map<Allele, Double> calculateReducedData(final Map<Allele, Histogram> perAlleleValues, final Allele ref) {
        final Map<Allele, Double> perAltRankSumResults = new HashMap<>();
        for (final Allele alt : perAlleleValues.keySet()) {
            if (!alt.equals(ref, false) && perAlleleValues.get(alt) != null) {
                perAltRankSumResults.put(alt, perAlleleValues.get(alt).median());
            }
        }
        return perAltRankSumResults;
    }

    /**
     * Parses expected raw data strings per-allele, which are in the form of run length encoded histograms and combines
     * them by summing the counts for each histogram bin, producing a new set histograms for each output allele.
     *
     * @param vcAlleles -- List of alleles over which to calculate RankSum for
     * @param annotationList -- Annotation
     * @return Single element with the annotation name as its key, and a separated list of run length encoded histograms
     * representing the raw data observed for each allele.
     */
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    public Map<String, Object> combineRawData(final List<Allele> vcAlleles, final List<ReducibleAnnotationData<?>> annotationList) {
        //VC already contains merged alleles from ReferenceConfidenceVariantContextMerger
        final ReducibleAnnotationData combinedData = initializeNewAnnotationData(vcAlleles);

        for (final ReducibleAnnotationData currentValue : annotationList) {
            parseRawDataString(currentValue);
            combineAttributeMap(currentValue, combinedData);

        }
        final String annotationString = makeCombinedAnnotationString(vcAlleles, combinedData.getAttributeMap());
        return Collections.singletonMap(getRawKeyName(), annotationString);
    }

    // Parses the raw data string into a Histogram and sets the inputs attribute map accordingly
    protected void parseRawDataString(final ReducibleAnnotationData<Histogram> myData) {
        final String rawDataString = myData.getRawData();
        String rawDataNoBrackets;
        final Map<Allele, Histogram> perAlleleValues = new HashMap<>();
        //Initialize maps
        for (final Allele current : myData.getAlleles()) {
            perAlleleValues.put(current, new Histogram());
        }
        //Map gives back list with []
        if (rawDataString.charAt(0) == '[') {
            rawDataNoBrackets = rawDataString.substring(1, rawDataString.length() - 1);
        }
        else {
            rawDataNoBrackets = rawDataString;
        }
        //TODO handle misformatted annotation field more gracefully
        //rawDataPerAllele is a per-sample list of the rank sum statistic for each allele
        final String[] rawDataPerAllele = rawDataNoBrackets.split(SPLIT_DELIM);
        for (int i=0; i<rawDataPerAllele.length; i++) {
            final String alleleData = rawDataPerAllele[i];
            final Histogram alleleList = perAlleleValues.get(myData.getAlleles().get(i));
            final String[] rawListEntriesAsStringVector = alleleData.split(RAW_DELIM);
            for (int j=0; j<rawListEntriesAsStringVector.length; j+=2) {
                if (!rawListEntriesAsStringVector[j].isEmpty()) {
                    Double value = Double.parseDouble(rawListEntriesAsStringVector[j].trim());
                    if (!value.isNaN() && (!rawListEntriesAsStringVector[j + 1].isEmpty())) {
                        int count = Integer.parseInt(rawListEntriesAsStringVector[j + 1].trim());
                        alleleList.add(value, count);
                    }
                }
            }
        }
        myData.setAttributeMap(perAlleleValues);
        myData.validateAllelesList();
    }


    protected void combineAttributeMap(final ReducibleAnnotationData<Histogram> toAdd, final ReducibleAnnotationData<Histogram> combined) {
        for (final Allele a : combined.getAlleles()) {
            if (toAdd.hasAttribute(a)) {
                final Histogram alleleData = combined.getAttribute(a);
                if (toAdd.getAttribute(a) != null) {
                    alleleData.add(toAdd.getAttribute(a));
                    combined.putAttribute(a, alleleData);
                }
            }
        }
    }

    private String makeReducedAnnotationString(final VariantContext vc, final Map<Allele,Double> perAltRankSumResults) {
        String annotationString = "";
        for (final Allele a : vc.getAlternateAlleles()) {
            if (!annotationString.isEmpty()) {
                annotationString += REDUCED_DELIM;
            }
            if (!perAltRankSumResults.containsKey(a)) {
                logger.warn("ERROR: VC allele not found in annotation alleles -- maybe there was trimming?");
            } else {
                annotationString += String.format("%.3f", perAltRankSumResults.get(a));
            }
        }
        return annotationString;
    }

    protected String makeCombinedAnnotationString(final List<Allele> vcAlleles, final Map<Allele, Histogram> perAlleleValues) {
        String annotationString = "";
        for (int i = 0; i< vcAlleles.size(); i++) {
            if (!vcAlleles.get(i).isReference()) {
                if (i != 0) {    //strings will always start with a printDelim because we won't have values for the reference allele, but keep this for consistency with other annotations
                    annotationString += PRINT_DELIM;
                }
                final Histogram alleleValue = perAlleleValues.get(vcAlleles.get(i));
                //can be null if there are no ref reads
                if (!alleleValue.isEmpty()) {
                    annotationString += alleleValue.toString();
                }
            }
        }
        return annotationString;
    }

    public Map<Allele,Double> calculateRankSum(final Map<Allele, CompressedDataList<Integer>> perAlleleValues, final Allele ref) {
        final Map<Allele, Double> perAltRankSumResults = new HashMap<>();
        //shortcut to not try to calculate rank sum if there are no reads that unambiguously support the ref
        if (perAlleleValues.get(ref).isEmpty())
            return perAltRankSumResults;
        for (final Allele alt : perAlleleValues.keySet()) {
            if (!alt.equals(ref, false)) {
                final MannWhitneyU mannWhitneyU = new MannWhitneyU();
                //load alts (series 1)
                final List<Double> alts = new ArrayList<>();
                for (final Number qual : perAlleleValues.get(alt)) {
                    alts.add((double) qual.intValue());
                }
                //load refs (series 2)
                final List<Double> refs = new ArrayList<>();
                for (final Number qual : perAlleleValues.get(ref)) {
                    refs.add((double) qual.intValue());
                }

                // we are testing that set1 (the alt bases) have lower quality scores than set2 (the ref bases)
                final MannWhitneyU.Result result = mannWhitneyU.test(ArrayUtils.toPrimitive(alts.toArray(new Double[alts.size()])),
                        ArrayUtils.toPrimitive(refs.toArray(new Double[refs.size()])),
                        MannWhitneyU.TestType.FIRST_DOMINATES);
                perAltRankSumResults.put(alt, result.getZ());
            }
        }
        return perAltRankSumResults;
    }

    public String outputSingletonValueAsHistogram(final Double rankSumValue) {
        Histogram h = new Histogram();
        h.add(rankSumValue);
        return h.toString();
    }


}
