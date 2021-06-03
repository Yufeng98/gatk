package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

import java.util.ArrayList;
import java.util.List;

/**
 * This class delegates genotyping to allele count- and ploidy-dependent {@link GenotypeLikelihoodCalculator}s
 * under the assumption that sample genotypes are independent conditional on their population frequencies.
 */
public final class IndependentSampleGenotypesModel {
    private static final int DEFAULT_CACHE_PLOIDY_CAPACITY = 10;
    private static final int DEFAULT_CACHE_ALLELE_CAPACITY = 50;

    private final int cacheAlleleCountCapacity;
    private final int cachePloidyCapacity;
    private GenotypeLikelihoodCalculator[][] likelihoodCalculators;
    private final GenotypeLikelihoodCalculators calculators;

    public IndependentSampleGenotypesModel() { this(DEFAULT_CACHE_PLOIDY_CAPACITY, DEFAULT_CACHE_ALLELE_CAPACITY); }

    /**
     *  Initialize model with given maximum allele count and ploidy for caching
     */
    public IndependentSampleGenotypesModel(final int calculatorCachePloidyCapacity, final int calculatorCacheAlleleCapacity) {
        cachePloidyCapacity = calculatorCachePloidyCapacity;
        cacheAlleleCountCapacity = calculatorCacheAlleleCapacity;
        likelihoodCalculators = new GenotypeLikelihoodCalculator[calculatorCachePloidyCapacity][calculatorCacheAlleleCapacity];
        calculators = new GenotypeLikelihoodCalculators();
    }
    //Original Method
    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(final AlleleList<A> genotypingAlleles, final GenotypingData<A> data) {
        Utils.nonNull(genotypingAlleles, "the allele cannot be null");
        Utils.nonNull(data, "the genotyping data cannot be null");

        final AlleleListPermutation<A> permutation = data.permutation(genotypingAlleles);
        final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper = AlleleLikelihoodMatrixMapper.newInstance(permutation);

        final int sampleCount = data.numberOfSamples();
        final PloidyModel ploidyModel = data.ploidyModel();
        final List<GenotypeLikelihoods> genotypeLikelihoods = new ArrayList<>(sampleCount);
        final int alleleCount = genotypingAlleles.numberOfAlleles();

        GenotypeLikelihoodCalculator likelihoodsCalculator = sampleCount > 0 ? getLikelihoodsCalculator(ploidyModel.samplePloidy(0), alleleCount) : null;
        for (int i = 0; i < sampleCount; i++) {
            final int samplePloidy = ploidyModel.samplePloidy(i);

            // get a new likelihoodsCalculator if this sample's ploidy differs from the previous sample's
            if (samplePloidy != likelihoodsCalculator.ploidy()) {
                likelihoodsCalculator = getLikelihoodsCalculator(samplePloidy, alleleCount);
            }

            final LikelihoodMatrix<A> sampleLikelihoods = alleleLikelihoodMatrixMapper.apply(data.readLikelihoods().sampleMatrix(i));
            genotypeLikelihoods.add(likelihoodsCalculator.genotypeLikelihoods(sampleLikelihoods));
        }
        return new GenotypingLikelihoods<>(genotypingAlleles, ploidyModel, genotypeLikelihoods);
    }
    //Prune Method
    public <A extends Allele> List<GenotypingLikelihoods<A>> calculateLikelihoods(final AlleleList<A> genotypingAlleles, final GenotypingData<A> data_lowerbound,final GenotypingData<A> data_upperbound,final GenotypingData<A> data_exact) {
        Utils.nonNull(genotypingAlleles, "the allele cannot be null");
        Utils.nonNull(data_lowerbound, "the genotyping data cannot be null");

        final AlleleListPermutation<A> permutation = data_lowerbound.permutation(genotypingAlleles);
        final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper = AlleleLikelihoodMatrixMapper.newInstance(permutation);

        final int sampleCount = data_lowerbound.numberOfSamples();
        final PloidyModel ploidyModel = data_lowerbound.ploidyModel();
        final List<GenotypeLikelihoods> genotypeLikelihoods_lo = new ArrayList<>(sampleCount);
        final List<GenotypeLikelihoods> genotypeLikelihoods_up = new ArrayList<>(sampleCount);
        final List<GenotypeLikelihoods> genotypeLikelihoods_exact = new ArrayList<>(sampleCount);
        final int alleleCount = genotypingAlleles.numberOfAlleles();

        GenotypeLikelihoodCalculator likelihoodsCalculator_lo = sampleCount > 0 ? getLikelihoodsCalculator(ploidyModel.samplePloidy(0), alleleCount) : null;
        GenotypeLikelihoodCalculator likelihoodsCalculator_up = sampleCount > 0 ? getLikelihoodsCalculator(ploidyModel.samplePloidy(0), alleleCount) : null;
        GenotypeLikelihoodCalculator likelihoodsCalculator_exact = sampleCount > 0 ? getLikelihoodsCalculator(ploidyModel.samplePloidy(0), alleleCount) : null;
        for (int i = 0; i < sampleCount; i++) {
            final int samplePloidy = ploidyModel.samplePloidy(i);

            // get a new likelihoodsCalculator if this sample's ploidy differs from the previous sample's
            if (samplePloidy != likelihoodsCalculator_lo.ploidy()) {
                likelihoodsCalculator_lo = getLikelihoodsCalculator(samplePloidy, alleleCount);
                likelihoodsCalculator_up = getLikelihoodsCalculator(samplePloidy, alleleCount);
                likelihoodsCalculator_exact = getLikelihoodsCalculator(samplePloidy, alleleCount);
            }

            final LikelihoodMatrix<A> sampleLikelihoods_lowerbound = alleleLikelihoodMatrixMapper.apply(data_lowerbound.readLikelihoods().sampleMatrix(i));
            final LikelihoodMatrix<A> sampleLikelihoods_upperbound = alleleLikelihoodMatrixMapper.apply(data_upperbound.readLikelihoods().sampleMatrix(i));
            final LikelihoodMatrix<A> sampleLikelihoods_exact = alleleLikelihoodMatrixMapper.apply(data_exact.readLikelihoods().sampleMatrix(i));
	        //Sanity check:
	        if(sampleLikelihoods_lowerbound.numberOfReads()!=sampleLikelihoods_upperbound.numberOfReads()||sampleLikelihoods_lowerbound.numberOfReads()!=sampleLikelihoods_exact.numberOfReads() ){
	            System.err.printf("Xiao: walkers/genotyper/IndependentSampleGenotypesModel.java/calculateLikelihoods: Readcount inconsistent\n");
	        }
	        if(sampleLikelihoods_lowerbound.numberOfAlleles()!=sampleLikelihoods_upperbound.numberOfAlleles()||sampleLikelihoods_lowerbound.numberOfAlleles()!=sampleLikelihoods_exact.numberOfAlleles() ){
	            System.err.printf("Xiao: walkers/genotyper/IndependentSampleGenotypesModel.java/calculateLikelihoods: Allelescount inconsistent\n");
	        }
            //end Sanity check
	        genotypeLikelihoods_lo.add(likelihoodsCalculator_lo.genotypeLikelihoods(sampleLikelihoods_lowerbound));
	        genotypeLikelihoods_up.add(likelihoodsCalculator_up.genotypeLikelihoods(sampleLikelihoods_upperbound));
	        genotypeLikelihoods_exact.add(likelihoodsCalculator_exact.genotypeLikelihoods(sampleLikelihoods_exact));
            
            
        }
        final  List<GenotypingLikelihoods<A>> result = new ArrayList<GenotypingLikelihoods<A>>();
        result.add(new GenotypingLikelihoods<>(genotypingAlleles, ploidyModel, genotypeLikelihoods_lo));
        result.add(new GenotypingLikelihoods<>(genotypingAlleles, ploidyModel, genotypeLikelihoods_up));
        result.add(new GenotypingLikelihoods<>(genotypingAlleles, ploidyModel, genotypeLikelihoods_exact));
        return result;
    }

    private GenotypeLikelihoodCalculator getLikelihoodsCalculator(final int samplePloidy, final int alleleCount) {
        if (samplePloidy >= cachePloidyCapacity || alleleCount >= cacheAlleleCountCapacity) {
            return calculators.getInstance(samplePloidy, alleleCount);
        }
        final GenotypeLikelihoodCalculator result = likelihoodCalculators[samplePloidy][alleleCount];
        if (result != null) {
            return result;
        } else {
            final GenotypeLikelihoodCalculator newOne = calculators.getInstance(samplePloidy, alleleCount);
            likelihoodCalculators[samplePloidy][alleleCount] = newOne;
            return newOne;
        }
    }
}
