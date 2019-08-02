package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;

import java.util.List;

public class ReblockingGVCFWriter extends GVCFWriter {

    public ReblockingGVCFWriter(final VariantContextWriter underlyingWriter, final List<? extends Number> gqPartitions,
                                final boolean floorBlocks, final CachingIndexedFastaSequenceFile referenceReader,
                                final ReblockingOptions reblockingOptions) {
        super(underlyingWriter, gqPartitions, floorBlocks);
        this.gvcfBlockCombiner = new ReblockingGVCFBlockCombiner(gqPartitions, floorBlocks, referenceReader, reblockingOptions);
    }

    @Override
    public void add(VariantContext vc) {
        gvcfBlockCombiner.submit(vc);
        output();
    }

    public int getVcfOutputEnd() {
        return ((ReblockingGVCFBlockCombiner)gvcfBlockCombiner).getVcfOutputEnd();
    }
}
