package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;

/**
 * Generic interface for SmithWaterman calculations
 *
 * This interface allows clients to use a generic SmithWaterman variable, without propagating the specific
 * implementation of SmithWaterman throughout their code:
 *
 * SmithWaterman sw = new SpecificSmithWatermanImplementation(ref, read, params)
 * sw.getCigar()
 * sw.getAlignmentStart2wrt1()
 *
 */
public interface SmithWaterman {

    /**
     * Get the cigar string for the alignment of this SmithWaterman class
     * @return a non-null cigar
     */
    public Cigar getCigar();

    /**
     * Get the starting position of the read sequence in the reference sequence
     * @return a positive integer >= 0
     */
    public int getAlignmentStart2wrt1();
}
