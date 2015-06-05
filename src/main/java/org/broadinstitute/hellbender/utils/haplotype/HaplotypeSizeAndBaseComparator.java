package org.broadinstitute.hellbender.utils.haplotype;

import java.io.Serializable;
import java.util.Comparator;

/**
 * Compares two haplotypes first by their lengths and then by lexicographic order of their bases.
 */
public final class HaplotypeSizeAndBaseComparator implements Comparator<Haplotype>, Serializable {
    @Override
    public int compare( final Haplotype hap1, final Haplotype hap2 ) {
        if (hap1.getBases().length < hap2.getBases().length)
            return -1;
        else if (hap1.getBases().length > hap2.getBases().length)
            return 1;
        else
            return hap1.getBaseString().compareTo(hap2.getBaseString());
    }
}
