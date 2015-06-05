package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import java.io.Serializable;
import java.util.Comparator;

/**
 * Sorts path costs.
 * <p>
 *     Path costs are first sorted by their path base sequence in alphanumerical order.
 * </p>
 *
 * <p>
 *     When these are the same, we consider their unique ids {@link ReadSegmentCost#uniqueId()} to break the tie.
 * </p>
 *
 */
final class ReadSegmentComparator implements Comparator<ReadSegmentCost>, Serializable {

    public static final Comparator<? super ReadSegmentCost> INSTANCE = new ReadSegmentComparator();
    private static final long serialVersionUID = 1l;

    @Override
    public int compare(final ReadSegmentCost o1, final ReadSegmentCost o2) {
        int minLength = Math.min(o1.bases.length, o2.bases.length);
        for (int i = 0; i < minLength; i++) {
            if (o1.bases[i] == o2.bases[i])
                continue;
            else if (o1.bases[i] < o2.bases[i]) {
                return -1;
            } else {
                return 1;
            }
        }
        if (o1.bases.length < o2.bases.length) {
            return -1;
        } else if (o1.bases.length > o2.bases.length) {
            return 1;
        } else {
            return Long.compare(o1.uniqueId(), o2.uniqueId());
        }

    }
}
