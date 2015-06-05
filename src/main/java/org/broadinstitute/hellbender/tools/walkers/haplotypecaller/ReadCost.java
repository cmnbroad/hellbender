package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Comparator;

/**
 * A pair read-likelihood (cost).
 */
public class ReadCost {

    /**
     * Reference to the read record this cost makes reference to.
     */
    public final SAMRecord read;

    /**
     * Holds the cost value. Public for convenience, please use with care.
     */
    private double cost;

    /**
     * Create a new read cost object provided the read and the gap extension penalty.
     *
     * @param r the read.
     * @param initialCost the initial cost for the read before any read-segment alignment.
     *
     * @throws NullPointerException if {@code r} is {@code null}.
     * @throws IllegalArgumentException if {@code initialCost} is not a valid likelihood.
     */
    public ReadCost(final SAMRecord r, final double initialCost) {
        if (r == null) throw new NullPointerException();
        if (!MathUtils.goodLog10Probability(initialCost))
            throw new IllegalArgumentException("initial cost must be a finite 0 or negative value (" + initialCost + ")");
        read = r;
        cost = initialCost;
    }


    /**
     * Comparator used to sort ReadCosts
     */
    public static final Comparator<ReadCost> COMPARATOR = new Comparator<ReadCost>() {
        @Override
        public int compare(final ReadCost o1, final ReadCost o2) {
            final String s1 = o1.read.getReadName() + (o1.read.getReadPairedFlag() ? (o1.read.getFirstOfPairFlag() ? "/1" : "/2") : "");
            final String s2 = o2.read.getReadName() + (o2.read.getReadPairedFlag() ? (o2.read.getFirstOfPairFlag() ? "/1" : "/2") : "");
            return s1.compareTo(s2);
        }
    };

    /**
     * Add to the cost.
     * @param value value to add.
     */
    public void addCost(final double value) {
        final double previousCost = cost;
        cost += value;
        if (!MathUtils.goodLog10Probability(cost))
            throw new IllegalArgumentException("invalid log10 likelihood value (" + cost + ") after adding (" + value + ") to (" + previousCost + ")");
    }

    /**
     * Return cost.
     * @return 0 or less.
     */
    public double getCost() {
        return cost;
    }
}
