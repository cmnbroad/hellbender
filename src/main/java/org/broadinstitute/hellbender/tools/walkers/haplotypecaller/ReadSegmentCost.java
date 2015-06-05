package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.Route;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.concurrent.atomic.AtomicLong;

/**
 * Path cost indicate the cost (alignment likelihood) of traversing a section of the graph using a segement of a read.
 *
 * <p>A path can be a whole haplotype path as well as just a smaller haplotype segment</p>.
 *
 * <p>We would generate many of this objects for each read. The final likelihood of a read vs each haplotype
 * would be the summation of the path-cost of that read along the corresponding haplotype path.</p>
 */
class ReadSegmentCost {

    /**
     * Reference to the path the cost makes reference to. Public and writable for convenience, notice that this is not
     * a public class.
     */
    public Route<MultiDeBruijnVertex, MultiSampleEdge> path;


    /**
     * Reference to the read the cost makes reference to. Public and writable for convenience, notice that this is not
     * a public class.
     */
    public SAMRecord read;

    /**
     * Holds the cost value. It public and non-final for convenience.
     */
    private double cost;

    /**
     * Caches the path bases (the haplotype segment bases).
     */
    protected byte[] bases;

    /**
     * Construct a new path cost.
     *
     * @param read the corresponding read.
     * @param path the corresponding path.
     * @param cost initial cost estimate. Might be updated later.
     */
    public ReadSegmentCost(final SAMRecord read,
                    final Route<MultiDeBruijnVertex, MultiSampleEdge> path, final double cost) {
        if (read == null)
            throw new IllegalArgumentException("the read provided cannot be null");
        if (path == null)
            throw new IllegalArgumentException("the path provided cannot be null");
        this.read = read;
        this.path = path;
        setCost(cost);
    }

    /**
     * Returns the cost of a read vs a haplotype.
     *
     * @return a valid log10 likelihood
     */
    public double getCost() {
        return cost;
    }

    /**
     * Changes the cost of the current path vs read combination.
     * @param value
     */
    public void setCost(final double value) {
        if (!MathUtils.goodLog10Probability(value))
            throw new IllegalArgumentException("infinite cost are not allowed");
        cost = value;
    }

    /**
     * Used to generate unique identifiers for path cost object.
     */
    private static final AtomicLong pathCostUniqueIdGenerator = new AtomicLong();

    /**
     * Holds the path cost unique identifier.
     */
    private Long uniqueId;

    /**
     * Returns the this path-cost unique identifier.
     * @return never {@code unique}.
     */
    public long uniqueId() {
        if (uniqueId == null)
            uniqueId = pathCostUniqueIdGenerator.incrementAndGet();
        return uniqueId;
    }
}
