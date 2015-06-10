package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.GenomeLoc;

/**
 * The state of an active region walker's isActive call at a specific locus in the genome
 */
public final class ActivityProfileState {
    final private GenomeLoc loc;
    public double isActiveProb;
    public Type resultState;
    public Number resultValue;

    public enum Type {
        NONE,
        HIGH_QUALITY_SOFT_CLIPS
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of isActiveProb
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param isActiveProb the probability of being active (between 0 and 1)
     */
    public ActivityProfileState(final GenomeLoc loc, final double isActiveProb) {
        this(loc, isActiveProb, Type.NONE, null);
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of isActiveProb that maintains some
     * information about the result state and value
     *
     * The only state value in use is HIGH_QUALITY_SOFT_CLIPS, and here the value is interpreted as the number
     * of bp affected by the soft clips.
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param isActiveProb the probability of being active (between 0 and 1)
     */
    public ActivityProfileState(final GenomeLoc loc, final double isActiveProb, final Type resultState, final Number resultValue) {
        // make sure the location of that activity profile is 1
        if ( loc.size() != 1 )
            throw new IllegalArgumentException("Location for an ActivityProfileState must have to size 1 bp but saw " + loc);
        if ( resultValue != null && resultValue.doubleValue() < 0 )
            throw new IllegalArgumentException("Result value isn't null and its < 0, which is illegal: " + resultValue);

        this.loc = loc;
        this.isActiveProb = isActiveProb;
        this.resultState = resultState;
        this.resultValue = resultValue;
    }

    /**
     * The offset of state w.r.t. our current region's start location
     * @param regionStartLoc the start of the region, as a genome loc
     * @return the position of this profile relative to the start of this region
     */
    public int getOffset(final GenomeLoc regionStartLoc) {
        return getLoc().getStart() - regionStartLoc.getStart();
    }


    /**
     * Get the genome loc associated with the ActivityProfileState
     * @return the location of this result
     */
    public GenomeLoc getLoc() {
        return loc;
    }

    @Override
    public String toString() {
        return "ActivityProfileState{" +
                "loc=" + loc +
                ", isActiveProb=" + isActiveProb +
                ", resultState=" + resultState +
                ", resultValue=" + resultValue +
                '}';
    }
}
