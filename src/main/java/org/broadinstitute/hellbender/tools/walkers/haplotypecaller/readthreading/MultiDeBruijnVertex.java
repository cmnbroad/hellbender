package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.DeBruijnVertex;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * A DeBruijnVertex that supports multiple copies of the same kmer
 *
 * This is implemented through the same mechanism as SeqVertex, where each
 * created MultiDeBruijnVertex has a unique id assigned upon creation.  Two
 * MultiDeBruijnVertex are equal iff they have the same ID
 */
public final class MultiDeBruijnVertex extends DeBruijnVertex {
    private final static boolean KEEP_TRACK_OF_READS = false;

    // Note that using an AtomicInteger is critical to allow multi-threaded HaplotypeCaller
    private static final AtomicInteger idCounter = new AtomicInteger(0);
    private int id = idCounter.getAndIncrement();

    private final List<String> reads = new LinkedList<>();

    /**
     * Create a new MultiDeBruijnVertex with kmer sequence
     * @param sequence the kmer sequence
     */
    MultiDeBruijnVertex(byte[] sequence) {
        super(sequence);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        MultiDeBruijnVertex that = (MultiDeBruijnVertex) o;

        return id == that.id;
    }

    @Override
    public String toString() {
        return "MultiDeBruijnVertex_id_" + id + "_seq_" + getSequenceString();
    }

    /**
     * Add name information to this vertex for debugging
     *
     * This information will be captured as a list of strings, and displayed in DOT if this
     * graph is written out to disk
     *
     * This functionality is only enabled when KEEP_TRACK_OF_READS is true
     *
     * @param name a non-null string
     */
    protected void addRead(final String name) {
        if ( name == null ) throw new IllegalArgumentException("name cannot be null");
        if ( KEEP_TRACK_OF_READS ) reads.add(name);
    }

    @Override
    public int hashCode() { return id; }

    @Override
    public String additionalInfo() {
        return super.additionalInfo() + (KEEP_TRACK_OF_READS ? (! reads.contains("ref") ? "__" + Utils.join(",", reads) : "") : "");
    }

     int getId() {
        return id;
    }
}
