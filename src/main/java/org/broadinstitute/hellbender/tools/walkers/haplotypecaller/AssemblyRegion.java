package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * Interface that describes the region of the genome that gets assembled by the local assembly engine.
 */
public final class AssemblyRegion {

    /**
     * The reads included in this active region.  May be empty upon creation, and expand / contract
     * as reads are added or removed from this region.
     */
    private final List<SAMRecord> reads = new ArrayList<>();

    /**
     * An ordered list (by genomic coordinate) of the ActivityProfileStates that went
     * into this active region.  May be empty, which says that no supporting states were
     * provided when this region was created.
     */
    private final List<ActivityProfileState> supportingStates;

    /**
     * The raw span of this active region, not including the active region extension
     */
    private final GenomeLoc activeRegionLoc;

    /**
     * The span of this active region on the genome, including the active region extension
     */
    private final GenomeLoc extendedLoc;

    /**
     * The extension, in bp, of this active region.
     */
    private final int extension;

    /**
     * A genomeLocParser so we can create genomeLocs
     */
    private final GenomeLocParser genomeLocParser;


    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
    private final boolean isActive;

    /**
     * The span of this active region, including the bp covered by all reads in this
     * region.  This union of extensionLoc and the loc of all reads in this region.
     *
     * Must be at least as large as extendedLoc, but may be larger when reads
     * partially overlap this region.
     */
    private GenomeLoc spanIncludingReads;


    /**
     * Indicates whether the active region has been finalized
     */
    private boolean hasBeenFinalized;

    /**
     * Create a new ActiveRegion containing no reads
     *
     * @param activeRegionLoc the span of this active region
     * @param supportingStates the states that went into creating this region, or null / empty if none are available.
     *                         If not empty, must have exactly one state for each bp in activeRegionLoc
     * @param isActive indicates whether this is an active region, or an inactve one
     * @param genomeLocParser a non-null parser to let us create new genome locs
     * @param extension the active region extension to use for this active region
     */
    public AssemblyRegion( final GenomeLoc activeRegionLoc, final List<ActivityProfileState> supportingStates, final boolean isActive, final GenomeLocParser genomeLocParser, final int extension ) {
        if ( activeRegionLoc == null ) throw new IllegalArgumentException("activeRegionLoc cannot be null");
        if ( activeRegionLoc.size() == 0 ) throw new IllegalArgumentException("Active region cannot be of zero size, but got " + activeRegionLoc);
        if ( genomeLocParser == null ) throw new IllegalArgumentException("genomeLocParser cannot be null");
        if ( extension < 0 ) throw new IllegalArgumentException("extension cannot be < 0 but got " + extension);

        this.activeRegionLoc = activeRegionLoc;
        this.supportingStates = supportingStates == null ? Collections.<ActivityProfileState>emptyList() : new ArrayList<>(supportingStates);
        this.isActive = isActive;
        this.genomeLocParser = genomeLocParser;
        this.extension = extension;
        this.extendedLoc = genomeLocParser.createGenomeLocOnContig(activeRegionLoc.getContig(), activeRegionLoc.getStart() - extension, activeRegionLoc.getStop() + extension);
        this.spanIncludingReads = extendedLoc;

        if ( ! this.supportingStates.isEmpty() ) {
            if ( this.supportingStates.size() != activeRegionLoc.size() )
                throw new IllegalArgumentException("Supporting states wasn't empty but it doesn't have exactly one state per bp in the active region: states " + this.supportingStates.size() + " vs. bp in region = " + activeRegionLoc.size());
            GenomeLoc lastStateLoc = null;
            for ( final ActivityProfileState state : this.supportingStates ) {
                if ( lastStateLoc != null ) {
                    if ( state.getLoc().getStart() != lastStateLoc.getStart() + 1 || state.getLoc().getContigIndex() != lastStateLoc.getContigIndex())
                        throw new IllegalArgumentException("Supporting state has an invalid sequence: last state was " + lastStateLoc + " but next state was " + state);
                }
                lastStateLoc = state.getLoc();
            }
        }
    }

    /**
     * Simple interface to create an active region that isActive without any profile state
     */
    public AssemblyRegion( final GenomeLoc activeRegionLoc, final GenomeLocParser genomeLocParser, final int extension ) {
        this(activeRegionLoc, Collections.<ActivityProfileState>emptyList(), true, genomeLocParser, extension);
    }

    @Override
    public String toString() {
        return "AssemblyRegion "  + activeRegionLoc.toString() + " active?=" + isActive + " nReads=" + reads.size();
    }

    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
    public boolean isActive() {
        return isActive;
    }

    /**
     * Get the span of this active region including the extension value
     * @return a non-null GenomeLoc
     */
    public GenomeLoc getExtendedSpan() { return extendedLoc; }

    /**
     * Get the raw span of this active region (excluding the extension)
     * @return a non-null genome loc
     */
    public GenomeLoc getSpan() { return activeRegionLoc; }

    /**
     * Get an unmodifiable list of reads currently in this active region.
     *
     * The reads are sorted by their coordinate position
    */
    public List<SAMRecord> getReads(){
        return Collections.unmodifiableList(reads);
    }

    /**
     * Trim this active to just the span, producing a new active region without any reads that has only
     * the extent of newExtend intersected with the current extent
     * @param span the new extend of the active region we want
     * @param extension the extension size we want for the newly trimmed active region
     * @return a non-null, empty active region
     */
    public AssemblyRegion trim(final GenomeLoc span, final int extension) {
        if ( span == null ) throw new IllegalArgumentException("Active region extent cannot be null");
        if ( extension < 0) throw new IllegalArgumentException("the extension size must be 0 or greater");
        final int extendStart = Math.max(1,span.getStart() - extension);
        final int maxStop = genomeLocParser.getSequenceDictionary().getSequence(span.getContigIndex()).getSequenceLength();
        final int extendStop = Math.min(span.getStop() + extension, maxStop);
        final GenomeLoc extendedSpan = genomeLocParser.createGenomeLoc(span.getContig(), extendStart, extendStop);
        return trim(span, extendedSpan);

//TODO - Inconsiste support of substates trimming. Check lack of consistency!!!!
//        final GenomeLoc subLoc = getLocation().intersect(span);
//        final int subStart = subLoc.getStart() - getLocation().getStart();
//        final int subEnd = subStart + subLoc.size();
//        final List<ActivityProfileState> subStates = supportingStates.isEmpty() ? supportingStates : supportingStates.subList(subStart, subEnd);
//        return new ActiveRegion( subLoc, subStates, isActive, genomeLocParser, extension );

    }

    public AssemblyRegion trim(final GenomeLoc span) {
        return trim(span, span);
    }

    /**
     * Trim this active to no more than the span, producing a new active region with properly trimmed reads that
     * attempts to provide the best possible representation of this active region covering the span.
     *
     * The challenge here is that span may (1) be larger than can be represented by this active region
     * + its original extension and (2) the extension must be symmetric on both sides.  This algorithm
     * therefore determines how best to represent span as a subset of the span of this
     * region with a padding value that captures as much of the span as possible.
     *
     * For example, suppose this active region is
     *
     * Active:    100-200 with extension of 50, so that the true span is 50-250
     * NewExtent: 150-225 saying that we'd ideally like to just have bases 150-225
     *
     * Here we represent the active region as a active region from 150-200 with 25 bp of padding.
     *
     * The overall constraint is that the active region can never exceed the original active region, and
     * the extension is chosen to maximize overlap with the desired region
     *
     * @param span the new extend of the active region we want
     * @return a non-null, empty active region
     */
    public AssemblyRegion trim(final GenomeLoc span, final GenomeLoc extendedSpan) {
        if ( span == null ) throw new IllegalArgumentException("Active region extent cannot be null");
        if ( extendedSpan == null ) throw new IllegalArgumentException("Active region extended span cannot be null");
        if ( ! extendedSpan.containsP(span)) {
            throw new IllegalArgumentException("The requested extended must fully contain the requested span");
        }

        final GenomeLoc subActive = getSpan().intersect(span);
        final int requiredOnRight = Math.max(extendedSpan.getStop() - subActive.getStop(), 0);
        final int requiredOnLeft = Math.max(subActive.getStart() - extendedSpan.getStart(), 0);
        final int requiredExtension = Math.min(Math.max(requiredOnLeft, requiredOnRight), getExtension());

        final AssemblyRegion result = new AssemblyRegion( subActive, Collections.<ActivityProfileState>emptyList(), isActive, genomeLocParser, requiredExtension );

        final List<SAMRecord> myReads = getReads();
        final GenomeLoc resultExtendedLoc = result.getExtendedSpan();
        final int resultExtendedLocStart = resultExtendedLoc.getStart();
        final int resultExtendedLocStop = resultExtendedLoc.getStop();

        final List<SAMRecord> trimmedReads = new ArrayList<>(myReads.size());
        for( final SAMRecord read : myReads ) {
            final SAMRecord clippedRead = ReadClipper.hardClipToRegion(read, resultExtendedLocStart, resultExtendedLocStop);
            if( result.readOverlapsRegion(clippedRead) && clippedRead.getReadLength() > 0 ) {
                trimmedReads.add(clippedRead);
            }
        }
        result.clearReads();
        trimmedReads.sort(new SAMRecordCoordinateComparator());
        result.addAll(trimmedReads);
        return result;
    }

    /**
     * Returns true if read would overlap the extended extent of this region
     * @param read the read we want to test
     * @return true if read can be added to this region, false otherwise
     */
    public boolean readOverlapsRegion(final SAMRecord read) {
        final GenomeLoc readLoc = genomeLocParser.createGenomeLoc( read );
        return readLoc.overlapsP(extendedLoc);
    }

    /**
     * Add read to this active region
     *
     * Read must have alignment start >= than the last read currently in this active region.
     *
     * @throws IllegalArgumentException if read doesn't overlap the extended region of this active region
     *
     * @param read a non-null GATKSAMRecord
     */
    public void add( final SAMRecord read ) {
        if ( read == null ) throw new IllegalArgumentException("Read cannot be null");

        final GenomeLoc readLoc = genomeLocParser.createGenomeLoc( read );
        if ( ! readOverlapsRegion(read) )
            throw new IllegalArgumentException("Read location " + readLoc + " doesn't overlap with active region extended span " + extendedLoc);

        spanIncludingReads = spanIncludingReads.union( readLoc );

        if ( ! reads.isEmpty() ) {
            final SAMRecord lastRead = reads.get(size() - 1);
            if ( ! lastRead.getReferenceIndex().equals(read.getReferenceIndex()) )
                throw new IllegalArgumentException("Attempting to add a read to ActiveRegion not on the same contig as other reads: lastRead " + lastRead + " attempting to add " + read);

            if ( read.getAlignmentStart() < lastRead.getAlignmentStart() )
                throw new IllegalArgumentException("Attempting to add a read to ActiveRegion out of order w.r.t. other reads: lastRead " + lastRead + " at " + lastRead.getAlignmentStart() + " attempting to add " + read + " at " + read.getAlignmentStart());
        }

        reads.add( read );
    }

    /**
     * Get the number of reads currently in this active region
     * @return an integer >= 0
     */
    public int size() { return reads.size(); }

    /**
     * Clear all of the reads currently in this active region
     */
    public void clearReads() {
        spanIncludingReads = extendedLoc;
        reads.clear();
    }

    /**
     * Add all reads to this active region
     * @param reads a collection of reads to add to this active region
     */
    public void addAll(final Collection<SAMRecord> reads){
        if ( reads == null ) throw new IllegalArgumentException("reads cannot be null");
        for ( final SAMRecord read : reads ) {
            add(read);
        }
    }

    /**
     * Get the extension applied to this region
     *
     * The extension is >= 0 bp in size, and indicates how much padding this art walker wanted for its regions
     *
     * @return the size in bp of the region extension
     */
    public int getExtension() { return extension; }
}
