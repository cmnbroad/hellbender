package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.Path;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.HaplotypeGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;

import java.util.*;

/**
 * Collects information as to how a read maps into the haplotype graph that is needed to calculate its likelihood
 * using the graph-based approach.
 */
public class ReadAnchoring {

    private static final Logger logger = LogManager.getLogger(ReadAnchoring.class);

    /** Holds a reference to the read itself */
    protected final SAMRecord read;
    protected final Map<MultiDeBruijnVertex,Integer> uniqueKmerOffsets;

    /**
     * Kmer offset on the read of the left anchor
     * <p>
     * {@code -1} if there is no left anchor.
     * </p>
     */
    protected int leftAnchorIndex;

    /**
     * Vertex in the graph where the left anchor falls.
     * <p>
     *     {@code null} if there is no left anchor.
     * </p>
     */
    protected MultiDeBruijnVertex leftAnchorVertex;

    /**
     * Kmer offset on the read of the right anchor.
     *
     * <p>
     *     {@code -1} if there is no right anchor.
     * </p>
     */
    protected int rightAnchorIndex;

    /**
     * Vertex in the graph where the right anchor falls.
     *
     * <p>
     *     {@code null} if there is no right anchor.
     * </p>
     */
    protected MultiDeBruijnVertex rightAnchorVertex;

    /**
     * Kmer sequence mapping information for the read sequence.
     *
     * never {@code null}.
     */
    protected final KmerSequenceGraphMap<MultiDeBruijnVertex, MultiSampleEdge> graphMap;

    /**
     * Alignment of read kmers on the reference haplotype kmers.
     *
     * <p>
     *     There is one entry for each base in the read.
     *
     * </p>
     * <p>
     *     The i-th entry indicates what kmer in the reference haplotype correspond to the kmer on the read starting
     *     at is i-th base.
     * </p>
     *
     * <p>
     *     {@code -1} means that there is no match.
     * </p>
     *
     * <p>
     *     The last kmerSize - 1 entry of the array are {@code -1}
     * </p>
     */
    protected final int[] referenceAlignment;

    /**
     * Maps between reference path vertex that are found between anchors and the kmer offset they map uniquely to
     * on the read.
     */
    protected final Map<MultiDeBruijnVertex, Integer> referenceWithinAnchorsMap;

    /**
     * Creates the read's anchoring information for the haplotype-graph.
     *
     * @param read the targeted read.
     * @param haplotypeGraph the targeted graph.
     *
     * @throws NullPointerException if any argument is {@code null}.
     * @throws IllegalArgumentException if elements in {@code anchorableVertices} are not vertex in {@code haplotypeGraph}
     */
    public ReadAnchoring(final SAMRecord read, final HaplotypeGraph haplotypeGraph) {
        this.read = read;
        final byte[] readBases = read.getReadBases();
        final KmerSequence readKmers = new KmerSequence(read, haplotypeGraph.getKmerSize());
        graphMap = new KmerSequenceGraphMap<>(haplotypeGraph, readKmers);
        final Map<MultiDeBruijnVertex, Integer> vertexOffset = graphMap.vertexOffset();
        referenceAlignment = calculateUniqueKmerAlignment(0, readBases.length, haplotypeGraph.getReferenceRoute(), vertexOffset, haplotypeGraph.getKmerSize());
        leftAnchorIndex = -1;
        leftAnchorVertex = null;
        for (int i = 0; i < readBases.length - haplotypeGraph.getKmerSize() + 1; i++) {
            if (referenceAlignment[i] == -1) continue;
            final MultiDeBruijnVertex candidate = haplotypeGraph.findKmer(readKmers.get(i));
            if (candidate != null && haplotypeGraph.getAnchorableVertices().contains(candidate)) {
                leftAnchorIndex = i;
                leftAnchorVertex = candidate;
                break;
            }
        }
        rightAnchorIndex = leftAnchorIndex;
        rightAnchorVertex = leftAnchorVertex;
        if (leftAnchorIndex != -1) {
            for (int i = readBases.length - haplotypeGraph.getKmerSize(); i > leftAnchorIndex; i--) {
                if (referenceAlignment[i] == -1) continue;
                final MultiDeBruijnVertex candidate = haplotypeGraph.findKmer(readKmers.get(i));
                if (candidate != null && haplotypeGraph.getAnchorableVertices().contains(candidate)) {
                    rightAnchorIndex = i;
                    rightAnchorVertex = candidate;
                    break;
                }
            }
        }
        referenceWithinAnchorsMap = buildReferenceWithinBoundariesMap(read, haplotypeGraph,
                vertexOffset, leftAnchorVertex, rightAnchorVertex);
        uniqueKmerOffsets = buildReadUniqueKmerOffsets(haplotypeGraph);
    }

    /**
     * For a given read, returns the set of reference path vertices that falls between the two anchor vertices.
     * <p/>
     * <p>
     * The resulting map has as key the reference vertices between those two boundaries (inclusive) and
     * the value is the corresponding offset in the kmer.
     * </p>
     *
     * @param read                 the target read.
     * @param readVertexKmerOffset map between vertices and their kmer offset on the read.
     * @param leftAnchorVertex     left anchor vertex.
     * @param rightAnchorVertex    right anchor vertex.
     * @return never {@code null}, but empty if the anchors are out of order in the reference.
     */
    private Map<MultiDeBruijnVertex, Integer> buildReferenceWithinBoundariesMap(
            final SAMRecord read, final HaplotypeGraph haplotypeGraph,
            final Map<MultiDeBruijnVertex, Integer> readVertexKmerOffset,
            final MultiDeBruijnVertex leftAnchorVertex, final MultiDeBruijnVertex rightAnchorVertex) {
        if (leftAnchorVertex == null)
            return Collections.emptyMap();

        final Map<MultiDeBruijnVertex, Integer> result = new HashMap<>();
        MultiDeBruijnVertex nextVertex = leftAnchorVertex;

        int leftAnchorOffset = 0;
        while (nextVertex != null) {
            result.put(nextVertex, leftAnchorOffset++);
            if (nextVertex == rightAnchorVertex)
                break;
            nextVertex = haplotypeGraph.getNextReferenceVertex(nextVertex);
        }
        if (nextVertex == null) {
            logger.warn("unexpected event kmers out of order between read anchor kmers: " + read.getReadString()
                    + " Offending kmer offsets: " + readVertexKmerOffset.get(leftAnchorVertex) + " " + readVertexKmerOffset.get(rightAnchorVertex)
                    + " sequences: " +
                    read.getReadString().substring(readVertexKmerOffset.get(leftAnchorVertex), haplotypeGraph.getKmerSize() + readVertexKmerOffset.get(leftAnchorVertex)) +
                    " " + read.getReadString().substring(readVertexKmerOffset.get(rightAnchorVertex), haplotypeGraph.getKmerSize() + readVertexKmerOffset.get(rightAnchorVertex)) +
                    " Reference haplotype: " + haplotypeGraph.getReferenceHaplotype().getBaseString());
            return Collections.emptyMap();
        }
        return result;
    }

    /**
     * Builds a map between unique kmers in the reference path and their kmer offset in the read.
     *
     * @param haplotypeGraph the anchoring graph.
     *
     * @return never {@code null}.
     */
    private Map<MultiDeBruijnVertex,Integer> buildReadUniqueKmerOffsets(final HaplotypeGraph haplotypeGraph) {
        if (!hasValidAnchors())
            return Collections.emptyMap();
        final Map<MultiDeBruijnVertex, Integer> vertexOffset = graphMap.vertexOffset();
        final Set<MultiDeBruijnVertex> readUniqueKmerVertices = new HashSet<>(vertexOffset.size());
        readUniqueKmerVertices.add(leftAnchorVertex);
        readUniqueKmerVertices.add(rightAnchorVertex);
        for (int i = leftAnchorIndex + 1; i < rightAnchorIndex; i++) {
            if (referenceAlignment[i] != -1) {
                readUniqueKmerVertices.add(haplotypeGraph.findKmer(graphMap.sequence.get(i)));
            }
        }
        final Map<MultiDeBruijnVertex, Integer> validVertexOffset = new HashMap<>(graphMap.vertexOffset());
        validVertexOffset.keySet().retainAll(readUniqueKmerVertices);
        return validVertexOffset;
    }

    /**
     * Checks whether it has some anchoring kmer and these are valid, i.e. the left anchor is the same or preceedes the right anchor in the reference path.
     * @return {@code true} iff so.
     */
    public boolean hasValidAnchors() {
        return referenceWithinAnchorsMap.size() >= 1;
    }

    /**
     * Calculates an array indicating for each kmer in the read what is the offset of that kmer in a path.
     * <p/>
     * <p>
     * The result is of the same length as the read. Position ith indicates the offset of the read kmer that
     * start at that position in the input path. Non matching kmers have -1 instead.
     * </p>
     *
     * @param readStart            inclusive first position of the read to consider.
     * @param readEnd              exclusive position after last to be considered.
     * @param path                 the path to which to align against.
     * @param readUniqueKmerOffset map of vertices to the kmer offset with the read.
     * @return never {@code null}.
     */
    private int[] calculateUniqueKmerAlignment(final int readStart, final int readEnd, final Path<MultiDeBruijnVertex, MultiSampleEdge> path,
                                               final Map<MultiDeBruijnVertex, Integer> readUniqueKmerOffset, final int kmerSize) {

        final int[] result = new int[readEnd - readStart];
        Arrays.fill(result, -1);
        int i = 0;
        for (final MultiDeBruijnVertex v : path.getVertices()) {
            final Integer kmerReadOffset = readUniqueKmerOffset.get(v);
            if (kmerReadOffset != null) {
                final int kro = kmerReadOffset;
                if (kro >= readStart && kro < readEnd - kmerSize + 1) {
                    result[kro - readStart] = i;
                }
            }
            i++;
        }
        // Now we remove conflicting mappings:
        // A conflicting mapping is when to kmer mapping suggest that
        // the same read position maps to two different bases in the path.
        maskOutConflictingKmerAlignments(result,kmerSize);
        return result;
    }

    /**
     * Mark with -1 those kmer matches that result in read base mapping conflicts.
     *
     * @param result in/out changed in-situ.
     */
    private void maskOutConflictingKmerAlignments(final int[] result, final int kmerSize) {
        int i;
        int lastKmer = -1;
        int lastKmerPos = -1;
        for (i = 0; i < result.length; i++) {
            final int kmer = result[i];
            if (kmer == -1)
                continue;
            if (lastKmer == -1) {
                lastKmer = kmer;
                lastKmerPos = i;
            } else if (lastKmer + kmerSize - 1 >= kmer && (i - lastKmerPos) != (kmer - lastKmer)) { // kmer overlap. fixing by eliminating offending kmers alignments.
                int iSkip = result.length;  // iSkip will contain the next position minus 1 to visit in the next iteration of the enclosing loop.
                for (int j = i; j < result.length; j++)
                    if (result[j] != -1) {
                        if (lastKmer + kmerSize - 1 >= result[j])
                            result[j] = -1;
                        else {
                            iSkip = j;
                            break;
                        }
                    }
                // then backwards and do the same.
                int j = lastKmerPos;
                lastKmer = -1;
                lastKmerPos = -1;
                for (; j >= 0; j--)
                    if (result[j] != -1) {
                        if (result[j] + kmerSize - 1 >= kmer)
                            result[j] = -1;
                        else {
                            lastKmer = result[j];
                            lastKmerPos = j;
                            break;
                        }
                    }
                i = iSkip;
            } else {
                lastKmer = kmer;
                lastKmerPos = i;
            }
        }
    }

    /**
     * Checks whether it is anchored at all.
     *
     * @return {@code true} iff so.
     */
    public boolean isAnchoredSomewhere() {
        return hasValidAnchors();
        //return hasValidAnchors();
    }

    /**
     * Whether the read is anchored perfectly, there are no non-aligned bases.
     *
     * @return {@code true} iff so.
     */
    public boolean isPerfectAnchoring() {
        return hasValidAnchors() && leftAnchorIndex == 0 && rightAnchorIndex == read.getReadLength() - graphMap.kmerSize &&
               !leftAnchorVertex.hasAmbiguousSequence();
    }

}
