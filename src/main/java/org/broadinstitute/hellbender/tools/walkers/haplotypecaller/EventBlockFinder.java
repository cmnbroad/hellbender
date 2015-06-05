package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.HaplotypeGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.utils.collections.CountSet;

import java.util.*;

/**
 * Encapsulates the graph traversals needed to find event-blocks.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class EventBlockFinder {

    private final HaplotypeGraph graph;

    private final Map<Pair<MultiDeBruijnVertex,MultiDeBruijnVertex>,EventBlock> eventBlockCache;

    /**
     * Constructs a new engine.
     *
     * @param graph the base haplotype graph to iterate over.
     */
    public EventBlockFinder(final HaplotypeGraph graph) {
        if (graph == null) throw new NullPointerException();
        this.graph = graph;
        eventBlockCache = new HashMap<>(20);
    }

    /**
     * Create a new traversal object based on a read anchoring.
     * @param anchoring
     * @return never {@code null}.
     */
    public Traversal traversal(final ReadAnchoring anchoring) {
        if (anchoring == null) throw new NullPointerException();
        return new Traversal(anchoring);
    }


    public class Traversal implements Iterable<EventBlock> {

        private final ReadAnchoring anchoring;

        private EventBlock lastEventBlock;


        private Traversal(final ReadAnchoring anchoring) {
            this.anchoring = anchoring;
            lastEventBlock = findLastEventBlock(anchoring);
        }

        @Override
        public java.util.Iterator<EventBlock> iterator() {
            return lastEventBlock == null ? Collections.EMPTY_SET.iterator() : new Iterator();
        }

        private class Iterator implements java.util.Iterator<EventBlock> {

            private MultiDeBruijnVertex currentVertex;

            private Iterator() {
                currentVertex = anchoring.leftAnchorVertex;
            }

            @Override
            public boolean hasNext() {
                return currentVertex != null;
            }

            @Override
            public EventBlock next() {
                final EventBlock result;
                if (currentVertex == null)
                    throw new NoSuchElementException("going beyond last event block");
                else if (currentVertex == lastEventBlock.getSource()) {
                    result = lastEventBlock;
                    currentVertex = null;
                } else {
                    final EventBlock candidate = findEventBlock(anchoring,false,currentVertex,lastEventBlock.getSource());
                    if (candidate == null) {
                        result = findEventBlock(anchoring,false,currentVertex,anchoring.rightAnchorVertex);
                        currentVertex = null;
                    } else {
                        result = candidate;
                        currentVertex = candidate.getSink();
                    }
                }
                return result;
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        }
    }

    /**
     * Finds the last event block.
     * <p>
     * It can do it forward or backwards.
     * </p>
     *
     * @param anchoring target read anchoring information.
     * @return {@code null} if there is no event block, depending on {@code backwards} before or after current
     */
    private EventBlock findLastEventBlock(
            final ReadAnchoring anchoring) {
        return findEventBlock(anchoring,true,anchoring.leftAnchorVertex,anchoring.rightAnchorVertex);
    }

    /**
     * Finds an event block forward or backwards along the reference route.
     * @param anchoring the read anchoring information.
     * @param backwards true if the block should be constructed from right to left.
     * @param leftVertex the left vertex
     * @param rightVertex the right vertex
     * @return {@code null} if there is no such a event block between these coordinates.
     */
    private EventBlock findEventBlock(
            final ReadAnchoring anchoring, final boolean backwards,
            final MultiDeBruijnVertex leftVertex, final MultiDeBruijnVertex rightVertex) {

        MultiDeBruijnVertex currentVertex = backwards ? rightVertex : leftVertex;
        boolean foundEvent = false;
        final CountSet pathSizes = new CountSet(10); // typically more than enough.
        pathSizes.setTo(0);

        // Map between reference vertices where there is some expected open alternative path rejoining and the
        // predicted length of paths rejoining at that point counting from the beginning of the block.
        final Map<MultiDeBruijnVertex, CountSet> expectedAlternativePathRejoins = new HashMap<>(4);

        // Keeps record of possible left-clipping veritces; those that are located before any event path furcation
        // has been found. The value indicates the blockLength at the time we traverse that node.
        final Deque<Pair<MultiDeBruijnVertex, Integer>> possibleClippingPoints = new LinkedList<>();

        // We keep the distance from the beggining of the block (leftVertex).
        int blockLength = 0;
        while (currentVertex != null) {
            int openingDegree = backwards ? graph.outDegreeOf(currentVertex) : graph.inDegreeOf(currentVertex);
            if (openingDegree > 1) {
                final CountSet joiningPathLengths = expectedAlternativePathRejoins.remove(currentVertex);
                if (joiningPathLengths != null)
                    pathSizes.addAll(joiningPathLengths);
            }
            final boolean isValidBlockEnd = isValidBlockEnd(anchoring, currentVertex, expectedAlternativePathRejoins);
            if (foundEvent && isValidBlockEnd) // !gotcha we found a valid block end.
                break;
            else if (!foundEvent && isValidBlockEnd) // if no event has been found yet, still is a good clipping point.
                possibleClippingPoints.addLast(new MutablePair<>(currentVertex, blockLength));

            // We reached the end:
            if (currentVertex == (backwards ? leftVertex : rightVertex))
                break;

            // process next vertices, the next one on the reference and also possible start of alternative paths,
            // updates traversal structures accordingly.
            currentVertex = advanceOnReferencePath(anchoring, backwards, currentVertex, pathSizes, expectedAlternativePathRejoins);
            foundEvent |= expectedAlternativePathRejoins.size() > 0;
            pathSizes.incAll(1);
            blockLength++;
        }

        // we have not found an event, thus there is no block to report:
        if (!foundEvent)
            return null;

        // We try to clip off as much as we can from the beginning of the block before any event, but at least
        // leaving enough block length to meet the shortest path unless all paths have the same size (SNPs only)
        final int maxClipping = pathSizes.size() <= 1 ? blockLength : pathSizes.min();
        MultiDeBruijnVertex clippingEnd = backwards ? anchoring.rightAnchorVertex : anchoring.leftAnchorVertex;
        while (!possibleClippingPoints.isEmpty()) {
            final Pair<MultiDeBruijnVertex, Integer> candidate = possibleClippingPoints.removeLast();
            if (candidate.getRight() <= maxClipping) {
                clippingEnd = candidate.getLeft();
                break;
            }
        }

        return resolveEventBlock(backwards ? new MutablePair<>(currentVertex, clippingEnd) : new MutablePair<>(clippingEnd, currentVertex));
    }

    /**
     * Gets or constructs a event-block through the cache.
     * @param borders the source and sink vertex pair for the requested event block.
     * @return never {@code null}
     */
    private EventBlock resolveEventBlock(final Pair<MultiDeBruijnVertex,MultiDeBruijnVertex> borders) {
        EventBlock result = eventBlockCache.get(borders);
        if (result == null)
            eventBlockCache.put(borders,result = new EventBlock(graph, borders.getLeft(),borders.getRight()));
        return result;
    }

    /**
     * Move on vertex along the reference path checking for the presence of new opening alternative paths.
     *
     * @param anchoring anchoring information on the targeted read.
     * @param backwards whether we are extending the block backwards or forwards.
     * @param currentVertex the current vertex.
     * @param pathSizes current block path sizes.
     * @param expectedAlternativePathRejoins information about location of vertices along the reference path where open alternative paths will rejoin.
     * @return the next current-vertex, never {@code null} unless there is a bug.
     */
    private MultiDeBruijnVertex advanceOnReferencePath(final ReadAnchoring anchoring, final boolean backwards, final MultiDeBruijnVertex currentVertex, final CountSet pathSizes, final Map<MultiDeBruijnVertex, CountSet> expectedAlternativePathRejoins) {
        final Set<MultiSampleEdge> nextEdges = backwards ? graph.incomingEdgesOf(currentVertex) : graph.outgoingEdgesOf(currentVertex);
        MultiDeBruijnVertex nextReferenceVertex = null;
        for (final MultiSampleEdge e : nextEdges) {
            final MultiDeBruijnVertex nextVertex = backwards ? graph.getEdgeSource(e) : graph.getEdgeTarget(e);
            if (e.isRef())
                nextReferenceVertex = nextVertex;
            else {
                final CountSet pathSizesPlusOne = pathSizes.clone();
                pathSizesPlusOne.incAll(1);
                graph.calculateRejoins(nextVertex, expectedAlternativePathRejoins, anchoring.referenceWithinAnchorsMap.keySet(), pathSizesPlusOne, true, backwards);
            }
        }
        return nextReferenceVertex;
    }

    /**
     * Check whether the current vertex is a valid block end.
     *
     * @param anchoring reads anchoring information necessary to make the evaluation.
     * @param currentVertex target potential block end
     * @param expectedAlternativePathRejoins traversal states regarding open alternative paths.
     *
     * @return {@code true} iff so.
     */
    private boolean isValidBlockEnd(final ReadAnchoring anchoring, final MultiDeBruijnVertex currentVertex, final Map<MultiDeBruijnVertex, CountSet> expectedAlternativePathRejoins) {
        final boolean isUniqueKmer = anchoring.uniqueKmerOffsets.containsKey(currentVertex);
        final boolean isAnchorable = graph.getAnchorableVertices().contains(currentVertex) && isUniqueKmer && expectedAlternativePathRejoins.size() == 0;
        return isUniqueKmer && isAnchorable;
    }
}
