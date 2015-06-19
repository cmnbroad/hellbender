package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.KMerCounter;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.*;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SWPairwiseAlignment;
import org.broadinstitute.hellbender.utils.smithwaterman.SWParameterSet;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWaterman;
import org.jgrapht.EdgeFactory;

import javax.rmi.CORBA.Util;
import java.io.File;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ReadThreadingGraph extends BaseGraph<MultiDeBruijnVertex, MultiSampleEdge> implements KmerSearchableGraph<MultiDeBruijnVertex,MultiSampleEdge> {

    private static final Logger logger = LogManager.getLogger(ReadThreadingGraph.class);

    private static final String ANONYMOUS_SAMPLE = "XXX_UNNAMED_XXX";
    private static final boolean WRITE_GRAPH = false;
    private static final boolean DEBUG_NON_UNIQUE_CALC = false;

    private static final int MAX_CIGAR_COMPLEXITY = 3;
    private static final long serialVersionUID = 1l;
    private int maxMismatchesInDanglingHead = -1;

    protected boolean alreadyBuilt;


    private boolean startThreadingOnlyAtExistingVertex = false;

    /** for debugging info printing */
    private static int counter = 0;

    /**
     * Sequences added for read threading before we've actually built the graph
     */
    private final Map<String, List<SequenceForKmers>> pending = new LinkedHashMap<>();

    /**
     * A set of non-unique kmers that cannot be used as merge points in the graph
     */
    protected Set<Kmer> nonUniqueKmers;

    /**
     * A map from kmers -> their corresponding vertex in the graph
     */
    protected Map<Kmer, MultiDeBruijnVertex> uniqueKmers = new LinkedHashMap<>();

    /**
     *
     */
    private final boolean debugGraphTransformations;
    final byte minBaseQualityToUseInAssembly;

    protected boolean increaseCountsBackwards = true;
    protected boolean increaseCountsThroughBranches = false; // this may increase the branches without bounds

    // --------------------------------------------------------------------------------
    // state variables, initialized in resetToInitialState()
    // --------------------------------------------------------------------------------
    private Kmer refSource;

    /**
     * Constructs an empty read-threading-grpah provided the kmerSize.
     * @param kmerSize 1 or greater.
     *
     * @throws IllegalArgumentException if (@code kmerSize) < 1.
     */
    public ReadThreadingGraph(final int kmerSize) {
        this(kmerSize, false, (byte)6, 1);
    }

    @VisibleForTesting
    final void setMaxMismatchesInDanglingHead(final int maxMismatchesInDanglingHead) {
        this.maxMismatchesInDanglingHead = maxMismatchesInDanglingHead;
    }

    /**
     * Return the collection of outgoing vertices that expand this vertex with a particular base.
     *
     * @param v original vertex.
     * @param b expanding base.
     * @return never null, but perhaps an empty set. You cannot assume that you can modify the result.
     */
    protected Set<MultiDeBruijnVertex> getNextVertices(final MultiDeBruijnVertex v, final byte b) {
        Utils.nonNull(v, "the input vertex cannot be null");
        if (!vertexSet().contains(v)) throw new IllegalArgumentException("the vertex must be present in the graph");
        final List<MultiDeBruijnVertex> result = new LinkedList<>();
        for (final MultiDeBruijnVertex w : outgoingVerticesOf(v)) {
            if (w.getSuffix() == b) {
                result.add(w);
            }
        }
        switch (result.size()) {
            case 0: return Collections.emptySet();
            case 1: return Collections.singleton(result.get(0));
            default:
                    return new HashSet<>(result);
        }
    }

    /**
     * Create a new ReadThreadingAssembler using kmerSize for matching
     * @param kmerSize must be >= 1
     */
    protected ReadThreadingGraph(final int kmerSize, final boolean debugGraphTransformations, final byte minBaseQualityToUseInAssembly, final int numPruningSamples) {
        super(kmerSize, new MyEdgeFactory(numPruningSamples));

        if ( kmerSize < 1 ) throw new IllegalArgumentException("bad minkKmerSize " + kmerSize);
        this.debugGraphTransformations = debugGraphTransformations;
        this.minBaseQualityToUseInAssembly = minBaseQualityToUseInAssembly;

        resetToInitialState();
    }

    /**
     * Reset this assembler to its initial state, so we can create another assembly with a different set of reads
     */
    private void resetToInitialState() {
        pending.clear();
        nonUniqueKmers = null;
        uniqueKmers.clear();
        refSource = null;
        alreadyBuilt = false;
    }

    /**
     * Add the all bases in sequence to the graph
     * @param sequence a non-null sequence
     * @param isRef is this the reference sequence?
     */
    protected void addSequence(final byte[] sequence, final boolean isRef) {
        addSequence("anonymous", sequence, isRef);
    }

    /**
     * Add all bases in sequence to this graph
     *
     * @see #addSequence(String, String, byte[], int, int, int, boolean) for full information
     */
    public void addSequence(final String seqName, final byte[] sequence, final boolean isRef) {
        addSequence(seqName, sequence, 1, isRef);
    }

    /**
     * Add all bases in sequence to this graph
     *
     * @see #addSequence(String, String, byte[], int, int, int, boolean) for full information
     */
    public void addSequence(final String seqName, final byte[] sequence, final int count, final boolean isRef) {
        addSequence(seqName, ANONYMOUS_SAMPLE, sequence, 0, sequence.length, count, isRef);
    }

    /**
     * Add bases in sequence to this graph
     *
     * @param seqName a useful seqName for this read, for debugging purposes
     * @param sequence non-null sequence of bases
     * @param start the first base offset in sequence that we should use for constructing the graph using this sequence, inclusive
     * @param stop the last base offset in sequence that we should use for constructing the graph using this sequence, exclusive
     * @param count the representative count of this sequence (to use as the weight)
     * @param isRef is this the reference sequence.
     */
    public void addSequence(final String seqName, final String sampleName, final byte[] sequence, final int start, final int stop, final int count, final boolean isRef) {
        // note that argument testing is taken care of in SequenceForKmers
        if ( alreadyBuilt ) throw new IllegalStateException("Graph already built");

        // get the list of sequences for this sample
        List<SequenceForKmers> sampleSequences = pending.get(sampleName);
        if ( sampleSequences == null ) { // need to create
            sampleSequences = new LinkedList<>();
            pending.put(sampleName, sampleSequences);
        }

        // add the new sequence to the list of sequences for sample
        sampleSequences.add(new SequenceForKmers(seqName, sequence, start, stop, count, isRef));
    }

    /**
     * Thread sequence seqForKmers through the current graph, updating the graph as appropriate
     * @param seqForKmers a non-null sequence
     */
    private void threadSequence(final SequenceForKmers seqForKmers) {
        final int uniqueStartPos = findStart(seqForKmers);
        if ( uniqueStartPos == -1 )
            return;

        final MultiDeBruijnVertex startingVertex = getOrCreateKmerVertex(seqForKmers.sequence, uniqueStartPos);

        // increase the counts of all edges incoming into the starting vertex supported by going back in sequence
        if ( increaseCountsBackwards )
            increaseCountsInMatchedKmers(seqForKmers, startingVertex, startingVertex.getSequence(), kmerSize - 2);

        if ( debugGraphTransformations ) startingVertex.addRead(seqForKmers.name);

        // keep track of information about the reference source
        if ( seqForKmers.isRef ) {
            if ( refSource != null ) throw new IllegalStateException("Found two refSources! prev: " + refSource + ", new: " + startingVertex);
            refSource = new Kmer(seqForKmers.sequence, seqForKmers.start, kmerSize);
        }

        // loop over all of the bases in sequence, extending the graph by one base at each point, as appropriate
        MultiDeBruijnVertex vertex = startingVertex;
        for ( int i = uniqueStartPos + 1; i <= seqForKmers.stop - kmerSize; i++ ) {
            vertex = extendChainByOne(vertex, seqForKmers.sequence, i, seqForKmers.count, seqForKmers.isRef);
            if ( debugGraphTransformations ) vertex.addRead(seqForKmers.name);
        }
    }

    /**
     * Find vertex and its position in seqForKmers where we should start assembling seqForKmers
     *
     * @param seqForKmers the sequence we want to thread into the graph
     * @return the position of the starting vertex in seqForKmer, or -1 if it cannot find one
     */
    protected int findStart(final SequenceForKmers seqForKmers) {
        if ( seqForKmers.isRef )
            return 0;

        for ( int i = seqForKmers.start; i < seqForKmers.stop - kmerSize; i++ ) {
            final Kmer kmer1 = new Kmer(seqForKmers.sequence, i, kmerSize);
            if ( isThreadingStart(kmer1) )
                return i;
        }

        return -1;
    }

    /**
     * Checks whether a kmer can be the threading start based on the current threading start location policy.
     *
     * @see #setThreadingStartOnlyAtExistingVertex(boolean)
     * @see #getThreadingStartOnlyAtExistingVertex()
     *
     * @param kmer the query kmer.
     * @return {@code true} if we can start thread the sequence at this kmer, {@code false} otherwise.
     */
    protected boolean isThreadingStart(final Kmer kmer) {
        Utils.nonNull(kmer);
        return startThreadingOnlyAtExistingVertex ? uniqueKmers.containsKey(kmer) : !nonUniqueKmers.contains(kmer);
    }

    /**
     * Changes the threading start location policy.
     *
     * @param value  {@code true} if threading will start only at existing vertices in the graph, {@code false} if
     *  it can start at any unique kmer.
     */
    public void setThreadingStartOnlyAtExistingVertex(final boolean value) {
        startThreadingOnlyAtExistingVertex = value;
    }

    /**
     * Indicates the threading start location policy.
     *
     * @return {@code true} if threading will start only at existing vertices in the graph, {@code false} if
     *  it can start at any unique kmer.
     */
    public boolean getThreadingStartOnlyAtExistingVertex() {
        return startThreadingOnlyAtExistingVertex;
    }

    /**
     * Build the read threaded assembly graph if it hasn't already been constructed from the sequences that have
     * been added to the graph.
     */
    public void buildGraphIfNecessary() {
        if ( alreadyBuilt ) return;

        // determine the kmer size we'll use, and capture the set of nonUniques for that kmer size
        final NonUniqueResult result = determineKmerSizeAndNonUniques(kmerSize, kmerSize);
        nonUniqueKmers = result.nonUniques;

        if ( DEBUG_NON_UNIQUE_CALC ) {
            logger.info("using " + kmerSize + " kmer size for this assembly with the following non-uniques");
        }

        // go through the pending sequences, and add them to the graph
        for ( final List<SequenceForKmers> sequencesForSample : pending.values() ) {
            for ( final SequenceForKmers sequenceForKmers : sequencesForSample ) {
                threadSequence(sequenceForKmers);
                if ( WRITE_GRAPH ) printGraph(new File("threading." + counter++ + "." + sequenceForKmers.name.replace(" ", "_") + ".dot"), 0);
            }

            // flush the single sample edge values from the graph
            for ( final MultiSampleEdge e : edgeSet() ) e.flushSingleSampleMultiplicity();
        }

        // clear
        pending.clear();
        alreadyBuilt = true;
        for (final MultiDeBruijnVertex v : uniqueKmers.values())
            v.setAdditionalInfo(v.additionalInfo() + "+");
    }


    @Override
    public boolean removeVertex(MultiDeBruijnVertex V) {
        final boolean result = super.removeVertex(V);
        if (result) {
            final byte[] sequence = V.getSequence();
            final Kmer kmer = new Kmer(sequence);
            uniqueKmers.remove(kmer);
        }
        return result;
    }


    public void removeSingletonOrphanVertices() {
        // Run through the graph and clean up singular orphaned nodes
        final List<MultiDeBruijnVertex> verticesToRemove = new LinkedList<>();
        for( final MultiDeBruijnVertex v : vertexSet() ) {
            if( inDegreeOf(v) == 0 && outDegreeOf(v) == 0 ) {
                verticesToRemove.add(v);
            }
        }
        this.removeVertex(null);
        removeAllVertices(verticesToRemove);
    }

    /**
     * Does the graph not have enough complexity?  We define low complexity as a situation where the number
     * of non-unique kmers is more than 20% of the total number of kmers.
     *
     * @return true if the graph has low complexity, false otherwise
     */
    public boolean isLowComplexity() {
        return nonUniqueKmers.size() * 4 > uniqueKmers.size();
    }

    @Override
    public ReadThreadingGraph clone() {
        return (ReadThreadingGraph) super.clone();
    }

    /**
     * Edge factory that encapsulates the numPruningSamples assembly parameter
     */
    private static final class MyEdgeFactory implements EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> {
        final int numPruningSamples;

        public MyEdgeFactory(int numPruningSamples) {
            this.numPruningSamples = numPruningSamples;
        }

        @Override
        public MultiSampleEdge createEdge(final MultiDeBruijnVertex sourceVertex, final MultiDeBruijnVertex targetVertex) {
            return new MultiSampleEdge(false, 1, numPruningSamples);
        }

        public MultiSampleEdge createEdge(final boolean isRef, final int multiplicity) {
            return new MultiSampleEdge(isRef, multiplicity, numPruningSamples);
        }
    }

    /**
     * Class to keep track of the important dangling chain merging data
     */
    static final class DanglingChainMergeHelper {
        final List<MultiDeBruijnVertex> danglingPath, referencePath;
        final byte[] danglingPathString, referencePathString;
        final Cigar cigar;

        public DanglingChainMergeHelper(final List<MultiDeBruijnVertex> danglingPath,
                                        final List<MultiDeBruijnVertex> referencePath,
                                        final byte[] danglingPathString,
                                        final byte[] referencePathString,
                                        final Cigar cigar) {
            this.danglingPath = danglingPath;
            this.referencePath = referencePath;
            this.danglingPathString = danglingPathString;
            this.referencePathString = referencePathString;
            this.cigar = cigar;
        }
    }

    /** structure that keeps track of the non-unique kmers for a given kmer size */
    private static final class NonUniqueResult {
        final Set<Kmer> nonUniques;
        final int kmerSize;

        private NonUniqueResult(Set<Kmer> nonUniques, int kmerSize) {
            this.nonUniques = nonUniques;
            this.kmerSize = kmerSize;
        }
    }

    /**
     * Keeps track of the information needed to add a sequence to the read threading assembly graph
     */
    static final class SequenceForKmers {
        final String name;
        final byte[] sequence;
        final int start, stop;
        final int count;
        final boolean isRef;

        /**
         * Create a new sequence for creating kmers
         */
        SequenceForKmers(final String name, final byte[] sequence, final int start, final int stop, final int count, final boolean ref) {
            Utils.nonNull(sequence, "Sequence is null ");
            if ( start < 0 ) throw new IllegalArgumentException("Invalid start " + start);
            if ( stop < start ) throw new IllegalArgumentException("Invalid stop " + stop);
            if ( count < 1 ) throw new IllegalArgumentException("Invalid count " + count);

            this.name = name;
            this.sequence = sequence;
            this.start = start;
            this.stop = stop;
            this.count = count;
            this.isRef = ref;
        }
    }

    /**
     * Try to recover dangling tails
     *
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     */
    public void recoverDanglingTails(final int pruneFactor, final int minDanglingBranchLength) {
        if ( ! alreadyBuilt )  throw new IllegalStateException("recoverDanglingTails requires the graph be already built");

        int attempted = 0;
        int nRecovered = 0;
        for ( final MultiDeBruijnVertex v : vertexSet() ) {
            if ( outDegreeOf(v) == 0 && ! isRefSink(v) ) {
                attempted++;
                nRecovered += recoverDanglingTail(v, pruneFactor, minDanglingBranchLength);
            }
        }

        logger.debug("Recovered " + nRecovered + " of " + attempted + " dangling tails");
    }


    /**
     * Try to recover dangling heads
     *
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     */
    public void recoverDanglingHeads(final int pruneFactor, final int minDanglingBranchLength) {
        if ( ! alreadyBuilt )  throw new IllegalStateException("recoverDanglingHeads requires the graph be already built");

        // we need to build a list of dangling heads because that process can modify the graph (and otherwise generate
        // a ConcurrentModificationException if we do it while iterating over the vertexes)
        final List<MultiDeBruijnVertex> danglingHeads = new ArrayList<>();

        int attempted = 0;
        int nRecovered = 0;
        for ( final MultiDeBruijnVertex v : vertexSet() ) {
            if ( inDegreeOf(v) == 0 && ! isRefSource(v) )
                danglingHeads.add(v);
        }

        // now we can try to recover the dangling heads
        for ( final MultiDeBruijnVertex v : danglingHeads ) {
            attempted++;
            nRecovered += recoverDanglingHead(v, pruneFactor, minDanglingBranchLength);
        }

        logger.debug("Recovered " + nRecovered + " of " + attempted + " dangling heads");
    }

    /**
     * Attempt to attach vertex with out-degree == 0 to the graph
     *
     * @param vertex the vertex to recover
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @return 1 if we successfully recovered the vertex and 0 otherwise
     */
    private int recoverDanglingTail(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {
        if ( outDegreeOf(vertex) != 0 ) throw new IllegalStateException("Attempting to recover a dangling tail for " + vertex + " but it has out-degree > 0");

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        final DanglingChainMergeHelper danglingTailMergeResult = generateCigarAgainstDownwardsReferencePath(vertex, pruneFactor, minDanglingBranchLength);

        // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
        if ( danglingTailMergeResult == null || ! cigarIsOkayToMerge(danglingTailMergeResult.cigar, false, true) )
            return 0;

        // merge
        return mergeDanglingTail(danglingTailMergeResult);
    }

    /**
     * Attempt to attach vertex with in-degree == 0, or a vertex on its path, to the graph
     *
     * @param vertex the vertex to recover
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @return 1 if we successfully recovered a vertex and 0 otherwise
     */
    private int recoverDanglingHead(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {
        if ( inDegreeOf(vertex) != 0 ) throw new IllegalStateException("Attempting to recover a dangling head for " + vertex + " but it has in-degree > 0");

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        final DanglingChainMergeHelper danglingHeadMergeResult = generateCigarAgainstUpwardsReferencePath(vertex, pruneFactor, minDanglingBranchLength);

        // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
        if ( danglingHeadMergeResult == null || ! cigarIsOkayToMerge(danglingHeadMergeResult.cigar, true, false) )
            return 0;

        // merge
        return mergeDanglingHead(danglingHeadMergeResult);
    }

    /**
     * Determine whether the provided cigar is okay to merge into the reference path
     *
     * @param cigar    the cigar to analyze
     * @param requireFirstElementM if true, require that the first cigar element be an M operator in order for it to be okay
     * @param requireLastElementM  if true, require that the last cigar element be an M operator in order for it to be okay
     * @return true if it's okay to merge, false otherwise
     */
    @VisibleForTesting
    final boolean cigarIsOkayToMerge(final Cigar cigar, final boolean requireFirstElementM, final boolean requireLastElementM) {

        final List<CigarElement> elements = cigar.getCigarElements();
        final int numElements = elements.size();

        // don't allow more than a couple of different ops
        if ( numElements == 0 || numElements > MAX_CIGAR_COMPLEXITY )
            return false;

        // the first element must be an M
        if ( requireFirstElementM && elements.get(0).getOperator() != CigarOperator.M )
            return false;

        // the last element must be an M
        if ( requireLastElementM && elements.get(numElements - 1).getOperator() != CigarOperator.M )
            return false;

        // note that there are checks for too many mismatches in the dangling branch later in the process

        return true;
    }

    /**
     * Actually merge the dangling tail if possible
     *
     * @param danglingTailMergeResult   the result from generating a Cigar for the dangling tail against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    @VisibleForTesting
    final int mergeDanglingTail(final DanglingChainMergeHelper danglingTailMergeResult) {

        final List<CigarElement> elements = danglingTailMergeResult.cigar.getCigarElements();
        final CigarElement lastElement = elements.get(elements.size() - 1);
        if ( lastElement.getOperator() != CigarOperator.M )
            throw new IllegalArgumentException("The last Cigar element must be an M");

        final int lastRefIndex = danglingTailMergeResult.cigar.getReferenceLength() - 1;
        final int matchingSuffix = Math.min(GraphUtils.longestSuffixMatch(danglingTailMergeResult.referencePathString, danglingTailMergeResult.danglingPathString, lastRefIndex), lastElement.getLength());
        if ( matchingSuffix == 0 )
            return 0;

        final int altIndexToMerge = Math.max(danglingTailMergeResult.cigar.getReadLength() - matchingSuffix - 1, 0);

        // there is an important edge condition that we need to handle here: Smith-Waterman correctly calculates that there is a
        // deletion, that deletion is left-aligned such that the LCA node is part of that deletion, and the rest of the dangling
        // tail is a perfect match to the suffix of the reference path.  In this case we need to push the reference index to merge
        // down one position so that we don't incorrectly cut a base off of the deletion.
        final boolean firstElementIsDeletion = elements.get(0).getOperator() == CigarOperator.D;
        final boolean mustHandleLeadingDeletionCase =  firstElementIsDeletion && (elements.get(0).getLength() + matchingSuffix == lastRefIndex + 1);
        final int refIndexToMerge = lastRefIndex - matchingSuffix + 1 + (mustHandleLeadingDeletionCase ? 1 : 0);

        // another edge condition occurs here: if Smith-Waterman places the whole tail into an insertion then it will try to
        // merge back to the LCA, which results in a cycle in the graph.  So we do not want to merge in such a case.
        if ( refIndexToMerge == 0 )
            return 0;

        // it's safe to merge now
        addEdge(danglingTailMergeResult.danglingPath.get(altIndexToMerge), danglingTailMergeResult.referencePath.get(refIndexToMerge), ((MyEdgeFactory)getEdgeFactory()).createEdge(false, 1));

        return 1;
    }

    /**
     * Actually merge the dangling head if possible
     *
     * @param danglingHeadMergeResult   the result from generating a Cigar for the dangling head against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    @VisibleForTesting
    final int mergeDanglingHead(final DanglingChainMergeHelper danglingHeadMergeResult) {

        final List<CigarElement> elements = danglingHeadMergeResult.cigar.getCigarElements();
        final CigarElement firstElement = elements.get(0);
        if ( firstElement.getOperator() != CigarOperator.M )
            throw new IllegalArgumentException("The first Cigar element must be an M");

        final int indexesToMerge = bestPrefixMatch(danglingHeadMergeResult.referencePathString, danglingHeadMergeResult.danglingPathString, firstElement.getLength());
        if ( indexesToMerge <= 0 )
            return 0;

        // we can't push back the reference path
        if ( indexesToMerge >= danglingHeadMergeResult.referencePath.size() - 1 )
            return 0;

        // but we can manipulate the dangling path if we need to
        if ( indexesToMerge >= danglingHeadMergeResult.danglingPath.size() &&
                ! extendDanglingPathAgainstReference(danglingHeadMergeResult, indexesToMerge - danglingHeadMergeResult.danglingPath.size() + 2) )
            return 0;

        addEdge(danglingHeadMergeResult.referencePath.get(indexesToMerge+1), danglingHeadMergeResult.danglingPath.get(indexesToMerge), ((MyEdgeFactory)getEdgeFactory()).createEdge(false, 1));

        return 1;
    }

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the sink) and the reference path.
     *
     * @param vertex   the sink of the dangling chain
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    @VisibleForTesting
    final DanglingChainMergeHelper generateCigarAgainstDownwardsReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {
        final int minTailPathLength = Math.max(1, minDanglingBranchLength); // while heads can be 0, tails absolutely cannot

        // find the lowest common ancestor path between this vertex and the diverging master path if available
        final List<MultiDeBruijnVertex> altPath = findPathUpwardsToLowestCommonAncestor(vertex, pruneFactor);
        if ( altPath == null || isRefSource(altPath.get(0)) || altPath.size() < minTailPathLength + 1 ) // add 1 to include the LCA
            return null;

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePath(altPath.get(0), TraversalDirection.downwards, Arrays.asList(incomingEdgeOf(altPath.get(1))));

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath, false);
        final byte[] altBases = getBasesForPath(altPath, false);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SmithWaterman alignment = new SWPairwiseAlignment(refBases, altBases, SWParameterSet.STANDARD_NGS, SWPairwiseAlignment.OVERHANG_STRATEGY.LEADING_INDEL);
        return new DanglingChainMergeHelper(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the source) and the reference path.
     *
     * @param vertex   the source of the dangling head
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    @VisibleForTesting
    final DanglingChainMergeHelper generateCigarAgainstUpwardsReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {

        // find the highest common descendant path between vertex and the reference source if available
        final List<MultiDeBruijnVertex> altPath = findPathDownwardsToHighestCommonDescendantOfReference(vertex, pruneFactor);
        if ( altPath == null || isRefSink(altPath.get(0)) || altPath.size() < minDanglingBranchLength + 1 ) // add 1 to include the LCA
            return null;

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePath(altPath.get(0), TraversalDirection.upwards, Collections.<MultiSampleEdge>emptyList());

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath, true);
        final byte[] altBases = getBasesForPath(altPath, true);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SmithWaterman alignment = new SWPairwiseAlignment(refBases, altBases, SWParameterSet.STANDARD_NGS, SWPairwiseAlignment.OVERHANG_STRATEGY.LEADING_INDEL);
        return new DanglingChainMergeHelper(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }

    /**
     * Finds the path upwards in the graph from this vertex to the first diverging node, including that (lowest common ancestor) vertex.
     * Note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex   the original vertex
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return the path if it can be determined or null if this vertex either doesn't merge onto another path or
     *  has an ancestor with multiple incoming edges before hitting the reference path
     */
    private List<MultiDeBruijnVertex> findPathUpwardsToLowestCommonAncestor(final MultiDeBruijnVertex vertex, final int pruneFactor) {
        final LinkedList<MultiDeBruijnVertex> path = new LinkedList<>();

        MultiDeBruijnVertex v = vertex;
        while ( inDegreeOf(v) == 1 && outDegreeOf(v) < 2 ) {
            final MultiSampleEdge edge = incomingEdgeOf(v);
            // if it has too low a weight, don't use it (or previous vertexes) for the path
            if ( edge.getPruningMultiplicity() < pruneFactor )
                path.clear();
                // otherwise it is safe to use
            else
                path.addFirst(v);
            v = getEdgeSource(edge);
        }
        path.addFirst(v);

        return outDegreeOf(v) > 1 ? path : null;
    }

    /**
     * Finds the path downwards in the graph from this vertex to the reference sequence, including the highest common descendant vertex.
     * However note that the path is reversed so that this vertex ends up at the end of the path.
     * Also note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex   the original vertex
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return the path if it can be determined or null if this vertex either doesn't merge onto the reference path or
     *  has a descendant with multiple outgoing edges before hitting the reference path
     */
    private List<MultiDeBruijnVertex> findPathDownwardsToHighestCommonDescendantOfReference(final MultiDeBruijnVertex vertex, final int pruneFactor) {
        final LinkedList<MultiDeBruijnVertex> path = new LinkedList<>();

        MultiDeBruijnVertex v = vertex;
        while ( ! isReferenceNode(v) && outDegreeOf(v) == 1 ) {
            final MultiSampleEdge edge = outgoingEdgeOf(v);
            // if it has too low a weight, don't use it (or previous vertexes) for the path
            if ( edge.getPruningMultiplicity() < pruneFactor )
                path.clear();
                // otherwise it is safe to use
            else
                path.addFirst(v);
            v = getEdgeTarget(edge);
        }
        path.addFirst(v);

        return isReferenceNode(v) ? path : null;
    }

    private enum TraversalDirection {
        downwards,
        upwards
    }

    /**
     * Finds the path in the graph from this vertex to the reference sink, including this vertex
     *
     * @param start   the reference vertex to start from
     * @param direction describes which direction to move in the graph (i.e. down to the reference sink or up to the source)
     * @param blacklistedEdges edges to ignore in the traversal down; useful to exclude the non-reference dangling paths
     * @return the path (non-null, non-empty)
     */
    private List<MultiDeBruijnVertex> getReferencePath(final MultiDeBruijnVertex start,
                                                       final TraversalDirection direction,
                                                       final Collection<MultiSampleEdge> blacklistedEdges) {

        final List<MultiDeBruijnVertex> path = new ArrayList<>();

        MultiDeBruijnVertex v = start;
        while ( v != null ) {
            path.add(v);
            v = (direction == TraversalDirection.downwards ? getNextReferenceVertex(v, true, blacklistedEdges) : getPrevReferenceVertex(v));
        }

        return path;
    }

    /**
     * The base sequence for the given path.
     *
     * @param path the list of vertexes that make up the path
     * @param expandSource if true and if we encounter a source node, then expand (and reverse) the character sequence for that node
     * @return  non-null sequence of bases corresponding to the given path
     */
    @VisibleForTesting
    final byte[] getBasesForPath(final List<MultiDeBruijnVertex> path, final boolean expandSource) {
        Utils.nonNull(path, "Path cannot be null");

        final StringBuilder sb = new StringBuilder();
        for ( final MultiDeBruijnVertex v : path ) {
            if ( expandSource && isSource(v) ) {
                final String seq = v.getSequenceString();
                sb.append(new StringBuilder(seq).reverse().toString());
            } else {
                sb.append((char)v.getSuffix());
            }
        }

        return sb.toString().getBytes();
    }

    /**
     * Finds the index of the best extent of the prefix match between the provided paths, for dangling head merging.
     * Assumes that path1.length >= maxIndex and path2.length >= maxIndex.
     *
     * @param path1  the first path
     * @param path2  the second path
     * @param maxIndex the maximum index to traverse (not inclusive)
     * @return the index of the ideal prefix match or -1 if it cannot find one, must be less than maxIndex
     */
    private int bestPrefixMatch(final byte[] path1, final byte[] path2, final int maxIndex) {
        final int maxMismatches = getMaxMismatches(maxIndex);
        int mismatches = 0;
        int index = 0;
        int lastGoodIndex = -1;
        while ( index < maxIndex ) {
            if ( path1[index] != path2[index] ) {
                if ( ++mismatches > maxMismatches )
                    return -1;
                lastGoodIndex = index;
            }
            index++;
        }
        // if we got here then we hit the max index
        return lastGoodIndex;
    }

    /**
     * Determine the maximum number of mismatches permitted on the branch.
     * Unless it's preset (e.g. by unit tests) it should be the length of the branch divided by the kmer size.
     *
     * @param lengthOfDanglingBranch  the length of the branch itself
     * @return positive integer
     */
    private int getMaxMismatches(final int lengthOfDanglingBranch) {
        return maxMismatchesInDanglingHead > 0 ? maxMismatchesInDanglingHead : Math.max(1, (lengthOfDanglingBranch / kmerSize));
    }

    private boolean extendDanglingPathAgainstReference(final DanglingChainMergeHelper danglingHeadMergeResult, final int numNodesToExtend) {

        final int indexOfLastDanglingNode = danglingHeadMergeResult.danglingPath.size() - 1;
        final int indexOfRefNodeToUse = indexOfLastDanglingNode + numNodesToExtend;
        if ( indexOfRefNodeToUse >= danglingHeadMergeResult.referencePath.size() )
            return false;

        final MultiDeBruijnVertex danglingSource = danglingHeadMergeResult.danglingPath.remove(indexOfLastDanglingNode);
        final StringBuilder sb = new StringBuilder();
        final byte[] refSourceSequence = danglingHeadMergeResult.referencePath.get(indexOfRefNodeToUse).getSequence();
        for ( int i = 0; i < numNodesToExtend; i++ )
            sb.append((char)refSourceSequence[i]);
        sb.append(danglingSource.getSequenceString());
        final byte[] sequenceToExtend = sb.toString().getBytes();

        // clean up the source and edge
        final MultiSampleEdge sourceEdge = outgoingEdgeOf(danglingSource);
        MultiDeBruijnVertex prevV = getEdgeTarget(sourceEdge);
        removeEdge(danglingSource, prevV);

        // extend the path
        for ( int i = numNodesToExtend; i > 0; i-- ) {
            final MultiDeBruijnVertex newV = new MultiDeBruijnVertex(Arrays.copyOfRange(sequenceToExtend, i, i + kmerSize));
            addVertex(newV);
            final MultiSampleEdge newE = addEdge(newV, prevV);
            newE.setMultiplicity(sourceEdge.getMultiplicity());
            danglingHeadMergeResult.danglingPath.add(newV);
            prevV = newV;
        }

        return true;
    }

    /**
     * Compute the smallest kmer size >= minKmerSize and <= maxKmerSize that has no non-unique kmers
     * among all sequences added to the current graph.  Will always return a result for maxKmerSize if
     * all smaller kmers had non-unique kmers.
     *
     * @param minKmerSize the minimum kmer size to consider when constructing the graph
     * @param maxKmerSize the maximum kmer size to consider
     * @return a non-null NonUniqueResult
     */
    protected NonUniqueResult determineKmerSizeAndNonUniques(final int minKmerSize, final int maxKmerSize) {
        final Collection<SequenceForKmers> withNonUniques = getAllPendingSequences();
        final Set<Kmer> nonUniqueKmers = new HashSet<>();

        // go through the sequences and determine which kmers aren't unique within each read
        int kmerSize = minKmerSize;
        for ( ; kmerSize <= maxKmerSize; kmerSize++) {
            // clear out set of non-unique kmers
            nonUniqueKmers.clear();

            // loop over all sequences that have non-unique kmers in them from the previous iterator
            final Iterator<SequenceForKmers> it = withNonUniques.iterator();
            while ( it.hasNext() ) {
                final SequenceForKmers sequenceForKmers = it.next();

                // determine the non-unique kmers for this sequence
                final Collection<Kmer> nonUniquesFromSeq = determineNonUniqueKmers(sequenceForKmers, kmerSize);
                if ( nonUniquesFromSeq.isEmpty() ) {
                    // remove this sequence from future consideration
                    it.remove();
                } else {
                    // keep track of the non-uniques for this kmerSize, and keep it in the list of sequences that have non-uniques
                    nonUniqueKmers.addAll(nonUniquesFromSeq);
                }
            }

            if ( nonUniqueKmers.isEmpty() )
                // this kmerSize produces no non-unique sequences, so go ahead and use it for our assembly
                break;
        }

        // necessary because the loop breaks with kmerSize = max + 1
        return new NonUniqueResult(nonUniqueKmers, Math.min(kmerSize, maxKmerSize));
    }

    /**
     * Get the collection of all sequences for kmers across all samples in no particular order
     * @return non-null Collection
     */
    private Collection<SequenceForKmers> getAllPendingSequences() {
        final LinkedList<SequenceForKmers> result = new LinkedList<>();
        for ( final List<SequenceForKmers> oneSampleWorth : pending.values() ) result.addAll(oneSampleWorth);
        return result;
    }

    /**
     * Get the collection of non-unique kmers from sequence for kmer size kmerSize
     * @param seqForKmers a sequence to get kmers from
     * @param kmerSize the size of the kmers
     * @return a non-null collection of non-unique kmers in sequence
     */
    protected static Collection<Kmer> determineNonUniqueKmers(final SequenceForKmers seqForKmers, final int kmerSize) {
        // count up occurrences of kmers within each read
        final KMerCounter counter = new KMerCounter(kmerSize);
        final int stopPosition = seqForKmers.stop - kmerSize;
        for ( int i = 0; i <= stopPosition; i++ ) {
            final Kmer kmer = new Kmer(seqForKmers.sequence, i, kmerSize);
            counter.addKmer(kmer, 1);
        }

        return counter.getKmersWithCountsAtLeast(2);
    }

    @Override
    public SeqGraph convertToSequenceGraph() {
        buildGraphIfNecessary();
        return super.convertToSequenceGraph();
    }

    private void increaseCountsInMatchedKmers(final SequenceForKmers seqForKmers,
                                              final MultiDeBruijnVertex vertex,
                                              final byte[] originalKmer,
                                              final int offset) {
        if ( offset == -1 ) return;

        for ( final MultiSampleEdge edge : incomingEdgesOf(vertex) ) {
            final MultiDeBruijnVertex prev = getEdgeSource(edge);
            final byte suffix = prev.getSuffix();
            final byte seqBase = originalKmer[offset];
//            logger.warn(String.format("Increasing counts for %s -> %s via %s at %d with suffix %s vs. %s",
//                    prev, vertex, edge, offset, (char)suffix, (char)seqBase));
            if ( suffix == seqBase && (increaseCountsThroughBranches || inDegreeOf(vertex) == 1) ) {
                edge.incMultiplicity(seqForKmers.count);
                increaseCountsInMatchedKmers(seqForKmers, prev, originalKmer, offset-1);
            }
        }
    }

    /**
     * Get the vertex for the kmer in sequence starting at start
     * @param sequence the sequence
     * @param start the position of the kmer start
     * @return a non-null vertex
     */
    private MultiDeBruijnVertex getOrCreateKmerVertex(final byte[] sequence, final int start) {
        final Kmer kmer = new Kmer(sequence, start, kmerSize);
        final MultiDeBruijnVertex vertex = getUniqueKmerVertex(kmer, true);
        return ( vertex != null ) ? vertex : createVertex(kmer);
    }

    /**
     * Get the unique vertex for kmer, or null if not possible.
     *
     * @param allowRefSource if true, we will allow kmer to match the reference source vertex
     * @return a vertex for kmer, or null if it's not unique
     */
    private MultiDeBruijnVertex getUniqueKmerVertex(final Kmer kmer, final boolean allowRefSource) {
        if ( ! allowRefSource && kmer.equals(refSource) ) return null;

        return uniqueKmers.get(kmer);
    }


    /**
     * Create a new vertex for kmer.  Add it to the uniqueKmers map if appropriate.
     *
     * kmer must not have a entry in unique kmers, or an error will be thrown
     *
     * @param kmer the kmer we want to create a vertex for
     * @return the non-null created vertex
     */
    private MultiDeBruijnVertex createVertex(final Kmer kmer) {
        final MultiDeBruijnVertex newVertex = new MultiDeBruijnVertex(kmer.bases());
        final int prevSize = vertexSet().size();
        addVertex(newVertex);

        // make sure we aren't adding duplicates (would be a bug)
        if ( vertexSet().size() != prevSize + 1) throw new IllegalStateException("Adding vertex " + newVertex + " to graph didn't increase the graph size");

        // add the vertex to the unique kmer map, if it is in fact unique
        if ( ! nonUniqueKmers.contains(kmer) && ! uniqueKmers.containsKey(kmer) ) // TODO -- not sure this last test is necessary
            uniqueKmers.put(kmer, newVertex);

        return newVertex;
    }

    /**
     * Workhorse routine of the assembler.  Given a sequence whose last vertex is anchored in the graph, extend
     * the graph one bp according to the bases in sequence.
     *
     * @param prevVertex a non-null vertex where sequence was last anchored in the graph
     * @param sequence the sequence we're threading through the graph
     * @param kmerStart the start of the current kmer in graph we'd like to add
     * @param count the number of observations of this kmer in graph (can be > 1 for GGA)
     * @param isRef is this the reference sequence?
     * @return a non-null vertex connecting prevVertex to in the graph based on sequence
     */
    private MultiDeBruijnVertex extendChainByOne(final MultiDeBruijnVertex prevVertex, final byte[] sequence, final int kmerStart, final int count, final boolean isRef) {
        final Set<MultiSampleEdge> outgoingEdges = outgoingEdgesOf(prevVertex);

        final int nextPos = kmerStart + kmerSize - 1;
        for ( final MultiSampleEdge outgoingEdge : outgoingEdges ) {
            final MultiDeBruijnVertex target = getEdgeTarget(outgoingEdge);
            if ( target.getSuffix() == sequence[nextPos] ) {
                // we've got a match in the chain, so simply increase the count of the edge by 1 and continue
                outgoingEdge.incMultiplicity(count);
                return target;
            }
        }

        // none of our outgoing edges had our unique suffix base, so we check for an opportunity to merge back in
        final Kmer kmer = new Kmer(sequence, kmerStart, kmerSize);
        final MultiDeBruijnVertex uniqueMergeVertex = getUniqueKmerVertex(kmer, false);

        if ( isRef && uniqueMergeVertex != null )
            throw new IllegalStateException("Found a unique vertex to merge into the reference graph " + prevVertex + " -> " + uniqueMergeVertex);

        // either use our unique merge vertex, or create a new one in the chain
        final MultiDeBruijnVertex nextVertex = uniqueMergeVertex == null ? createVertex(kmer) : uniqueMergeVertex;
        addEdge(prevVertex, nextVertex, ((MyEdgeFactory)getEdgeFactory()).createEdge(isRef, count));
        return nextVertex;
    }

    /**
     * Add the given read to the sequence graph.  Ultimately the read will get sent through addSequence(), but first
     * this method ensures we only use high quality bases and accounts for reduced reads, etc.
     *
     * @param read a non-null read
     */
    protected void addRead(final SAMRecord read) {
        final byte[] sequence = read.getReadBases();
        final byte[] qualities = read.getBaseQualities();

        int lastGood = -1; // the index of the last good base we've seen
        for( int end = 0; end <= sequence.length; end++ ) {
            if ( end == sequence.length || ! baseIsUsableForAssembly(sequence[end], qualities[end]) ) {
                // the first good base is at lastGood, can be -1 if last base was bad
                final int start = lastGood;
                // the stop base is end - 1 (if we're not at the end of the sequence)
                final int len = end - start;

                if ( start != -1 && len >= kmerSize ) {
                    // if the sequence is long enough to get some value out of, add it to the graph
                    final String name = read.getReadName() + "_" + start + "_" + end;
                    addSequence(name, read.getReadGroup().getSample(), read.getReadBases(), start, end, 1, false);
                }

                lastGood = -1; // reset the last good base
            } else if ( lastGood == -1 ) {
                lastGood = end; // we're at a good base, the last good one is us
            }
        }
    }

    /**
     * Determines whether a base can safely be used for assembly.
     * Currently disallows Ns and/or those with low quality
     *
     * @param base  the base under consideration
     * @param qual  the quality of that base
     * @return true if the base can be used for assembly, false otherwise
     */
    protected boolean baseIsUsableForAssembly(final byte base, final byte qual) {
        return base != BaseUtils.Base.N.base && qual >= minBaseQualityToUseInAssembly;
    }

    /**
     * Get the set of non-unique kmers in this graph.  For debugging purposes
     * @return a non-null set of kmers
     */
    protected Set<Kmer> getNonUniqueKmers() {
        return nonUniqueKmers;
    }

    @Override
    public String toString() {
        return "ReadThreadingAssembler{" +
                "kmerSize=" + kmerSize +
                '}';
    }


    @Override
    public MultiDeBruijnVertex findKmer(final Kmer k) {
        return uniqueKmers.get(k);
    }

    /*************************************************************
     * Simple string representation support for testing purposes *
     *************************************************************/

    private static final Pattern PROPERTIES_PATTERN = Pattern.compile("^\\s*\\[[^\\]]*\\]");
    private static final Pattern PATH_PATTERN = Pattern.compile("\\{((\\S+):)?([^\\}]*)\\}");
    private static final Pattern KMERSIZE_EXTRACTOR_PATTERN = Pattern.compile("^\\s*\\[[^\\]]*(ks|kmerSize)\\s*=\\s*(\\d+)\\s*[,\\]]");


    /**
     * Constructs a read-threading-graph for a string representation.
     *
     * <p>
     *     Note: only used for testing.
     *     Checkout {@see HaplotypeGraphUnitTest} for examples.
     * </p>
     * @param s the string representation of the graph {@code null}.
     */
    public ReadThreadingGraph(final String s) {
        super(kmerSizeFromString(s),new MyEdgeFactory(1));
        debugGraphTransformations = false;
        minBaseQualityToUseInAssembly = 0;
        applyString(s);
        alreadyBuilt = true;
    }

    /**
     * Obtain the kmer size for the string representation.
     * @param str the source string representation.
     * @return 1 or greater.
     * @throws IllegalArgumentException if {@code} str does not contain a valid representation.
     */
    private static int kmerSizeFromString(final String str) {
        final Matcher matcher = KMERSIZE_EXTRACTOR_PATTERN.matcher(str);
        if (matcher.find()) {
            return Integer.parseInt(matcher.group(2));
        } else
            throw new IllegalArgumentException("the input graph spec does not indicate the kmerSize");
    }

    /**
     * Apply description string into the graph.
     *
     * <p>
     *     Note: this is done just for testing purposes.
     *     Checkout {@see HaplotypeGraphUnitTest} for examples.
     * </p>
     * @param str the string representation.
     */
    private void applyString(final String str) {
        final Matcher propertiesSectionMatcher = PROPERTIES_PATTERN.matcher(str);
        final int pathStart = propertiesSectionMatcher.find() ? propertiesSectionMatcher.end() : 0;

        final String pathString = str.substring(pathStart);
        final Matcher pathMatcher = PATH_PATTERN.matcher(pathString);

        boolean referenceFound = false;
        final Map<String,MultiDeBruijnVertex> vertexById = new HashMap<>();

        // Loop between path strings and add them one by one.
        while (pathMatcher.find()) {
            final String label = pathMatcher.group(2);
            final boolean isReference = (label != null && label.equals("REF"));
            if (referenceFound) {
                if (isReference) {
                    throw new IllegalArgumentException("there are two reference paths");
                }
            } else if ( isReference ) {
                referenceFound = true;
            }

            // Divide each path into its elements getting a list of sequences and labels if applies:
            final String elementsString = pathMatcher.group(3);
            final String[] elements = elementsString.split("\\s*->\\s*");
            if (elements.length == 0) {
                throw new IllegalArgumentException("empty path not allowed");
            }
            final String[] seqs = new String[elements.length];
            final String[] ids = new String[elements.length];
            for (int i = 0; i < elements.length; i++) {
                ids[i] = pathElementId(elements[i]);
                seqs[i] = pathElementSeq(elements[i]);
                if (seqs[i].isEmpty() && ids[i] == null) {
                    throw new IllegalArgumentException("path with empty element without an id");
                }
            }
            final boolean isSource =  ids[0] == null || !vertexById.containsKey(ids[0]);
            if (isSource && seqs[0].length() != kmerSize) {
                throw new IllegalArgumentException("source sequence length must be the same as the kmerSize "
                        + ids[0] + " " + seqs[0] + " " + pathMatcher.group());
            }
            final MultiDeBruijnVertex firstVertex;
            if (ids[0] != null && vertexById.containsKey(ids[0])) {
                firstVertex = vertexById.get(ids[0]);
            } else {
                firstVertex = new MultiDeBruijnVertex(seqs[0].getBytes());
                addVertex(firstVertex);
                if (ids[0] != null) {
                    vertexById.put(ids[0], firstVertex);
                }
            }
            if (!seqs[0].isEmpty() &&
                    ((isSource && !firstVertex.getSequenceString().equals(seqs[0]))
                            || (!isSource && firstVertex.getSuffix() != seqs[0].getBytes()[0]))) {
                throw new IllegalArgumentException("mismatched first element sequence");
            }

            MultiDeBruijnVertex lastVertex = firstVertex;
            for (int i = 1; i < elements.length; i++) {
                if (seqs[i].length() > 1) {
                    throw new IllegalArgumentException("non-source vertex sequence must have length 1");
                }
                final MultiDeBruijnVertex nextVertex;
                if (ids[i] == null || !vertexById.containsKey(ids[i])) {
                    final Set<MultiDeBruijnVertex> nextVertices = getNextVertices(lastVertex,seqs[i].getBytes()[0]);
                    if (nextVertices.size() == 0) {
                        nextVertex = new MultiDeBruijnVertex(extendSequence(lastVertex.getSequence(),seqs[i].getBytes()[0]));
                        addVertex(nextVertex);
                    } else {
                        nextVertex = nextVertices.iterator().next();
                    }
                    if (ids[i] != null)
                        vertexById.put(ids[i],nextVertex);
                } else {
                    nextVertex = vertexById.get(ids[i]);
                }
                final MultiSampleEdge edge = addEdge(lastVertex,nextVertex);
                if (isReference) {
                    edge.setIsRef(true);
                }
                lastVertex = nextVertex;
            }
        }
    }

    private static String pathElementId(final String element) {
        final int openBracketPosition = element.indexOf('(');

        if (openBracketPosition == -1) {
            return null;
        }

        final int closeBracketPosition = element.lastIndexOf(')');
        if (closeBracketPosition == -1) {
            throw new IllegalArgumentException("non-closed id parantesys found in element: " + element);
        }
        final String result = element.substring(openBracketPosition + 1,closeBracketPosition).trim();
        if (result.isEmpty()) {
            throw new IllegalArgumentException("empty id found in element: " + element);
        }
        return result;
    }

    /**
     * Returns the lenght of a path element in the string representation.
     * @param element the query element.
     * @return 0 or greater.
     */
    private static String pathElementSeq(final String element) {
        final int parentesysPos = element.indexOf('(');

        if (parentesysPos == -1)
            return element.trim();

        return element.substring(0,parentesysPos).trim();
    }

    /**
     * Add a base to the end of a byte sequence.
     * @param sequence sequence where to add the base to.
     * @param b base to add.
     * @return never {@code null}, a new array each time.
     */
    private static byte[] extendSequence(final byte[] sequence, final byte b) {
        final byte[] result = new byte[sequence.length];
        System.arraycopy(sequence, 1, result, 0, sequence.length - 1);
        result[result.length - 1] = b;
        return result;
    }
}