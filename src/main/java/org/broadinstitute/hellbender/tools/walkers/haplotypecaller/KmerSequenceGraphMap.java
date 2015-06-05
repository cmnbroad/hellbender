package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseVertex;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.KmerSearchableGraph;

import java.util.*;

/**
 * Contains information as to how a kmer sequence maps to an (assembly) graph.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.com&gt;
 */
public class KmerSequenceGraphMap<V extends BaseVertex,E extends BaseEdge> {

    protected final KmerSequence sequence;
    protected final KmerSearchableGraph<V,E> graph;
    protected final int kmerSize;

    private List<V> vertexList;

    private List<V> vertexMatchOnlyList;

    private Set<V> vertexSet;

    private Map<V,Integer> vertexOffset;

    //private List<VertexSegment<V>> vertexSegmentList;

    /**
     * Constructs a new Kmer sequence graph map give the graph and sequence.
     *
     * @param g the graph to map to.
     * @param s the sequence to map.
     * @throws NullPointerException if either the graph or the sequence is <code>null</code>.
     * @throws IllegalArgumentException if the kmer sizes of the input graph and sequence are not the same.
     */
    public KmerSequenceGraphMap(final KmerSearchableGraph<V, E> g, final KmerSequence s) {
        if (s.kmerSize() != g.getKmerSize()) {
            throw new IllegalArgumentException("kmer size for the graph (" + g.getKmerSize() + ") and the sequence (" + s.kmerSize() + ") are different");
        }
        sequence = s;
        graph = g;
        kmerSize = s.kmerSize();
    }

    /**
     * Vertices that form part of the kmer sequence path along the graph.
     *
     * <p>The ith position in the resulting list corresponds to the ith kmer in the sequence</p>.
     *
     * <p>
     *     The resulting list will contain null values for those kmers where there is no unique kmer match in the
     *     graph.
     * </p>
     *
     * @return never {@code null}
     */
    public List<V> vertexList() {
        if (vertexList == null)
            buildVertexCollections();
        return vertexList;
    }

    /**
     * Vertices that form part of the kmer sequence path along the graph.
     *
     * <p> Only contains unique kmer vertices where the non-unique ones have been sliced out from the list</p>
     *
     * @return never {@code null}
     */
    public List<V> vertexMatchOnlyList() {
        if (vertexMatchOnlyList == null) {
            buildVertexCollections();
        }
        return vertexMatchOnlyList;
    }


    /**
     * Return a map from vertices to their kmer offset in the kmer sequence.
     * @return never {@code null}
     */
    public Map<V,Integer> vertexOffset() {
        if (vertexOffset == null) {
            buildVertexCollections();
        }
        return vertexOffset;
    }

    /**
     * Set of all vertices with unique kmers in the kmer sequence.
     * <p>
     *     This structure is more appropriate to query whether a vertex belong or not to such a set.
     * </p>
     * @return never {@code null}.
     */
    public Set<V> vertexSet() {
        if (vertexSet == null) {
            buildVertexCollections();
        }
        return vertexSet;
    }

    /**
     * Updates vertex structures.
     */
    protected void buildVertexCollections() {
        @SuppressWarnings("unchecked")
        final V[] result = (V[]) new BaseVertex[sequence.size()];
        final Set<V> set = new HashSet<>(sequence.size());
        final Map<V,Integer> posMap = new HashMap<>(sequence.size());
        @SuppressWarnings("unchecked")
        final V[] matchOnly = (V[]) new BaseVertex[sequence.size()];
        int next = 0;
        int matchOnlyNext = 0;
        for (int i = 0; i < sequence.size(); i++) {
            final Kmer k = sequence.get(i);
            final V v = graph.findKmer(k);
            if (v != null) {
                set.add(v);
                posMap.put(v,i);
                matchOnly[matchOnlyNext++] = v;
            }
            result[next++] = v;
        }
        vertexList = Arrays.asList(result);
        vertexMatchOnlyList = Arrays.asList(Arrays.copyOf(matchOnly, matchOnlyNext));
        vertexSet = Collections.unmodifiableSet(set);
        vertexOffset = Collections.unmodifiableMap(posMap);
    }

    /**
     * Returns the list of kmers in the sequence that do not have a unique mapping on the graph.
     * @return never {@code null}
     */
    @SuppressWarnings("unused")
    public List<Kmer> missingKmers() {
        if (vertexList == null) {
            buildVertexCollections();
        }
        if (vertexList.size() == vertexMatchOnlyList.size()) {
            return Collections.emptyList();
        } else {
            final List<Kmer> result = new ArrayList<>(vertexList.size() - vertexMatchOnlyList.size());
            final int size = sequence.size();
            for (int i = 0; i < vertexList.size(); i++) {
                if (vertexList.get(i) == null) {
                    result.add(sequence.get(i));
                }
            }
            return result;
        }
    }

}
