package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.Route;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.HaplotypeGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;

import java.util.*;

/**
 * Represents an event block in the graph.
 *
 * <p>
 *     Event block is defined as the non-trivial section of the haplotype-graph between two vertices along the
 *     reference route, that has at least one alternative route between those two vertices.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class EventBlock {

    private final HaplotypeGraph graph;

    private final MultiDeBruijnVertex source;

    private final int sourcePosition;

    private final MultiDeBruijnVertex sink;

    private final int sinkPosition;

    private Set<Route<MultiDeBruijnVertex,MultiSampleEdge>> routesAcross;


    /**
     * Constructs a event block given the base haplotype graph and the source and sink vertice (both included in the block)
     * @param graph the base haplotype graph.
     * @param source the starting vertex.
     * @param sink the ending vertex.
     *
     * @throws NullPointerException if any of the input is {@code null}.
     * @throws IllegalArgumentException if {@code source} or {@code sink} are not part of the graphs reference route,
     *      such a route does not exists, any of the vertices is not part of such a route or they are out of order.
     */
    public EventBlock(final HaplotypeGraph graph, final MultiDeBruijnVertex source, final MultiDeBruijnVertex sink) {
        if (graph == null) throw new NullPointerException("the graph cannot be null");
        if (source == null) throw new NullPointerException("the source vertex is null");
        if (sink == null) throw new NullPointerException("the sink node is null");
        this.graph = graph;
        this.source = source;
        this.sink = sink;
        final HaplotypeRoute route = graph.getReferenceRoute();
        if (route == null)
            throw new IllegalArgumentException("there is reference route in the graph");
        this.sourcePosition = route.getVertexPosition(source);
        this.sinkPosition = route.getVertexPosition(sink);
        if (sourcePosition == -1)
            throw new IllegalArgumentException("the source vertex does not belong to the reference route");
        if (sinkPosition == -1)
            throw new IllegalArgumentException("the sink vertex does not belong to the reference route");
        if (sourcePosition > sinkPosition)
            throw new IllegalArgumentException("source and sink vertices are out of order in reference route");
    }

    /**
     * Returns a reference to the event block graph.
     *
     * @return never {@code null}.
     */
    public HaplotypeGraph getGraph() {
        return graph;
    }

    /**
     * Returns a reference to the block starting vertex.
     *
     * @return never {@code null}.
     */
    public MultiDeBruijnVertex getSource() {
        return source;
    }

    /**
     * Returns the reference ot the end block vertex.
     *
     * @return never {@code null}.
     */
    public MultiDeBruijnVertex getSink() {
        return sink;
    }

    /**
     * Returns all possible routes between the event block start and end vertices.
     * @return never {@code null}, and unmodifiable route set.
     */
    public Set<Route<MultiDeBruijnVertex,MultiSampleEdge>> getRoutesAcross() {
        // catching:
        if (routesAcross != null) return routesAcross;

        final Set<Route<MultiDeBruijnVertex, MultiSampleEdge>> result = new HashSet<>(10); // 10 is rather generous.

        // bread-first iterative search for all paths.
        final Queue<Route<MultiDeBruijnVertex, MultiSampleEdge>> queue = new LinkedList<>();

        queue.add(new Route<>(source, graph)); // the seed is the empty route at the start vertex.

        final HaplotypeRoute referenceRoute = graph.getReferenceRoute();

        while (!queue.isEmpty()) {
            final Route<MultiDeBruijnVertex, MultiSampleEdge> route = queue.remove();
            final MultiDeBruijnVertex routeEndVertex = route.getLastVertex();

            if (routeEndVertex == sink) // bingo!!!
                result.add(route);
            else { // only queue promising extension of this route.
                final int routeEndPosition = referenceRoute.getVertexPosition(routeEndVertex);
                if (routeEndPosition == -1 || (routeEndPosition >= sourcePosition && routeEndPosition < sinkPosition))
                    for (final MultiSampleEdge e : graph.outgoingEdgesOf(routeEndVertex))
                        queue.add(new Route<>(route, e));
            }
        }
        return routesAcross = Collections.unmodifiableSet(result);
    }
}
