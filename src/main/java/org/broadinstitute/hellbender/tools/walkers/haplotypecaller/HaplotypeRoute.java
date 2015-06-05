package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.Route;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;

import java.util.*;

/**
 * Graph route that represent an haplotype on the haplotype assembly graph.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.com&gt;
 */
public class HaplotypeRoute extends Route<MultiDeBruijnVertex,MultiSampleEdge> {

    protected final Set<MultiDeBruijnVertex> vertexSet;

    protected final Map<MultiDeBruijnVertex,Integer> vertexOrder;

    protected final Set<MultiDeBruijnVertex> forkAndJoins;

    /**
     * Constructs a HaplotypeRoute given its route.
     *
     * @param route the haplotype route.
     */
    public HaplotypeRoute(final Route<MultiDeBruijnVertex, MultiSampleEdge> route) {
        super(route);
        vertexOrder = new LinkedHashMap<>(route.length() + 1);
        int nextOrder = 0;
        vertexOrder.put(getFirstVertex(),nextOrder++);
        for (final MultiSampleEdge edge : edgesInOrder)
            vertexOrder.put(graph.getEdgeTarget(edge), nextOrder++);
        Route<MultiDeBruijnVertex,MultiSampleEdge> currentRoute = this;
        forkAndJoins = new HashSet<>(route.length());
        while (currentRoute != null) {
          if (currentRoute.lastVertexIsForkOrJoin())
              forkAndJoins.add(currentRoute.getLastVertex());
          currentRoute = currentRoute.getPrefixRouteWithLastVertexThatIsForkOrJoin();
        }
        vertexSet = Collections.unmodifiableSet(new HashSet<>(vertexOrder.keySet()));
    }



    @SuppressWarnings("unused")
    public Route<MultiDeBruijnVertex,MultiSampleEdge> subRoute(final MultiDeBruijnVertex start, final MultiDeBruijnVertex end) {
       final Integer startOrder = vertexOrder.get(start);
       final Integer endOrder = vertexOrder.get(end);
       if (startOrder == null || endOrder == null)
           return null;
       else if (startOrder > endOrder)
           return null;
       else {
           Route<MultiDeBruijnVertex,MultiSampleEdge> result = new Route<>(start,graph);
           for (final MultiSampleEdge edge : edgesInOrder.subList(startOrder,endOrder))
                result = new Route(result,edge);
           return result;
       }
    }

    /**
     * Returns the set of vertex on the route.
     * @return read only, never {@code null} vertex set.
     */
    public Set<MultiDeBruijnVertex> vertexSet() {
        return vertexSet;
    }


    /**
     * Returns the position of the vertex in the route.
     *
     * @param vertex the query vertex.
     *
     * @throws NullPointerException if {@code vertex} is {@code null}.
     *
     * @return -1 if there is no such a vertex in the route, otherwise a number between 0 and {@link #length()} - 1.
     */
    public int getVertexPosition(final MultiDeBruijnVertex vertex) {
        final Integer result = vertexOrder.get(vertex);
        return result == null ? -1 : result;
    }
}
