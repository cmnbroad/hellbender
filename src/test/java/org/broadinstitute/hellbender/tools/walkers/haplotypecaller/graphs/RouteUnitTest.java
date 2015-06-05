package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.jgrapht.EdgeFactory;
import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class RouteUnitTest extends BaseTest {

       @Test(dataProvider="slicePrefixTestData")
       public void testSplicePrefix(final Route<BaseVertex,BaseEdge> route) {
           final int routeLength = route.length();
           for (int i = 0; i < routeLength; i++) {
               final Route<BaseVertex,BaseEdge> spliced = route.splicePrefix(i);
               Assert.assertEquals(spliced.length(), route.length() - i);
               final List<BaseEdge> routeEdges = route.getEdges();
               final List<BaseEdge> expectedSlicedEdges = routeEdges.subList(i,routeLength);
               Assert.assertEquals(spliced.getEdges(), expectedSlicedEdges);
           }
       }

       @Test(dataProvider="isSuffixTestData")
       public void testIsSuffix(final Route<BaseVertex,BaseEdge> route, final Path<BaseVertex,BaseEdge> path, final boolean expectedResult) {
          Assert.assertEquals(route.isSuffix(path), expectedResult);
       }

       @DataProvider(name="isSuffixTestData")
       public Iterator<Object[]> isSuffixTestData() {
           return IS_SUFFIX_TEST_DATA.iterator();
       }

       @DataProvider(name="slicePrefixTestData")
       public Iterator<Object[]> slicePrefixTestData() {
           return Arrays.asList(SLICE_PREFIX_TEST_DATA).iterator();
       }

       private static final int[] TEST_EDGE_PAIRS1 = new int[] {
                      3 , 4,
                          4 , 5,
                              5, 7,
                                 7, 8,
                                    8, 9,
                          4 , 6,
                              6,       9,
                                       9,     11,
                                              11, 12,
       };

       private static final int[] TEST_EDGE_PAIRS = new int[] {
              1 , 2,
                  2 , 3,
                      3 , 4,
                          4 , 5,
                              5, 7,
                                 7, 8,
                                    8, 9,
                          4 , 6,
                              6,       9,
                                       9, 10,
                                          10, 11,
                                              11, 12,
                  2,          5,
                              5,                  12,

                      3,                             13,
                                                     13, 14,
                                                         14, 15
       };

    public static final EdgeFactory<BaseVertex, BaseEdge> TEST_GRAPH_EDGE_FACTORY = new EdgeFactory<BaseVertex, BaseEdge>() {
        @Override
        public BaseEdge createEdge(final BaseVertex baseVertex, final BaseVertex baseVertex2) {
            return new BaseEdge(false, 0);
        }
    };


    private static Map<Integer, BaseVertex> vertexByInteger = new HashMap<>();
    private static final BaseGraph<BaseVertex, BaseEdge> TEST_GRAPH = new BaseGraph<>(1, TEST_GRAPH_EDGE_FACTORY);
    private static final List<Object[]> IS_SUFFIX_TEST_DATA;

    private static final Object[][] SLICE_PREFIX_TEST_DATA;

    static {
        for (int i = 0; i < TEST_EDGE_PAIRS.length; i += 2) {
            final int sourceInteger = TEST_EDGE_PAIRS[i];
            final int targetInteger = TEST_EDGE_PAIRS[i + 1];
            final BaseVertex sourceVertex = resolveVertexByInteger(sourceInteger);
            final BaseVertex targetVertex = resolveVertexByInteger(targetInteger);
            TEST_GRAPH.addEdge(sourceVertex, targetVertex);
        }
        Assert.assertEquals(1, TEST_GRAPH.getSources().size());
        final Deque<Path<BaseVertex,BaseEdge>> pendingPaths = new LinkedList<>();
        final Deque<Route<BaseVertex,BaseEdge>> pendingRoutes = new LinkedList<>();
        final List<Path<BaseVertex,BaseEdge>> allPossiblePaths = new LinkedList<>();
        final List<Route<BaseVertex,BaseEdge>> allPossibleRoutes = new LinkedList<>();
        for (final BaseVertex vertex : TEST_GRAPH.vertexSet()) {
            pendingPaths.add(new Path(vertex, TEST_GRAPH));
            pendingRoutes.add(new Route(vertex,TEST_GRAPH));
        }
        while (!pendingPaths.isEmpty()) { // !pendingRoutes.isEmpty();
            final Path<BaseVertex,BaseEdge> path = pendingPaths.remove();
            final Route<BaseVertex,BaseEdge> route = pendingRoutes.remove();
            final BaseVertex lastVertex = path.getLastVertex();
            allPossiblePaths.add(path);
            allPossibleRoutes.add(route);

            if (allPossiblePaths.size() % 100 == 0)
                Reporter.log("" + allPossiblePaths.size(), true);
            for (final BaseEdge edge : TEST_GRAPH.outgoingEdgesOf(lastVertex))
                pendingPaths.add(new Path<>(path,edge));
            for (final BaseEdge edge : TEST_GRAPH.outgoingEdgesOf(lastVertex))
                pendingRoutes.add(new Route<>(route,edge));
        }

        final int numberOfPaths = allPossiblePaths.size();
        final boolean[][] isSuffix = buildIsSuffixMatrix(allPossiblePaths, numberOfPaths);
        IS_SUFFIX_TEST_DATA = createTestData(allPossiblePaths,allPossibleRoutes,isSuffix);
        SLICE_PREFIX_TEST_DATA = createSlicePrefixTestData(allPossibleRoutes);
    }

    private static Object[][] createSlicePrefixTestData(List<Route<BaseVertex, BaseEdge>> allPossibleRoutes) {
        final Object[][] result = new Object[allPossibleRoutes.size()][1];
        final Object[] routes = allPossibleRoutes.toArray();
        for (int i = 0; i < result.length; i++)
            result[i][0] = routes[i];
        return result;
    }

    private static boolean[][] buildIsSuffixMatrix(final List<Path<BaseVertex, BaseEdge>> allPossiblePaths, final int numberOfPaths) {
        final boolean[][] isSuffix = new boolean[numberOfPaths][numberOfPaths];
        final ListIterator<Path<BaseVertex,BaseEdge>> iIterator = allPossiblePaths.listIterator();
        for (int i = 0; i < numberOfPaths; i++) {
            isSuffix[i][i] = true;
            final ListIterator<Path<BaseVertex,BaseEdge>> jIterator = allPossiblePaths.listIterator(i + 1);
            final Path<BaseVertex,BaseEdge> iPath = iIterator.next();
            for (int j = i + 1; j < numberOfPaths; j++) {
                final Path<BaseVertex,BaseEdge> jPath = jIterator.next();
                if (iPath.getLastVertex() != jPath.getLastVertex()) {
                    isSuffix[i][j] = isSuffix[j][i] = false;
                } else {
                    isSuffix[i][j] = isSuffix[j][i] = true; // let assume they are suffix of each other by default.
                    final Path<BaseVertex,BaseEdge> shortPath;
                    final Path<BaseVertex,BaseEdge> longPath;
                    if (iPath.getEdges().size() <= jPath.getEdges().size()) {
                        shortPath = iPath;
                        longPath = jPath;
                    } else {
                        longPath = iPath;
                        shortPath = jPath;
                    }
                    final ListIterator<BaseEdge> longPathEdgesIterator = longPath.getEdges().listIterator(longPath.getEdges().size());
                    final ListIterator<BaseEdge> shortPathEdgesIterator = shortPath.getEdges().listIterator(shortPath.getEdges().size());

                    while (shortPathEdgesIterator.hasPrevious()) {
                        final BaseEdge shortEdge = shortPathEdgesIterator.previous();
                        final BaseEdge longEdge = longPathEdgesIterator.previous();
                        if (shortEdge != longEdge) {
                           isSuffix[i][j] = isSuffix[j][i] = false;
                            break;
                        }
                    }
                    if (isSuffix[i][j]) {
                        if (longPathEdgesIterator.hasPrevious()) {
                            if (longPath == iPath)
                                isSuffix[j][i] = false;
                            else
                                isSuffix[i][j] = false;
                        }
                    }
                }

            }
        }
        return isSuffix;
    }

    @SuppressWarnings({"rawtypes","unchecked"})
    private static List<Object[]> createTestData(final List<Path<BaseVertex, BaseEdge>> allPossiblePaths, final List<Route<BaseVertex, BaseEdge>> allPossibleRoutes, final boolean[][] isSuffix) {
        final List<Object[]> result = new ArrayList<>(allPossiblePaths.size() * allPossiblePaths.size() * 2 );

        final Path<BaseVertex,BaseEdge>[] allPaths = allPossiblePaths.toArray(new Path[allPossiblePaths.size()]);
        final Route<BaseVertex,BaseEdge>[] allRoutes = allPossibleRoutes.toArray(new Route[allPossibleRoutes.size()]);
        final int numberOfPaths = allPaths.length;
        for (int i = 0; i < numberOfPaths; i++)
            for (int j = 0; j < numberOfPaths; j++) {
                result.add(new Object[] { allRoutes[i], allPaths[j], isSuffix[i][j] });
                result.add(new Object[] { allRoutes[i], allRoutes[j], isSuffix[i][j] });
                result.add(new Object[] { allRoutes[i], inverseRebuild(allRoutes[j]), isSuffix[i][j]});
            }

        return result;
    }

    private static Route<BaseVertex,BaseEdge> inverseRebuild(final Route<BaseVertex,BaseEdge> original) {
        final ListIterator<BaseEdge> it = original.getEdges().listIterator(original.length());
        Route<BaseVertex,BaseEdge> result = new Route<>(original.getLastVertex(),original.getGraph());
        while (it.hasPrevious()) {
            result = new Route<>(it.previous(),result);
        }
        return result;
    }

    private static BaseVertex resolveVertexByInteger(final int targetInteger) {
        if (vertexByInteger.containsKey(targetInteger))
            return vertexByInteger.get(targetInteger);
        else {
            int value = targetInteger;
            final StringBuffer stringBuffer = new StringBuffer();
            while (value > 0) {
               int c = value % 4;
               switch (c) {
                   case 0: stringBuffer.append('A'); break;
                   case 1: stringBuffer.append('C'); break;
                   case 2: stringBuffer.append('G'); break;
                   case 3: stringBuffer.append('T'); break;
               }
               value = value / 4;
            }
            if (stringBuffer.length() == 0) stringBuffer.append('A');
            final byte[] sequence = stringBuffer.reverse().toString().getBytes();
            final BaseVertex result = new BaseVertex(sequence);
            vertexByInteger.put(targetInteger, result);
            TEST_GRAPH.addVertex(result);
            return result;
        }

    }


}
