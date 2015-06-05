package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class BaseEdgeUnitTest extends BaseTest {
    @DataProvider(name = "EdgeCreationData")
    public Object[][] makeMyDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        for ( final int multiplicity : Arrays.asList(1, 2, 3) ) {
            for ( final boolean isRef : Arrays.asList(true, false) ) {
                tests.add(new Object[]{isRef, multiplicity});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "EdgeCreationData")
    public void testBasic(final boolean isRef, final int mult) {
        final BaseEdge e = new BaseEdge(isRef, mult);
        Assert.assertEquals(e.isRef(), isRef);
        Assert.assertEquals(e.getMultiplicity(), mult);

        e.setIsRef(!isRef);
        Assert.assertEquals(e.isRef(), !isRef);

        e.setMultiplicity(mult + 1);
        Assert.assertEquals(e.getMultiplicity(), mult + 1);

        e.incMultiplicity(2);
        Assert.assertEquals(e.getMultiplicity(), mult + 3);

        final BaseEdge copy = e.copy();
        Assert.assertEquals(copy.isRef(), e.isRef());
        Assert.assertEquals(copy.getMultiplicity(), e.getMultiplicity());
    }

    @Test
    public void testEdgeWeightComparator() {
        final BaseEdge e10 = new BaseEdge(false, 10);
        final BaseEdge e5 = new BaseEdge(true, 5);
        final BaseEdge e2 = new BaseEdge(false, 2);
        final BaseEdge e1 = new BaseEdge(false, 1);

        final List<BaseEdge> edges = new ArrayList<>(Arrays.asList(e1, e2, e5, e10));
        Collections.sort(edges, new BaseEdge.EdgeWeightComparator());
        Assert.assertEquals(edges.get(0), e10);
        Assert.assertEquals(edges.get(1), e5);
        Assert.assertEquals(edges.get(2), e2);
        Assert.assertEquals(edges.get(3), e1);
    }

    @Test
    public void testMax() {
        for ( final boolean firstIsRef : Arrays.asList(true, false) ) {
            for ( final boolean secondIsRef : Arrays.asList(true, false) ) {
                for ( final int firstMulti : Arrays.asList(1, 4) ) {
                    for ( final int secondMulti : Arrays.asList(2, 3) ) {
                        final BaseEdge expected = new BaseEdge(firstIsRef || secondIsRef, Math.max(firstMulti, secondMulti));
                        final BaseEdge actual = new BaseEdge(firstIsRef, firstMulti).max(new BaseEdge(secondIsRef, secondMulti));
                        Assert.assertEquals(actual.getMultiplicity(), expected.getMultiplicity());
                        Assert.assertEquals(actual.isRef(), expected.isRef());
                    }
                }
            }
        }
    }
}
