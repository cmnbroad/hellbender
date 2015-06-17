package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.KBestHaplotype;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.KBestHaplotypeFinder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public final class DanglingChainMergingGraphUnitTest extends BaseTest {

    public static byte[] getBytes(final String alignment) {
        return alignment.replace("-","").getBytes();
    }

    @DataProvider(name = "DanglingTails")
    public Object[][] makeDanglingTailsData() {
        List<Object[]> tests = new ArrayList<>();

        // add 1M to the expected CIGAR because it includes the previous (common) base too
        tests.add(new Object[]{"AAAAAAAAAA", "CAAA", "5M", true, 3});                  // incomplete haplotype
        tests.add(new Object[]{"AAAAAAAAAA", "CAAAAAAAAAA", "1M1I10M", true, 10});     // insertion
        tests.add(new Object[]{"CCAAAAAAAAAA", "AAAAAAAAAA", "1M2D10M", true, 10});    // deletion
        tests.add(new Object[]{"AAAAAAAA", "CAAAAAAA", "9M", true, 7});                // 1 snp
        tests.add(new Object[]{"AAAAAAAA", "CAAGATAA", "9M", true, 2});                // several snps
        tests.add(new Object[]{"AAAAA", "C", "1M4D1M", false, -1});                    // funky SW alignment
        tests.add(new Object[]{"AAAAA", "CA", "1M3D2M", false, 1});                    // very little data
        tests.add(new Object[]{"AAAAAAA", "CAAAAAC", "8M", true, -1});                 // ends in mismatch
        tests.add(new Object[]{"AAAAAA", "CGAAAACGAA", "1M2I4M2I2M", false, 0});       // alignment is too complex
        tests.add(new Object[]{"AAAAA", "XXXXX", "1M5I", false, -1});                  // insertion

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "DanglingTails")
    public void testDanglingTails(final String refEnd,
                                  final String altEnd,
                                  final String cigar,
                                  final boolean cigarIsGood,
                                  final int mergePointDistanceFromSink) {

        final int kmerSize = 15;

        // construct the haplotypes
        final String commonPrefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
        final String ref = commonPrefix + refEnd;
        final String alt = commonPrefix + altEnd;

        // create the graph and populate it
        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(kmerSize);
        rtgraph.addSequence("ref", ref.getBytes(), true);
        final SAMRecord read = ArtificialSAMUtils.createArtificialRead(alt.getBytes(), Utils.dupBytes((byte) 30, alt.length()), alt.length() + "M");
        rtgraph.addRead(read);
        rtgraph.buildGraphIfNecessary();

        // confirm that we have just a single dangling tail
        MultiDeBruijnVertex altSink = null;
        for ( final MultiDeBruijnVertex v : rtgraph.vertexSet() ) {
            if ( rtgraph.isSink(v) && !rtgraph.isReferenceNode(v) ) {
                Assert.assertTrue(altSink == null, "We found more than one non-reference sink");
                altSink = v;
            }
        }

        Assert.assertTrue(altSink != null, "We did not find a non-reference sink");

        // confirm that the SW alignment agrees with our expectations
        final ReadThreadingGraph.DanglingChainMergeHelper result = rtgraph.generateCigarAgainstDownwardsReferencePath(altSink, 0, 4);

        if ( result == null ) {
            Assert.assertFalse(cigarIsGood);
            return;
        }

        Assert.assertTrue(cigar.equals(result.cigar.toString()), "SW generated cigar = " + result.cigar.toString());

        // confirm that the goodness of the cigar agrees with our expectations
        Assert.assertEquals(rtgraph.cigarIsOkayToMerge(result.cigar, false, true), cigarIsGood);

        // confirm that the tail merging works as expected
        if ( cigarIsGood ) {
            final int mergeResult = rtgraph.mergeDanglingTail(result);
            Assert.assertTrue(mergeResult == 1 || mergePointDistanceFromSink == -1);

            // confirm that we created the appropriate edge
            if ( mergePointDistanceFromSink >= 0 ) {
                MultiDeBruijnVertex v = altSink;
                for ( int i = 0; i < mergePointDistanceFromSink; i++ ) {
                    if ( rtgraph.inDegreeOf(v) != 1 )
                        Assert.fail("Encountered vertex with multiple edges");
                    v = rtgraph.getEdgeSource(rtgraph.incomingEdgeOf(v));
                }
                Assert.assertTrue(rtgraph.outDegreeOf(v) > 1);
            }
        }
    }

    @Test
    public void testWholeTailIsInsertion() {
        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(10);
        final ReadThreadingGraph.DanglingChainMergeHelper result = new ReadThreadingGraph.DanglingChainMergeHelper(null, null, "AXXXXX".getBytes(), "AAAAAA".getBytes(), new TextCigarCodec().decode("5I1M"));
        final int mergeResult = rtgraph.mergeDanglingTail(result);
        Assert.assertEquals(mergeResult, 0);
    }

    @Test
    public void testGetBasesForPath() {

        final int kmerSize = 4;
        final String testString = "AATGGGGCAATACTA";

        final ReadThreadingGraph graph = new ReadThreadingGraph(kmerSize);
        graph.addSequence(testString.getBytes(), true);
        graph.buildGraphIfNecessary();

        final List<MultiDeBruijnVertex> vertexes = new ArrayList<>();
        MultiDeBruijnVertex v = graph.getReferenceSourceVertex();
        while ( v != null ) {
            vertexes.add(v);
            v = graph.getNextReferenceVertex(v);
        }

        final String resultForTails = new String(graph.getBasesForPath(vertexes, false));
        Assert.assertEquals(resultForTails, testString.substring(kmerSize - 1));
        final String resultForHeads = new String(graph.getBasesForPath(vertexes, true));
        Assert.assertEquals(resultForHeads, "GTAAGGGCAATACTA");  // because the source node will be reversed
    }

    @DataProvider(name = "DanglingHeads")
    public Object[][] makeDanglingHeadsData() {
        List<Object[]> tests = new ArrayList<>();

        // add 1M to the expected CIGAR because it includes the last (common) base too
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", true});        // 1 snp
        tests.add(new Object[]{"XXXAACCGGTTACGT", "XAAACCGGTTACGT", "7M", false});         // 1 snp
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "XAACGGTTACGT", "4M1D4M", false});   // deletion
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AYYCGGTTACGT", "8M", true});        // 2 snps
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGTAA", "AYCYGGTTACGTAA", "9M", true});    // 2 snps
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AYCGGTTACGT", "7M", true});         // very little data
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "YCCGGTTACGT", "6M", true});         // begins in mismatch

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "DanglingHeads")
    public void testDanglingHeads(final String ref,
                                  final String alt,
                                  final String cigar,
                                  final boolean shouldBeMerged) {

        final int kmerSize = 5;

        // create the graph and populate it
        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(kmerSize);
        rtgraph.addSequence("ref", ref.getBytes(), true);
        final SAMRecord read = ArtificialSAMUtils.createArtificialRead(alt.getBytes(), Utils.dupBytes((byte) 30, alt.length()), alt.length() + "M");
        rtgraph.addRead(read);
        rtgraph.setMaxMismatchesInDanglingHead(10);
        rtgraph.buildGraphIfNecessary();

        // confirm that we have just a single dangling head
        MultiDeBruijnVertex altSource = null;
        for ( final MultiDeBruijnVertex v : rtgraph.vertexSet() ) {
            if ( rtgraph.isSource(v) && !rtgraph.isReferenceNode(v) ) {
                Assert.assertTrue(altSource == null, "We found more than one non-reference source");
                altSource = v;
            }
        }

        Assert.assertTrue(altSource != null, "We did not find a non-reference source");

        // confirm that the SW alignment agrees with our expectations
        final ReadThreadingGraph.DanglingChainMergeHelper result = rtgraph.generateCigarAgainstUpwardsReferencePath(altSource, 0, 1);

        if ( result == null ) {
            Assert.assertFalse(shouldBeMerged);
            return;
        }

        Assert.assertTrue(cigar.equals(result.cigar.toString()), "SW generated cigar = " + result.cigar.toString());

        // confirm that the tail merging works as expected
        final int mergeResult = rtgraph.mergeDanglingHead(result);
        Assert.assertTrue(mergeResult > 0 || !shouldBeMerged);

        // confirm that we created the appropriate bubble in the graph only if expected
        rtgraph.cleanNonRefPaths();
        final SeqGraph seqGraph = rtgraph.convertToSequenceGraph();
        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(seqGraph, seqGraph.getReferenceSourceVertex(), seqGraph.getReferenceSinkVertex());
        Assert.assertEquals(paths.size(), shouldBeMerged ? 2 : 1);
    }
}
