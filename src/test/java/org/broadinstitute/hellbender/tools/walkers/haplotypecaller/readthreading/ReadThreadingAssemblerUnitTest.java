package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.KBestHaplotype;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.KBestHaplotypeFinder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class ReadThreadingAssemblerUnitTest extends BaseTest {

    private final static boolean DEBUG = false;

    private static class TestAssembler {
        final ReadThreadingAssembler assembler;

        Haplotype refHaplotype;
        final List<SAMRecord> reads = new LinkedList<>();

        private TestAssembler(final int kmerSize) {
            this.assembler = new ReadThreadingAssembler(100000, Arrays.asList(kmerSize));
            assembler.setJustReturnRawGraph(true);
            assembler.setPruneFactor(0);
        }

        public void addSequence(final byte[] bases, final boolean isRef) {
            if ( isRef ) {
                refHaplotype = new Haplotype(bases, true);
            } else {
                final SAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, Utils.dupBytes((byte) 30, bases.length), bases.length + "M");
                reads.add(read);
            }
        }

        public SeqGraph assemble() {
            assembler.removePathsNotConnectedToRef = false; // needed to pass some of the tests
            assembler.setRecoverDanglingBranches(false); // needed to pass some of the tests
            assembler.setDebugGraphTransformations(true);
            final SeqGraph graph = assembler.assemble(reads, refHaplotype, Collections.<Haplotype>emptyList()).get(0).getGraph();
            if ( DEBUG ) graph.printGraph(new File("test.dot"), 0);
            return graph;
        }
    }

    private void assertLinearGraph(final TestAssembler assembler, final String seq) {
        final SeqGraph graph = assembler.assemble();
        graph.simplifyGraph();
        Assert.assertEquals(graph.vertexSet().size(), 1);
        Assert.assertEquals(graph.vertexSet().iterator().next().getSequenceString(), seq);
    }

    private void assertSingleBubble(final TestAssembler assembler, final String one, final String two) {
        final SeqGraph graph = assembler.assemble();
        graph.simplifyGraph();
        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(graph);
        Assert.assertEquals(paths.size(), 2);
        final Set<String> expected = new HashSet<>(Arrays.asList(one, two));
        for ( final KBestHaplotype path : paths ) {
            final String seq = new String(path.bases());
            Assert.assertTrue(expected.contains(seq));
            expected.remove(seq);
        }
    }

    @Test(enabled = ! DEBUG)
    public void testRefCreation() {
        final String ref = "ACGTAACCGGTT";
        final TestAssembler assembler = new TestAssembler(3);
        assembler.addSequence(ref.getBytes(), true);
        assertLinearGraph(assembler, ref);
    }

    @Test(enabled = ! DEBUG)
    public void testRefNonUniqueCreation() {
        final String ref = "GAAAAT";
        final TestAssembler assembler = new TestAssembler(3);
        assembler.addSequence(ref.getBytes(), true);
        assertLinearGraph(assembler, ref);
    }

    @Test(enabled = ! DEBUG)
    public void testRefAltCreation() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref = "ACAACTGA";
        final String alt = "ACAGCTGA";
        assembler.addSequence(ref.getBytes(), true);
        assembler.addSequence(alt.getBytes(), false);
        assertSingleBubble(assembler, ref, alt);
    }

    @Test(enabled = ! DEBUG)
    public void testPartialReadsCreation() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref  = "ACAACTGA";
        final String alt1 = "ACAGCT";
        final String alt2 =    "GCTGA";
        assembler.addSequence(ref.getBytes(), true);
        assembler.addSequence(alt1.getBytes(), false);
        assembler.addSequence(alt2.getBytes(), false);
        assertSingleBubble(assembler, ref, "ACAGCTGA");
    }

    @Test(enabled = ! DEBUG)
    public void testMismatchInFirstKmer() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref = "ACAACTGA";
        final String alt =   "AGCTGA";
        assembler.addSequence(ref.getBytes(), true);
        assembler.addSequence(alt.getBytes(), false);

        final SeqGraph graph = assembler.assemble();
        graph.simplifyGraph();
        graph.removeSingletonOrphanVertices();
        final Set<SeqVertex> sources = graph.getSources();
        final Set<SeqVertex> sinks = graph.getSinks();

        Assert.assertEquals(sources.size(), 1);
        Assert.assertEquals(sinks.size(), 1);
        Assert.assertNotNull(graph.getReferenceSourceVertex());
        Assert.assertNotNull(graph.getReferenceSinkVertex());

        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(graph);
        Assert.assertEquals(paths.size(), 1);
    }

    @Test(enabled = ! DEBUG)
    public void testStartInMiddle() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref  = "CAAAATG";
        final String read =   "AAATG";
        assembler.addSequence(ref.getBytes(), true);
        assembler.addSequence(read.getBytes(), false);
        assertLinearGraph(assembler, ref);
    }

    @Test(enabled = ! DEBUG)
    public void testStartInMiddleWithBubble() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref  = "CAAAATGGGG";
        final String read =   "AAATCGGG";
        assembler.addSequence(ref.getBytes(), true);
        assembler.addSequence(read.getBytes(), false);
        assertSingleBubble(assembler, ref, "CAAAATCGGG");
    }

    @Test(enabled = ! DEBUG)
    public void testNoGoodStarts() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref  = "CAAAATGGGG";
        final String read =   "AAATCGGG";
        assembler.addSequence(ref.getBytes(), true);
        assembler.addSequence(read.getBytes(), false);
        assertSingleBubble(assembler, ref, "CAAAATCGGG");
    }


    @Test(enabled = !DEBUG)
    public void testCreateWithBasesBeforeRefSource() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref  =  "ACTG";
        final String read =   "CTGGGACT";
        assembler.addSequence(ReadThreadingGraphUnitTest.getBytes(ref), true);
        assembler.addSequence(ReadThreadingGraphUnitTest.getBytes(read), false);
        assertLinearGraph(assembler, "ACTGGGACT");
    }

    @Test(enabled = !DEBUG)
    public void testSingleIndelAsDoubleIndel3Reads() {
        final TestAssembler assembler = new TestAssembler(25);
        // The single indel spans two repetitive structures
        final String ref   = "GTTTTTCCTAGGCAAATGGTTTCTATAAAATTATGTGTGTGTGTCTCTCTCTGTGTGTGTGTGTGTGTGTGTGTGTATACCTAATCTCACACTCTTTTTTCTGG";
        final String read1 = "GTTTTTCCTAGGCAAATGGTTTCTATAAAATTATGTGTGTGTGTCTCT----------GTGTGTGTGTGTGTGTGTATACCTAATCTCACACTCTTTTTTCTGG";
        final String read2 = "GTTTTTCCTAGGCAAATGGTTTCTATAAAATTATGTGTGTGTGTCTCT----------GTGTGTGTGTGTGTGTGTATACCTAATCTCACACTCTTTTTTCTGG";
        assembler.addSequence(ReadThreadingGraphUnitTest.getBytes(ref), true);
        assembler.addSequence(ReadThreadingGraphUnitTest.getBytes(read1), false);
        assembler.addSequence(ReadThreadingGraphUnitTest.getBytes(read2), false);

        final SeqGraph graph = assembler.assemble();
        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(graph);
        Assert.assertEquals(paths.size(), 2);
        final byte[] refPath = paths.get(0).bases().length == ref.length() ? paths.get(0).bases() : paths.get(1).bases();
        final byte[] altPath = paths.get(0).bases().length == ref.length() ? paths.get(1).bases() : paths.get(0).bases();
        Assert.assertEquals(refPath, ReadThreadingGraphUnitTest.getBytes(ref));
        Assert.assertEquals(altPath, ReadThreadingGraphUnitTest.getBytes(read1));
    }
}
