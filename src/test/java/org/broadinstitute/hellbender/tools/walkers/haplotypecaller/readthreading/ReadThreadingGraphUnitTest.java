package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.KBestHaplotypeFinder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class ReadThreadingGraphUnitTest extends BaseTest {
    private final static boolean DEBUG = false;

    public static byte[] getBytes(final String alignment) {
        return alignment.replace("-","").getBytes();
    }

    private void assertNonUniques(final ReadThreadingGraph assembler, String... nonUniques) {
        final Set<String> actual = new HashSet<>();
        assembler.buildGraphIfNecessary();
        for ( final Kmer kmer : assembler.getNonUniqueKmers() ) actual.add(kmer.baseString());
        final Set<String> expected = new HashSet<>(Arrays.asList(nonUniques));
        Assert.assertEquals(actual, expected);
    }

    @Test
    public void testSimpleHaplotypeRethreading() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(11);
        final String ref   = "CATGCACTTTAAAACTTGCCTTTTTAACAAGACTTCCAGATG";
        final String alt   = "CATGCACTTTAAAACTTGCCGTTTTAACAAGACTTCCAGATG";
        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(alt), false);
        assembler.buildGraphIfNecessary();
        Assert.assertNotEquals(ref.length() - 11 + 1, assembler.vertexSet().size(), "the number of vertex in the graph is the same as if there was no alternative sequence");
        Assert.assertEquals(ref.length() - 11 + 1 + 11, assembler.vertexSet().size(), "the number of vertex in the graph is not the same as if there is an alternative sequence");
        MultiDeBruijnVertex startAlt = assembler.findKmer(new Kmer(alt.getBytes(),20,11));
        Assert.assertNotNull(startAlt);
    }

    @Test(enabled = ! DEBUG)
    public void testNonUniqueMiddle() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(3);
        final String ref   = "GACACACAGTCA";
        final String read1 = "GACAC---GTCA";
        final String read2 =   "CAC---GTCA";
        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(read1), false);
        assembler.addSequence(getBytes(read2), false);
        assertNonUniques(assembler, "ACA", "CAC");
    }

    @Test(enabled = ! DEBUG)
    public void testReadsCreateNonUnique() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(3);
        final String ref   = "GCAC--GTCA"; // CAC is unique
        final String read1 = "GCACACGTCA"; // makes CAC non unique because it has a duplication
        final String read2 =    "CACGTCA"; // shouldn't be allowed to match CAC as start
        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(read1), false);
        assembler.addSequence(getBytes(read2), false);
//        assembler.convertToSequenceGraph().printGraph(new File("test.dot"), 0);

        assertNonUniques(assembler, "CAC");
        //assertSingleBubble(assembler, ref, "CAAAATCGGG");
    }

    @Test(enabled = ! DEBUG)
         public void testCountingOfStartEdges() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(3);
        final String ref   = "NNNGTCAAA"; // ref has some bases before start
        final String read1 =    "GTCAAA"; // starts at first non N base

        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(read1), false);
        assembler.buildGraphIfNecessary();
//        assembler.printGraph(new File("test.dot"), 0);

        for ( final MultiSampleEdge edge : assembler.edgeSet() ) {
            final MultiDeBruijnVertex source = assembler.getEdgeSource(edge);
            final MultiDeBruijnVertex target = assembler.getEdgeTarget(edge);
            final boolean headerVertex = source.getSuffix() == 'N' || target.getSuffix() == 'N';
            if ( headerVertex ) {
                Assert.assertEquals(edge.getMultiplicity(), 1, "Bases in the unique reference header should have multiplicity of 1");
            } else {
                Assert.assertEquals(edge.getMultiplicity(), 2, "Should have multiplicity of 2 for any edge outside the ref header but got " + edge + " " + source + " -> " + target);
            }
        }
    }

    @Test(enabled = !DEBUG)
    public void testCountingOfStartEdgesWithMultiplePrefixes() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(3);
        assembler.increaseCountsThroughBranches = true;
        final String ref   = "NNNGTCAXX"; // ref has some bases before start
        final String alt1  = "NNNCTCAXX"; // alt1 has SNP right after N
        final String read  =     "TCAXX"; // starts right after SNP, but merges right before branch

        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(alt1), false);
        assembler.addSequence(getBytes(read), false);
        assembler.buildGraphIfNecessary();
        assembler.printGraph(new File("test.dot"), 0);

        final List<String> oneCountVertices = Arrays.asList("NNN", "NNG", "NNC", "NGT", "NCT");
        final List<String> threeCountVertices = Arrays.asList("CAX", "AXX");

        for ( final MultiSampleEdge edge : assembler.edgeSet() ) {
            final MultiDeBruijnVertex source = assembler.getEdgeSource(edge);
            final MultiDeBruijnVertex target = assembler.getEdgeTarget(edge);
            final int expected = oneCountVertices.contains(target.getSequenceString()) ? 1 : (threeCountVertices.contains(target.getSequenceString()) ? 3 : 2);
            Assert.assertEquals(edge.getMultiplicity(), expected, "Bases at edge " + edge + " from " + source + " to " + target + " has bad multiplicity");
        }
    }

    @Test(enabled = !DEBUG)
    public void testCyclesInGraph() {

        // b37 20:12655200-12655850
        final String ref = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTATACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";

        // SNP at 20:12655528 creates a cycle for small kmers
        final String alt = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTGTACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";

        final List<SAMRecord> reads = new ArrayList<>();
        for ( int index = 0; index < alt.length() - 100; index += 20 )
            reads.add(ArtificialSAMUtils.createArtificialRead(Arrays.copyOfRange(alt.getBytes(), index, index + 100), Utils.dupBytes((byte) 30, 100), 100 + "M"));

        // test that there are cycles detected for small kmer
        final ReadThreadingGraph rtgraph25 = new ReadThreadingGraph(25);
        rtgraph25.addSequence("ref", ref.getBytes(), true);
        for ( final SAMRecord read : reads )
            rtgraph25.addRead(read);
        rtgraph25.buildGraphIfNecessary();
        Assert.assertTrue(rtgraph25.hasCycles());

        // test that there are no cycles detected for large kmer
        final ReadThreadingGraph rtgraph75 = new ReadThreadingGraph(75);
        rtgraph75.addSequence("ref", ref.getBytes(), true);
        for ( final SAMRecord read : reads )
            rtgraph75.addRead(read);
        rtgraph75.buildGraphIfNecessary();
        Assert.assertFalse(rtgraph75.hasCycles());
    }

    @Test(enabled = !DEBUG)
    public void testNsInReadsAreNotUsedForGraph() {

        final int length = 100;
        final byte[] ref = Utils.dupBytes((byte) 'A', length);

        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(25);
        rtgraph.addSequence("ref", ref, true);

        // add reads with Ns at any position
        for ( int i = 0; i < length; i++ ) {
            final byte[] bases = ref.clone();
            bases[i] = 'N';
            final SAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, Utils.dupBytes((byte) 30, length), length + "M");
            rtgraph.addRead(read);
        }
        rtgraph.buildGraphIfNecessary();

        final SeqGraph graph = rtgraph.convertToSequenceGraph();
        Assert.assertEquals(new KBestHaplotypeFinder(graph, graph.getReferenceSourceVertex(), graph.getReferenceSinkVertex()).size(), 1);
    }

// TODO -- update to use determineKmerSizeAndNonUniques directly
//    @DataProvider(name = "KmerSizeData")
//    public Object[][] makeKmerSizeDataProvider() {
//        List<Object[]> tests = new ArrayList<Object[]>();
//
//        // this functionality can be adapted to provide input data for whatever you might want in your data
//        tests.add(new Object[]{3, 3, 3, Arrays.asList("ACG"), Arrays.asList()});
//        tests.add(new Object[]{3, 4, 3, Arrays.asList("CAGACG"), Arrays.asList()});
//
//        tests.add(new Object[]{3, 3, 3, Arrays.asList("AAAAC"), Arrays.asList("AAA")});
//        tests.add(new Object[]{3, 4, 4, Arrays.asList("AAAAC"), Arrays.asList()});
//        tests.add(new Object[]{3, 5, 4, Arrays.asList("AAAAC"), Arrays.asList()});
//        tests.add(new Object[]{3, 4, 3, Arrays.asList("CAAA"), Arrays.asList()});
//        tests.add(new Object[]{3, 4, 4, Arrays.asList("CAAAA"), Arrays.asList()});
//        tests.add(new Object[]{3, 5, 4, Arrays.asList("CAAAA"), Arrays.asList()});
//        tests.add(new Object[]{3, 5, 5, Arrays.asList("ACGAAAAACG"), Arrays.asList()});
//
//        for ( int maxSize = 3; maxSize < 20; maxSize++ ) {
//            for ( int dupSize = 3; dupSize < 20; dupSize++ ) {
//                final int expectedSize = Math.min(maxSize, dupSize);
//                final String dup = Utils.dupString("C", dupSize);
//                final List<String> nonUnique = dupSize > maxSize ? Arrays.asList(Utils.dupString("C", maxSize)) : Collections.<String>emptyList();
//                tests.add(new Object[]{3, maxSize, expectedSize, Arrays.asList("ACGT", "A" + dup + "GT"), nonUnique});
//                tests.add(new Object[]{3, maxSize, expectedSize, Arrays.asList("A" + dup + "GT", "ACGT"), nonUnique});
//            }
//        }
//
//        return tests.toArray(new Object[][]{});
//    }
//
//    /**
//     * Example testng test using MyDataProvider
//     */
//    @Test(dataProvider = "KmerSizeData")
//    public void testDynamicKmerSizing(final int min, final int max, final int expectKmer, final List<String> seqs, final List<String> expectedNonUniques) {
//        final ReadThreadingGraph assembler = new ReadThreadingGraph(min, max);
//        for ( String seq : seqs ) assembler.addSequence(seq.getBytes(), false);
//        assembler.buildGraphIfNecessary();
//        Assert.assertEquals(assembler.getKmerSize(), expectKmer);
//        assertNonUniques(assembler, expectedNonUniques.toArray(new String[]{}));
//    }

}
