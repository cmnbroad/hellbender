package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import htsjdk.samtools.Cigar;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.Path;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class PathUnitTest extends BaseTest {
    @Test(enabled = true)
    public void testAlignReallyLongDeletion() {
        final String ref = "ATGGTGGCTCATACCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGAACATCACCTGAGGCCAGGAGTTCAAAACCAGCCTGGCTAACATAGCAAAACCCCATCTCTAATGAAAATACAAAAATTAGCTGGGTGTGGTGGTGTCCGCCTGTAGTCCCAGCTACTCAGGAGACTAAGGCATGAGAATCACTTGAACCCAGGATGCAGAGGCTGTAGTGAGCCGAGATTGCACCACGGCTGCACTCCAGCCTGGGCAACAGAGCGAGACTCTGTCTCAAATAAAATAGCGTAACGTAACATAACATAACATAACATAACATAACATAACATAACATAACATAACATAACATAACACAACAACAAAATAAAATAACATAAATCATGTTGTTAGGAAAAAAATCAGTTATGCAGCTACATGCTATTTACAAGAGATATACCTTAAAATATAAGACACAGAGGCCGGGCGCGGTAGCTCATGCCTGTAATCCCAGCACTTTGGGAGGCTGAGGCAAGCGGATCATGAGGTCAGGAGATCGAGACCATCC";
        final String hap = "ATGGTGGCTCATACCTGTAATCCCAGCACTTTGGGAGGCTGAGGCAAGCGGATCATGAGGTCAGGAGATCGAGACCATCCT";

        final SeqGraph graph = new SeqGraph(11);
        final SeqVertex v = new SeqVertex(hap);
        graph.addVertex(v);
        final Path<SeqVertex,BaseEdge> path = new Path<>(v, graph);
        final Cigar cigar = path.calculateCigar(ref.getBytes());
        Assert.assertNull(cigar, "Should have failed gracefully");
    }
}
