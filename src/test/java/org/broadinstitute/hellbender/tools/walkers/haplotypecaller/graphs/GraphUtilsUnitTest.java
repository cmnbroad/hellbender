package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.collections.PrimitivePair;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import static org.broadinstitute.hellbender.utils.Utils.reverse;

public final class GraphUtilsUnitTest extends BaseTest {
    @DataProvider(name = "findLongestUniqueMatchData")
    public Object[][] makefindLongestUniqueMatchData() {
        List<Object[]> tests = new ArrayList<>();

        { // test all edge conditions
            final String ref = "ACGT";
            for ( int start = 0; start < ref.length(); start++ ) {
                for ( int end = start + 1; end <= ref.length(); end++ ) {
                    final String kmer = ref.substring(start, end);
                    tests.add(new Object[]{ref, kmer, end - 1, end - start});
                    tests.add(new Object[]{ref, "N" + kmer, end - 1, end - start});
                    tests.add(new Object[]{ref, "NN" + kmer, end - 1, end - start});
                    tests.add(new Object[]{ref, kmer + "N", -1, 0});
                    tests.add(new Object[]{ref, kmer + "NN", -1, 0});
                }
            }
        }

        { // multiple matches
            final String ref = "AACCGGTT";
            for ( final String alt : Arrays.asList("A", "C", "G", "T") )
                tests.add(new Object[]{ref, alt, -1, 0});
            tests.add(new Object[]{ref, "AA", 1, 2});
            tests.add(new Object[]{ref, "CC", 3, 2});
            tests.add(new Object[]{ref, "GG", 5, 2});
            tests.add(new Object[]{ref, "TT", 7, 2});
        }

        { // complex matches that have unique substrings of lots of parts of kmer in the ref
            final String ref = "ACGTACGTACGT";
            tests.add(new Object[]{ref, "ACGT", -1, 0});
            tests.add(new Object[]{ref, "TACGT", -1, 0});
            tests.add(new Object[]{ref, "GTACGT", -1, 0});
            tests.add(new Object[]{ref, "CGTACGT", -1, 0});
            tests.add(new Object[]{ref, "ACGTACGT", -1, 0});
            tests.add(new Object[]{ref, "TACGTACGT", 11, 9});
            tests.add(new Object[]{ref, "NTACGTACGT", 11, 9});
            tests.add(new Object[]{ref, "GTACGTACGT", 11, 10});
            tests.add(new Object[]{ref, "NGTACGTACGT", 11, 10});
            tests.add(new Object[]{ref, "CGTACGTACGT", 11, 11});
        }

        return tests.toArray(new Object[][]{});
    }

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "findLongestUniqueMatchData")
    public void testfindLongestUniqueMatch(final String seq, final String kmer, final int start, final int length) {
        // adaptor this code to do whatever testing you want given the arguments start and size
        final PrimitivePair.Int actual = GraphUtils.findLongestUniqueSuffixMatch(seq.getBytes(), kmer.getBytes());
        if ( start == -1 )
            Assert.assertNull(actual);
        else {
            Assert.assertNotNull(actual);
            Assert.assertEquals(actual.first, start);
            Assert.assertEquals(actual.second, length);
        }
    }

    @Test
    public void compPrefixLen1(){
        final List<byte[]> bytes = Arrays.asList("ABC".getBytes(), "CDE".getBytes());
        final int pref = GraphUtils.compPrefixLen(bytes, 3);
        Assert.assertEquals(pref, 0);
    }

    @Test
    public void compPrefixLen2(){
        final List<byte[]> bytes = Arrays.asList("ABC".getBytes(), "ABD".getBytes());
        final int pref = GraphUtils.compPrefixLen(bytes, 3);
        Assert.assertEquals(pref, 2);
    }

    @Test
    public void compPrefixLen3(){
        final List<byte[]> bytes = Arrays.asList("ABC".getBytes(), "ABD".getBytes());
        final int pref = GraphUtils.compPrefixLen(bytes, 1);
        Assert.assertEquals(pref, 1);
    }

    @Test
    public void compSuffixLen1(){
        final List<byte[]> bytes = Arrays.asList(reverse("ABC".getBytes()), reverse("CDE".getBytes()));
        final int pref = GraphUtils.compSuffixLen(bytes, 3);
        Assert.assertEquals(pref, 0);
    }

    @Test
    public void compSuffixLen2(){
        final List<byte[]> bytes = Arrays.asList(reverse("ABC".getBytes()), reverse("ABD".getBytes()));
        final int pref = GraphUtils.compSuffixLen(bytes, 3);
        Assert.assertEquals(pref, 2);
    }

    @Test
    public void compSuffixLen3(){
        final List<byte[]> bytes = Arrays.asList(reverse("ABC".getBytes()), reverse("ABD".getBytes()));
        final int pref = GraphUtils.compSuffixLen(bytes, 1);
        Assert.assertEquals(pref, 1);
    }

    @Test
    public void kMers(){
        final SeqVertex v1 = new SeqVertex("fred");
        final SeqVertex v2 = new SeqVertex("frodo");
        final Collection<SeqVertex> vertices = Arrays.asList(v1, v2);
        final List<byte[]> kmers = GraphUtils.getKmers(vertices);
        Assert.assertEquals(kmers.size(), 2);
        Assert.assertTrue(Arrays.equals(kmers.get(0), "fred".getBytes()));
        Assert.assertTrue(Arrays.equals(kmers.get(1), "frodo".getBytes()));
    }

    @Test
    public void minKmerLength(){
        final SeqVertex v1 = new SeqVertex("fred");
        final SeqVertex v2 = new SeqVertex("frodo");
        final Collection<SeqVertex> vertices = Arrays.asList(v1, v2);
        final List<byte[]> kmers = GraphUtils.getKmers(vertices);
        Assert.assertEquals(GraphUtils.minKmerLength(kmers), 4);
    }

}
