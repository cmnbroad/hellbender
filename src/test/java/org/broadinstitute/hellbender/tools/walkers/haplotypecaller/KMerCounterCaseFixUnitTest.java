package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.HashSet;
import java.util.Set;

public class KMerCounterCaseFixUnitTest extends BaseTest {
    @Test
	public void testMyData() {
        final KMerCounter counter = new KMerCounter(3);

        Assert.assertNotNull(counter.toString());

        counter.addKmers(
			 "ATG", "ATG", "ATG", "ATG",
			 "ACC", "ACC", "ACC",
			 "AAA", "AAA",
			 "CTG",
			 "NNA",
                "CCC"
			 );

        testCounting(counter, "ATG", 4);
        testCounting(counter, "ACC", 3);
        testCounting(counter, "AAA", 2);
        testCounting(counter, "CTG", 1);
        testCounting(counter, "NNA", 1);
        testCounting(counter, "CCC", 1);
        testCounting(counter, "NNN", 0);
        testCounting(counter, "NNC", 0);

        Assert.assertNotNull(counter.toString());

        assertCounts(counter, 5);
        assertCounts(counter, 4, "ATG");
        assertCounts(counter, 3, "ATG", "ACC");
        assertCounts(counter, 2, "ATG", "ACC", "AAA");
        assertCounts(counter, 1, "ATG", "ACC", "AAA", "CTG", "NNA", "CCC");
    }

    private void assertCounts(final KMerCounter counter, final int minCount, final String... expecteds) {
        final Set<Kmer> expected = new HashSet<>();
        for ( final String one : expecteds ) expected.add(new Kmer(one));
        Assert.assertEquals(new HashSet<>(counter.getKmersWithCountsAtLeast(minCount)), expected);
    }

    private void testCounting(final KMerCounter counter, final String in, final int expectedCount) {
        Assert.assertEquals(counter.getKmerCount(new Kmer(in)), expectedCount);
    }
}
