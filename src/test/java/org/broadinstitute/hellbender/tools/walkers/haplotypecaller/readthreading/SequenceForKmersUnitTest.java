package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class SequenceForKmersUnitTest extends BaseTest {
    @Test
    public void testNoCount() {
        final byte[] seq = "ACGT".getBytes();
        final SequenceForKmers sk = new SequenceForKmers("foo", seq, 0, seq.length, 1, true);
        Assert.assertEquals(sk.name, "foo");
        Assert.assertEquals(sk.sequence, seq);
        Assert.assertEquals(sk.start, 0);
        Assert.assertEquals(sk.stop, seq.length);
        Assert.assertEquals(sk.count, 1);
        Assert.assertEquals(sk.isRef, true);
    }
}
