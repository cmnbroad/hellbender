package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class DeBruijnVertexUnitTest extends BaseTest {
    @Test
    public void testBasic() {
        final byte[] bases = "ACT".getBytes();
        final DeBruijnVertex v = new DeBruijnVertex(bases);
        Assert.assertEquals(v.getSequence(), bases);
        Assert.assertEquals(v.getSequenceString(), new String(bases));
        Assert.assertEquals(v.length(), bases.length);
        Assert.assertEquals(v.getSuffix(), (byte) 'T');
        Assert.assertEquals(v.getSuffixString(), "T");

        Assert.assertEquals(v.getAdditionalSequence(true), bases);
        Assert.assertEquals(v.getAdditionalSequence(false).length, 1);
        Assert.assertEquals(v.getAdditionalSequence(false)[0], (byte) 'T');
    }
}
