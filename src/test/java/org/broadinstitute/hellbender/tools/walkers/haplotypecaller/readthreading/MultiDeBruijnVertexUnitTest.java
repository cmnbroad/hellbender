package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import org.junit.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

public final class MultiDeBruijnVertexUnitTest {

   @Test
    public void test(){
       final MultiDeBruijnVertex v1 = new MultiDeBruijnVertex("fred".getBytes());
       final MultiDeBruijnVertex v2 = new MultiDeBruijnVertex("fred".getBytes());
       Assert.assertNotEquals(v1, v2);

       Assert.assertEquals(v1.getKmerSize(), 4);
       Assert.assertEquals(v1.getSuffix(), (byte)'d');
       Assert.assertTrue(Arrays.equals(v1.getSequence(), "fred".getBytes()));
       Assert.assertEquals(v1.getSuffixString(), "d");
       Assert.assertTrue(Arrays.equals(v1.getAdditionalSequence(true), "fred".getBytes()));
       Assert.assertTrue(Arrays.equals(v1.getAdditionalSequence(false), "d".getBytes()));

       Assert.assertNotEquals(v1.getId(), v2.getId());
       Assert.assertNotEquals(v1.hashCode(), v2.hashCode());
       Assert.assertEquals(v1.getKmerSize(), v2.getKmerSize());
       Assert.assertEquals(v1.additionalInfo(), v2.additionalInfo());
       Assert.assertTrue(Arrays.equals(v1.getSequence(), v2.getSequence()));
       Assert.assertEquals(v1.hasAmbiguousSequence(), v2.hasAmbiguousSequence());
       v1.toString();//not blow up - we dont check the string
   }


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullRead(){
        final MultiDeBruijnVertex v1 = new MultiDeBruijnVertex("fred".getBytes());
        v1.addRead(null);
    }

}
