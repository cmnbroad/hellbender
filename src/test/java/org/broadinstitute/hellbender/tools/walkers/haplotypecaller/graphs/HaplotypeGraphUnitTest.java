package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.HaplotypeGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Tests for {@link HaplotypeGraph}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class HaplotypeGraphUnitTest extends BaseTest {

   @Test(dataProvider="buildByStringDataProvider")
   public void testBuildByString(final String string, final int kmerSize, final int vertexCount, final int edgeCount) {
       final HaplotypeGraph g = new HaplotypeGraph(string);
       Assert.assertEquals(g.getKmerSize(), kmerSize, g.toString());
       Assert.assertEquals(g.vertexSet().size(), vertexCount, g.toString());
       Assert.assertEquals(g.edgeSet().size(), edgeCount, g.toString());
   }

   @Test(dataProvider="equalTestDataProvider")
   public void testEquals(final HaplotypeGraph one, final HaplotypeGraph two, final boolean expected) {
      Assert.assertEquals(one.equals(two), expected);
   }

   @Test(dataProvider="equalTestDataProvider")
   public void testEqualReciprocal(final HaplotypeGraph one, final HaplotypeGraph two, final boolean expected) {
       Assert.assertEquals(two.equals(one), expected);
   }


   @Test(dataProvider="equalTestDataProvider")
   public void testReflexibeEquals(final HaplotypeGraph one, final HaplotypeGraph two,
                                   @SuppressWarnings("unused") final boolean expected) {
       Assert.assertTrue(one.equals(one));
       Assert.assertTrue(two.equals(two));
   }

   @Test(dataProvider="mergingCommonChainsDataProvider")
   public void testMergingCommonChains(final HaplotypeGraph actual, HaplotypeGraph expected) {


        final Map<Kmer,MultiDeBruijnVertex> beforeMap = new HashMap<>(actual.uniqueKmerMap());
        actual.mergeCommonChains();
        final Map<Kmer,MultiDeBruijnVertex> afterMap = new HashMap<>(actual.uniqueKmerMap());
        final Map<Kmer,MultiDeBruijnVertex> mergedMap = new HashMap<>(expected.uniqueKmerMap());

        Assert.assertEquals(actual, expected, "" + actual.vertexSet() + " EDGES " + actual.edgeSet());
        Assert.assertEquals(beforeMap.size(), afterMap.size());
        Assert.assertEquals(afterMap.size(), mergedMap.size());
        for (final Kmer k : beforeMap.keySet()) {
            Assert.assertTrue(afterMap.containsKey(k));
            Assert.assertTrue(mergedMap.containsKey(k));
            final byte[] seq1 = beforeMap.get(k).getSequence();
            final byte[] seq2 = afterMap.get(k).getSequence();
            final byte[] seq3 = mergedMap.get(k).getSequence();
            Assert.assertEquals(seq1.length, seq2.length);
            Assert.assertEquals(seq2.length, seq3.length);
            for (int i = 0; i < seq3.length; i++) {
                final byte bk = k.base(i);
                final byte b1 = seq1[i];
                final byte b2 = seq2[i];
                final byte b3 = seq3[i];
                final byte theByte = b1 == 'N' || b2 == 'N' || b3 == 'N' ? (byte)'N' : b1;
                if (theByte == 'N') continue;
                Assert.assertEquals(b1, b2);
                Assert.assertEquals(b2, b3);
                Assert.assertEquals(bk, b1);
            }
        }
   }


   @DataProvider(name="mergingCommonChainsDataProvider")
   public Iterator<Object[]> mergingCommonChainsDataProvider() {
      final List<Object[]> list = new LinkedList<>();
      for (int i = 0; i < MERGING_COMMON_CHAINS_DATA.length; i += 2) {
          final HaplotypeGraph before = new HaplotypeGraph(MERGING_COMMON_CHAINS_DATA[i]);
          final HaplotypeGraph after = new HaplotypeGraph(MERGING_COMMON_CHAINS_DATA[i+1]);
          list.add(new Object[] { before , after});
      }
      return list.iterator();
   }

   @DataProvider(name="equalTestDataProvider")
   public Iterator<Object[]> equalsTestDataProvider() {
      final List<Object[]> result = new LinkedList<>();
      for (int i = 0; i < EQUAL_TEST_DATA.length; i += 3) {
          final HaplotypeGraph g1 = new HaplotypeGraph(EQUAL_TEST_DATA[i]);
          final HaplotypeGraph g2 = new HaplotypeGraph(EQUAL_TEST_DATA[i+1]);
          final boolean outcome = Boolean.parseBoolean(EQUAL_TEST_DATA[i + 2]);
          result.add(new Object[] { g1, g2, outcome});
      }
      return result.iterator();
   }

   @DataProvider(name="buildByStringDataProvider")
   public Iterator<Object[]> buildByStringDataProvider() {
      return Arrays.asList(BUILD_BY_STRING_TEST_DATA).iterator();
   }

   private static final Object[][] BUILD_BY_STRING_TEST_DATA = new Object[][] {
           {"[ks=3]{REF: ACT}",3,1,0},
           {"[ks=3]{REF: ACT(3) -> T(1)      ->      G(2) -> A}" +
                         "{ (3) -> A -> G ->          (2) }" +
                                "{  (1) -> A -> G ->  (2) }",3,8,9},
           {"[ks=3]{REF: ACT -> C(1) -> G}{ACT -> C(1) -> G}{ACT -> C(1) -> G}",3,5,4} ,
           {"[ks=3]{REF: ACT -> A(1) -> G -> A(2) -> C -> G -> T }" +
                              "{A(1) -> T -> A(2) }",3,8,8}  ,
           {"[ks=3]{REF: ACT -> A -> T(2) -> C -> A -> G -> T -> A -> C -> G -> T -> A(1) -> T}" +
                   "{ ACT -> A -> T(2) -> C -> A -> G -> T -> A -> C -> G -> T -> A(1) -> T}",3,15,14} ,
           {"[ks=3]{REF: ACT -> A -> T    -> C -> A -> G -> T -> A -> C -> G -> T -> A    -> T}",3,13,12},
           {"[ks=3]{REF: ACT -> A -> T(1) }" +
                   "{ ACT -> A -> T(1) }",3,5,4},
           {"[ks=3]{REF: TTT -> A(1) -> C -> T(2)}{ A(1) -> T(2) } ",3,4,4}
   };

   private static final String[] EQUAL_TEST_DATA = new String[] {
           "[ks=3]{REF: ACT}","[ks=3]{REF: ACT}", "true",
           "[ks=3]{REF: TCA}","[ks=3]{REF: ACT}", "false",
           "[ks=4]{REF: ACTG}","[ks=3]{REF: ACT}", "false",
           "[ks=3]{REF: ACT(3) -> T(1)      ->      G(2) -> A}" +
                        "{ (3) -> A -> G ->          (2) }" +
                               "{  (1) -> A -> G ->  (2) }"
                   ,"[ks=3]{REF: ACT(3) -> T(1)      ->      G(2) -> A}" +
                                                                 "{  (1) -> A -> G ->  (2) }" +
                                                          "{ (3) -> A -> G ->          (2) }", "true",
           "[ks=3]{REF: ACT(3) -> T(1)      ->      G(2) -> A}" +
                        "{ (3) -> A -> T ->          (2) }" +
                              "{  (1) -> A -> G ->  (2) }"
                    ,"[ks=3]{REF: ACT(3) -> T(1)      ->      G(2) -> A}" +
                     "{  (1) -> A -> G ->  (2) }" +
                     "{ (3) -> A -> G ->          (2) }", "true",
           "[ks=3]{REF: ACT -> G -> C(2) }{ ACT -> T -> C(2) }","[ks=3]{REF: ACT -> T -> C(2) }{ ACT -> G -> C(2) }","false",

   };

   private static final String[] MERGING_COMMON_CHAINS_DATA = new String[] {  // pairs before and after.
           "[ks=3]{REF: ACT -> A(1) -> G -> A -> G(2) -> T }" +
                             "{A(1) -> T -> A -> G(2) }",
           "[ks=3]{REF: ACT -> A(1) -> G -> A(2) -> G -> T }" +
                             "{A(1) -> T -> A(2) }",

           "[ks=3]{REF: ACT -> A(1) -> G -> A -> C -> G(2) -> T }" +
                             "{A(1) -> T -> A -> C -> G(2) }",
           "[ks=3]{REF: ACT -> A(1) -> G -> A(2) -> C -> G -> T }" +
                             "{A(1) -> T -> A(2) }",

           "[ks=3]{REF: ACT -> A -> T(1) -> C -> A -> G -> T -> A -> C(2) -> G -> T -> A}" +
                                 "{ T(1) -> A -> A -> G -> T -> A -> C(2) }",
           "[ks=3]{REF: ACT -> A -> T(1) -> C -> A(2) -> G -> T -> A -> C -> G -> T -> A}" +
                                 "{ T(1) -> A -> A(2) } ",

//           "[ks=3]{REF: ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A(1)}" +
//                     "{ ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A(1)}" ,
//           "[ks=3]{REF: ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A}" ,

           "[ks=3]{REF: ACT -> A -> T(1) }" +
                     "{ AGT -> A -> T(1) }" ,
           "[ks=3]{REF: ACT -> A(1) -> T }" +
                   "{ AGT -> A(1)  }"  ,
           "[ks=3]{REF: ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A -> T}" ,
           "[ks=3]{REF: ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A -> T}" ,
           "[ks=3]{REF: ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A -> T(1)}" + "{ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A -> T(1)}" ,
           "[ks=3]{REF: ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A -> T}" ,

           "[ks=3]{REF: TTT -> T -> T -> T -> T -> T -> T -> T -> T -> T -> T(1) -> T -> T -> T(2) -> T -> T}"
                  + "{  TTT -> T -> T -> T -> T -> T -> T -> T -> T -> T -> T(1) -> G -> T -> T -> T -> T -> T(2) -> T -> T}",
           "[ks=3]{REF: TTT -> T -> T -> T -> T -> T -> T -> T -> T -> T -> T(1) -> T(2) -> T -> T -> T -> T}"
                  +                                                      "{ T(1) -> G -> T -> T -> T(2) }",

           "[ks=3]{REF: TTT -> T -> G(1) -> A -> C -> C -> T(2)}" +
                                 "{ G(1) -> T -> C -> C -> T(2)}" +
                                 "{ G(1) -> G -> C -> C -> T(2)}" +
                                 "{ G(1) -> C -> T(2)} ",
           "[ks=3]{REF: TTT -> T -> G(1) -> A -> C(2) -> C(3) -> T }" +
                                 "{  G(1) -> T -> C(2) }" +
                                 "{  G(1) -> G -> C(2) }" +
                                 "{  G(1) -> C(3) }",

           "[ks=3]{REF: TTT -> T -> G(1) -> A -> C -> G}{ TTT -> T -> G(1) -> G -> C -> G}",
           "[ks=3]{REF: TTT -> T -> G(1) -> A -> C -> G}{ G(1) -> G -> C -> G}",

   };

}
