package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.appengine.repackaged.com.google.common.collect.Lists;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;
import java.util.UUID;

public final class SelfKeyReadsTest {

    @DataProvider(name = "keyedVariantShardsReads")
    public Object[][] keyedReads(){
        DataflowTestData testData = new DataflowTestData();

        List<MutableGATKRead> reads = testData.getReads();
        List<KV<UUID, MutableGATKRead>> expected = Lists.newArrayList(
                KV.of(new UUID(1, 1), reads.get(0)),
                KV.of(new UUID(2, 2), reads.get(1)),
                KV.of(new UUID(3, 3), reads.get(2)),
                KV.of(new UUID(4, 4), reads.get(3))
        );
        return new Object[][]{
                {reads, expected},
        };
    }

    @Test(dataProvider = "keyedVariantShardsReads")
    public void selfKeyReadsTest(List<MutableGATKRead> reads, List<KV<UUID, MutableGATKRead>> expected) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<MutableGATKRead> pReads = DataflowTestUtils.PCollectionCreateAndVerify(p, reads);

        PCollection<KV<UUID, MutableGATKRead>> kReads = pReads.apply(new SelfKeyReads());
        DataflowAssert.that(kReads).containsInAnyOrder(expected);
        p.run();
    }
}