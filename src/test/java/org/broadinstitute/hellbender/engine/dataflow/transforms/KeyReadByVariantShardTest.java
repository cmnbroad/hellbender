package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.appengine.repackaged.com.google.common.collect.Lists;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantShard;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;
import java.util.UUID;

public final class KeyReadByVariantShardTest {

    @DataProvider(name = "keyedVariantShardsReads")
    public Object[][] keyedVariantShardsReads(){
        DataflowTestData testData = new DataflowTestData();

        List<MutableGATKRead> reads = testData.getReads();
        List<KV<VariantShard, MutableGATKRead>> kvVariantShardRead = testData.getKvVariantShardRead();

        return new Object[][]{
                {reads, kvVariantShardRead},
        };
    }

    @Test(dataProvider = "keyedVariantShardsReads")
    public void keyReadsByVariantShardTest(List<MutableGATKRead> reads, List<KV<VariantShard, MutableGATKRead>> kvVariantShardRead) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<MutableGATKRead> pRead = DataflowTestUtils.PCollectionCreateAndVerify(p, reads);

        PCollection<KV<VariantShard, MutableGATKRead>> kVariant = pRead.apply(new KeyReadByVariantShard());
        DataflowAssert.that(kVariant).containsInAnyOrder(kvVariantShardRead);
        p.run();
    }
}