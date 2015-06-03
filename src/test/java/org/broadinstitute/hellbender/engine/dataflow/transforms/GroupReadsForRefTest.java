package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class GroupReadsForRefTest {

    @DataProvider(name = "refShards")
    public Object[][] refShards(){

        DataflowTestData testData = new DataflowTestData();

        List<MutableGATKRead> inputs = testData.getReads();
        List<KV<ReferenceShard, Iterable<MutableGATKRead>>> kvs = testData.getKvRefShardiReads();

        return new Object[][]{
                {inputs, kvs},
        };
    }

    @Test(dataProvider = "refShards")
    public void groupReadsForRefTest(List<MutableGATKRead> reads, List<KV<ReferenceShard, Iterable<MutableGATKRead>>> expectedResult) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<MutableGATKRead> pReads = DataflowTestUtils.PCollectionCreateAndVerify(p, reads);
        PCollection<KV<ReferenceShard, Iterable<MutableGATKRead>>> grouped = pReads.apply(new GroupReadsForRef());

        DataflowAssert.that(grouped).containsInAnyOrder(expectedResult);
        p.run();
    }
}