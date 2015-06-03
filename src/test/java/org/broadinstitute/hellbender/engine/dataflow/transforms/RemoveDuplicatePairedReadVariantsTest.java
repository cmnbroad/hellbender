package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.appengine.repackaged.com.google.common.collect.Lists;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;
import org.broadinstitute.hellbender.utils.variant.SkeletonVariant;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class RemoveDuplicatePairedReadVariantsTest {

    @DataProvider(name = "dupedPairedReadsAndVariants")
    public Object[][] dupedPairedReadsAndVariants() {
        DataflowTestData testData = new DataflowTestData();
        List<KV<MutableGATKRead, Variant>> dupes = testData.getKvReadVariant();

        List<KV<UUID, UUID>> kvUUIDUUID = Lists.newArrayList(
                KV.of(dupes.get(0).getKey().getUUID(), dupes.get(0).getValue().getUUID()),
                KV.of(dupes.get(1).getKey().getUUID(), dupes.get(1).getValue().getUUID()),
                KV.of(dupes.get(2).getKey().getUUID(), dupes.get(2).getValue().getUUID()),
                KV.of(dupes.get(3).getKey().getUUID(), dupes.get(3).getValue().getUUID()),
                KV.of(dupes.get(4).getKey().getUUID(), dupes.get(4).getValue().getUUID()));

        Iterable<UUID> uuids0 = Lists.newArrayList(dupes.get(1).getValue().getUUID(), dupes.get(0).getValue().getUUID());
        Iterable<UUID> uuids2 = Lists.newArrayList(dupes.get(2).getValue().getUUID());
        Iterable<UUID> uuids3 = Lists.newArrayList(dupes.get(3).getValue().getUUID());
        ArrayList<KV<UUID, Iterable<UUID>>> kvUUIDiUUID = Lists.newArrayList(
                KV.of(dupes.get(0).getKey().getUUID(), uuids0),
                KV.of(dupes.get(2).getKey().getUUID(), uuids2),
                KV.of(dupes.get(3).getKey().getUUID(), uuids3)
        );



        Iterable<Variant> variant01 = Lists.newArrayList(dupes.get(1).getValue(), dupes.get(0).getValue());
        Iterable<Variant> variant2 = Lists.newArrayList(dupes.get(2).getValue());
        Iterable<Variant> variant3 = Lists.newArrayList(dupes.get(3).getValue());
        List<KV<UUID, Iterable<Variant>>> kvUUIDiVariant = Lists.newArrayList(
                KV.of(dupes.get(0).getKey().getUUID(), variant01),
                KV.of(dupes.get(2).getKey().getUUID(), variant2),
                KV.of(dupes.get(3).getKey().getUUID(), variant3));

        List<KV<MutableGATKRead, Iterable<Variant>>> finalExpected = testData.getKvReadiVariant();

        return new Object[][]{
                {dupes, kvUUIDUUID, kvUUIDiUUID, kvUUIDiVariant, finalExpected},
        };
    }

    @Test(dataProvider = "dupedPairedReadsAndVariants")
    public void removeDupedReadVariantsTest(List<KV<MutableGATKRead, Variant>> dupes, List<KV<UUID, UUID>> kvUUIDUUID,
                                            ArrayList<KV<UUID, Iterable<UUID>>> kvUUIDiUUID, List<KV<UUID, Iterable<Variant>>> kvUUIDiVariant,
                                            List<KV<MutableGATKRead, Iterable<Variant>>> finalExpected) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<KV<MutableGATKRead, Variant>> pKVs = DataflowTestUtils.PCollectionCreateAndVerify(p, dupes);

        PCollection<KV<UUID, UUID>> uuids = pKVs.apply(ParDo.of(new stripToUUIDs()));
        DataflowAssert.that(uuids).containsInAnyOrder(kvUUIDUUID);

        PCollection<KV<UUID, Iterable<UUID>>> readsIterableVariants = uuids.apply(RemoveDuplicates.<KV<UUID, UUID>>create()).apply(GroupByKey.<UUID, UUID>create());
        DataflowAssert.that(readsIterableVariants).containsInAnyOrder(kvUUIDiUUID);

        // Now add the variants back in.
        // Now join the stuff we care about back in.
        // Start by making KV<UUID, KV<UUID, Variant>> and joining that to KV<UUID, Iterable<UUID>>
        PCollection<KV<UUID, Iterable<Variant>>> matchedVariants = RemoveReadVariantDupesUtility.addBackVariants(pKVs, readsIterableVariants);
        DataflowAssert.that(matchedVariants).containsInAnyOrder(kvUUIDiVariant);

        PCollection<KV<MutableGATKRead, Iterable<Variant>>> finalResult = RemoveReadVariantDupesUtility.addBackReads(pKVs, matchedVariants);
        DataflowAssert.that(finalResult).containsInAnyOrder(finalExpected);

        p.run();
    }

    @Test(dataProvider = "dupedPairedReadsAndVariants")
    public void fullRemoveDupesTest(List<KV<MutableGATKRead, Variant>> dupes, List<KV<UUID, UUID>> kvUUIDUUID,
                                            ArrayList<KV<UUID, Iterable<UUID>>> kvUUIDiUUID, List<KV<UUID, Iterable<Variant>>> kvUUIDiVariant,
                                            List<KV<MutableGATKRead, Iterable<Variant>>> finalExpected) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<KV<MutableGATKRead, Variant>> pKVs = DataflowTestUtils.PCollectionCreateAndVerify(p, dupes);

        PCollection<KV<MutableGATKRead, Iterable<Variant>>> result = pKVs.apply(new RemoveDuplicatePairedReadVariants());
        DataflowAssert.that(result).containsInAnyOrder(finalExpected);

        p.run();
    }
}