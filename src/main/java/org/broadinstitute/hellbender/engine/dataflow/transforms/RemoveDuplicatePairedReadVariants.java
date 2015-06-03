package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.api.client.util.Sets;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.transforms.join.CoGbkResult;
import com.google.cloud.dataflow.sdk.transforms.join.CoGroupByKey;
import com.google.cloud.dataflow.sdk.transforms.join.KeyedPCollectionTuple;
import com.google.cloud.dataflow.sdk.values.*;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.util.List;
import java.util.Set;
import java.util.UUID;


public class RemoveDuplicatePairedReadVariants extends PTransform<PCollection<KV<MutableGATKRead,Variant>>, PCollection<KV<MutableGATKRead, Iterable<Variant>>>> {
    private static final long serialVersionUID = 1L;

    @Override
    public PCollection<KV<MutableGATKRead, Iterable<Variant>>> apply(PCollection<KV<MutableGATKRead, Variant>> input) {
        PCollection<KV<UUID, UUID>> readVariantUuids = input.apply(ParDo.of(new stripToUUIDs())).setName("RemoveDuplicatePairedReadVariants_stripToUUIDs");

        PCollection<KV<UUID, Iterable<UUID>>> readsIterableVariants = readVariantUuids.apply(RemoveDuplicates.<KV<UUID, UUID>>create()).apply(GroupByKey.<UUID, UUID>create());

        // Now join the stuff we care about back in.
        // Start by making KV<UUID, KV<UUID, Variant>> and joining that to KV<UUID, Iterable<UUID>>
        PCollection<KV<UUID, Iterable<Variant>>> matchedVariants =
                RemoveReadVariantDupesUtility.addBackVariants(input, readsIterableVariants);

        // And now, we do the same song and dance to get the Reads back in.
        return RemoveReadVariantDupesUtility.addBackReads(input, matchedVariants);
    }
}

class stripToUUIDs extends DoFn<KV<MutableGATKRead, Variant>, KV<UUID, UUID>> {
    @Override
    public void processElement(ProcessContext c) throws Exception {
        c.output(KV.of(c.element().getKey().getUUID(), c.element().getValue().getUUID()));
    }
}

class RemoveReadVariantDupesUtility {
    public static PCollection<KV<UUID, Iterable<Variant>>> addBackVariants(PCollection<KV<MutableGATKRead, Variant>> readVariants, PCollection<KV<UUID, Iterable<UUID>>> readsIterableVariants) {
        // Now join the stuff we care about back in.
        // Start by making KV<UUID, KV<UUID, Variant>> and joining that to KV<UUID, Iterable<UUID>>
        PCollection<KV<UUID, KV<UUID, Variant>>> variantToBeJoined =
                readVariants.apply(ParDo.of(new DoFn<KV<MutableGATKRead, Variant>, KV<UUID, KV<UUID, Variant>>>() {
                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        c.output(KV.of(c.element().getKey().getUUID(), KV.of(c.element().getValue().getUUID(), c.element().getValue())));
                    }
                })).setName("RemoveReadVariantDupesUtility_VariantToBeJoined");

        // Another coGroupBy...
        final TupleTag<Iterable<UUID>> iterableUuidTag = new TupleTag<>();
        final TupleTag<KV<UUID, Variant>> variantUuidTag = new TupleTag<>();

        PCollection<KV<UUID, CoGbkResult>> coGbkAgain = KeyedPCollectionTuple
                .of(iterableUuidTag, readsIterableVariants)
                .and(variantUuidTag, variantToBeJoined).apply(CoGroupByKey.<UUID>create());

        // For each Read UUID, get all of the variants that have a variant UUID
        return coGbkAgain.apply(ParDo.of(
                new DoFn<KV<UUID, CoGbkResult>, KV<UUID, Iterable<Variant>>>() {
                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        Iterable<KV<UUID,Variant>> kVariants = c.element().getValue().getAll(variantUuidTag);
                        Iterable<Iterable<UUID>> kUUIDss = c.element().getValue().getAll(iterableUuidTag);
                        // For every UUID that's left, keep the variant.
                        Set<UUID> uuidHashSet = Sets.newHashSet();
                        for (Iterable<UUID> uuids : kUUIDss) {
                            for (UUID uuid : uuids) {
                                uuidHashSet.add(uuid);
                            }
                        }
                        Set<Variant> iVariants = Sets.newHashSet();
                        for (KV<UUID,Variant> uVariants : kVariants) {
                            if (uuidHashSet.contains(uVariants.getKey())) {
                                iVariants.add(uVariants.getValue());
                            }
                        }
                        c.output(KV.of(c.element().getKey(), iVariants));
                    }
                })).setName("RemoveReadVariantDupesUtility_CoGroupBy");

    }

    public static PCollection<KV<MutableGATKRead, Iterable<Variant>>> addBackReads(PCollection<KV<MutableGATKRead, Variant>> readVariants, PCollection<KV<UUID, Iterable<Variant>>> matchedVariants) {
        // And now, we do the same song and dance to get the Reads back in.
        final TupleTag<MutableGATKRead> justReadTag = new TupleTag<>();
        final TupleTag<Iterable<Variant>> iterableVariant = new TupleTag<>();

        PCollection<KV<UUID, MutableGATKRead>> kReads = readVariants.apply(Keys.<MutableGATKRead>create()).apply(new SelfKeyReads());
        PCollection<KV<UUID, CoGbkResult>> coGbkLast = KeyedPCollectionTuple
                .of(justReadTag, kReads)
                .and(iterableVariant, matchedVariants).apply(CoGroupByKey.<UUID>create());

        return coGbkLast.apply(ParDo.of(new DoFn<KV<UUID, CoGbkResult>, KV<MutableGATKRead, Iterable<Variant>>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                Iterable<MutableGATKRead> iReads = c.element().getValue().getAll(justReadTag);
                // We only care about the first read (the rest are the same.
                Iterable<Iterable<Variant>> variants = c.element().getValue().getAll(iterableVariant);
                List<MutableGATKRead> reads = Lists.newArrayList();
                for (MutableGATKRead r : iReads) {
                    reads.add(r);
                }
                if (reads.size() < 1) {
                    throw new GATKException("no read found");
                }

                for (Iterable<Variant> v : variants) {
                        c.output(KV.of(reads.get(0), v));
                }
            }
        })).setName("RemoveDuplicatePairedReadVariants_addBackReads");
    }
}
