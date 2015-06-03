package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.values.KV;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.engine.dataflow.datasources.FakeReferenceSource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.SkeletonVariant;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.util.List;
import java.util.UUID;

public class DataflowTestData {
    private final List<KV<Integer, Integer>> readStartLength;

    private final List<MutableGATKRead> reads;
    private final List<KV<ReferenceShard, Iterable<MutableGATKRead>>> kvRefShardiReads;
    private final List<SimpleInterval> readIntervals;
    private final List<SimpleInterval> allIntervals;
    private final List<KV<ReferenceBases, Iterable<MutableGATKRead>>> kvRefBasesiReads;
    private final List<KV<VariantShard, MutableGATKRead>> kvVariantShardRead;
    private final List<Variant> variants;
    private final List<KV<VariantShard, Variant>> kvVariantShardVariant;
    private final List<KV<MutableGATKRead, ReferenceBases>> kvReadsRefBases;
    private final List<KV<MutableGATKRead, Variant>> kvReadVariant;
    private final List<KV<MutableGATKRead, Iterable<Variant>>> kvReadiVariant;
    private final List<KV<MutableGATKRead, ReadContextData>> kvReadContextData;


    public DataflowTestData() {
        readStartLength = Lists.newArrayList(KV.of(100, 50), KV.of(140, 100), KV.of(100000, 10), KV.of(299999, 10));

        // TODO Make reads construction more general.
        reads = Lists.newArrayList(
                makeRead(readStartLength.get(0), 1),
                makeRead(readStartLength.get(1), 2),
                makeRead(readStartLength.get(2), 3),
                makeRead(readStartLength.get(3), 4));

        kvRefShardiReads =  Lists.newArrayList(
                KV.of(new ReferenceShard(0, "1"), Lists.newArrayList(reads.get(1), reads.get(0))),
                KV.of(new ReferenceShard(1, "1"), Lists.newArrayList(reads.get(2))),
                KV.of(new ReferenceShard(2, "1"), Lists.newArrayList(reads.get(3))));

        readIntervals = Lists.newArrayList(
                makeInterval(readStartLength.get(0)),
                makeInterval(readStartLength.get(1)),
                makeInterval(readStartLength.get(2)),
                makeInterval(readStartLength.get(3)));

        // The first two reads are mapped onto the same reference shard. The ReferenceBases returned should
        // be from the start of the first read [rStartLength.get(0).getKey()] to the end
        // the second [rStartLength.get(1).getKey() + rStartLength.get(1).getValue()-1].
        SimpleInterval spannedReadInterval =
                new SimpleInterval("1", readStartLength.get(0).getKey(), readStartLength.get(1).getKey() + readStartLength.get(1).getValue()-1);

        allIntervals = Lists.newArrayList(readIntervals.iterator());
        allIntervals.add(spannedReadInterval);

        kvRefBasesiReads = Lists.newArrayList(
                        KV.of(FakeReferenceSource.bases(spannedReadInterval), Lists.newArrayList(reads.get(1), reads.get(0))),
                        KV.of(FakeReferenceSource.bases(readIntervals.get(2)), Lists.newArrayList(reads.get(2))),
                        KV.of(FakeReferenceSource.bases(readIntervals.get(3)), Lists.newArrayList(reads.get(3))));

        kvReadsRefBases = Lists.newArrayList(
                KV.of(reads.get(0), getBases(reads.get(0).getStart(), reads.get(0).getEnd())),
                KV.of(reads.get(1), getBases(reads.get(1).getStart(), reads.get(1).getEnd())),
                KV.of(reads.get(2), getBases(reads.get(2).getStart(), reads.get(2).getEnd())),
                KV.of(reads.get(3), getBases(reads.get(3).getStart(), reads.get(3).getEnd()))
        );

        variants = Lists.newArrayList(
                new SkeletonVariant(new SimpleInterval("1", 170, 180), true, false, new UUID(1001, 1001)),
                new SkeletonVariant(new SimpleInterval("1", 210, 220), false, true, new UUID(1002, 1002)),
                new SkeletonVariant(new SimpleInterval("1", 100000, 100000), true, false, new UUID(1003, 1003)),
                new SkeletonVariant(new SimpleInterval("1", 299998, 300002), false, true, new UUID(1004, 1004))
        );

        kvVariantShardRead = Lists.newArrayList(
                KV.of(new VariantShard(0, "1"), reads.get(0)),
                KV.of(new VariantShard(0, "1"), reads.get(1)),
                KV.of(new VariantShard(1, "1"), reads.get(2)),
                KV.of(new VariantShard(2, "1"), reads.get(3)),    // The last read spans
                KV.of(new VariantShard(3, "1"), reads.get(3)));   // two shards.

        kvVariantShardVariant = Lists.newArrayList(
                KV.of(new VariantShard(0, "1"), variants.get(0)),
                KV.of(new VariantShard(0, "1"), variants.get(1)),
                KV.of(new VariantShard(1, "1"), variants.get(2)),
                KV.of(new VariantShard(2, "1"), variants.get(3)),    // The last variant spans
                KV.of(new VariantShard(3, "1"), variants.get(3)));   // two shards.

        kvReadVariant = Lists.newArrayList(
                KV.of(reads.get(1), variants.get(0)),
                KV.of(reads.get(1), variants.get(1)),
                KV.of(reads.get(2), variants.get(2)),
                KV.of(reads.get(3), variants.get(3)),    // The read and variant span two variant shards, that's
                KV.of(reads.get(3), variants.get(3))     // why there are two of them (2,3).
        );

        Iterable<Variant> variant01 = Lists.newArrayList(kvReadVariant.get(1).getValue(), kvReadVariant.get(0).getValue());
        Iterable<Variant> variant2 = Lists.newArrayList(kvReadVariant.get(2).getValue());
        Iterable<Variant> variant3 = Lists.newArrayList(kvReadVariant.get(3).getValue());

        kvReadiVariant = Lists.newArrayList(
                KV.of(kvReadVariant.get(0).getKey(), variant01),
                KV.of(kvReadVariant.get(2).getKey(), variant2),
                KV.of(kvReadVariant.get(3).getKey(), variant3)
        );

        kvReadContextData = Lists.newArrayList(
                KV.of(kvReadsRefBases.get(0).getKey(), new ReadContextData(kvReadsRefBases.get(0).getValue(), com.google.appengine.repackaged.com.google.common.collect.Lists.newArrayList())),
                KV.of(kvReadsRefBases.get(1).getKey(), new ReadContextData(kvReadsRefBases.get(1).getValue(), kvReadiVariant.get(0).getValue())),
                KV.of(kvReadsRefBases.get(2).getKey(), new ReadContextData(kvReadsRefBases.get(2).getValue(), kvReadiVariant.get(1).getValue())),
                KV.of(kvReadsRefBases.get(3).getKey(), new ReadContextData(kvReadsRefBases.get(3).getValue(), kvReadiVariant.get(2).getValue()))
        );
    }

    public MutableGATKRead makeRead(KV<Integer, Integer> startLength, int i) {
        return makeRead(startLength.getKey(), startLength.getValue(), i);
    }

    public MutableGATKRead makeRead(int start, int length, int i) {
        return ArtificialReadUtils.createRandomRead(start, length, i);
    }

    private SimpleInterval makeInterval(KV<Integer, Integer> startLength) {
        return new SimpleInterval("1", startLength.getKey(), startLength.getKey() + startLength.getValue() - 1);
    }

    private ReferenceBases getBases(int start, int end) {
        return FakeReferenceSource.bases(new SimpleInterval("1", start, end));
    }

    public final List<KV<Integer, Integer>> getReadStartLength() {
        return readStartLength;
    }

    public List<KV<ReferenceShard, Iterable<MutableGATKRead>>> getKvRefShardiReads() {
        return kvRefShardiReads;
    }

    public List<SimpleInterval> getReadIntervals() {
        return readIntervals;
    }

    public List<SimpleInterval> getAllIntervals() {
        return allIntervals;
    }

    public List<KV<ReferenceBases, Iterable<MutableGATKRead>>> getKvRefBasesiReads() {
        return kvRefBasesiReads;
    }

    public List<MutableGATKRead> getReads() {
        return reads;
    }

    public List<KV<MutableGATKRead, ReferenceBases>> getKvReadsRefBases() {
        return kvReadsRefBases;
    }

    public List<KV<MutableGATKRead, Iterable<Variant>>> getKvReadiVariant() {
        return kvReadiVariant;
    }

    public List<KV<MutableGATKRead, Variant>> getKvReadVariant() {
        return kvReadVariant;
    }


    public List<Variant> getVariants() {
        return variants;
    }

    public List<KV<MutableGATKRead, ReadContextData>> getKvReadContextData() {
        return kvReadContextData;
    }

    public List<KV<VariantShard, MutableGATKRead>> getKvVariantShardRead() {
        return kvVariantShardRead;
    }

    public List<KV<VariantShard, Variant>> getKvVariantShardVariant() {
        return kvVariantShardVariant;
    }
}
