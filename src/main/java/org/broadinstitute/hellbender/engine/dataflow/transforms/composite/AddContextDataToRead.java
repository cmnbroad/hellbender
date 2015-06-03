package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.join.CoGbkResult;
import com.google.cloud.dataflow.sdk.transforms.join.CoGroupByKey;
import com.google.cloud.dataflow.sdk.transforms.join.KeyedPCollectionTuple;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.TupleTag;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.engine.dataflow.datasources.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.util.List;
import java.util.UUID;


public class AddContextDataToRead {
    public static PCollection<KV<MutableGATKRead, ReadContextData>> Add(PCollection<MutableGATKRead> pReads,  RefAPISource refAPISource, RefAPIMetadata refAPIMetadata, VariantsFileSource variantsFileSource) {
        PCollection<Variant> pVariants = variantsFileSource.getAllVariants();
        PCollection<KV<MutableGATKRead, Iterable<Variant>>> kvReadVariants = KeyVariantsByRead.Key(pVariants, pReads);
        PCollection<KV<MutableGATKRead, ReferenceBases>> kvReadRefBases =
                APIToRefBasesKeyedByRead.addBases(pReads, refAPISource, refAPIMetadata);
        return Join(pReads, kvReadRefBases, kvReadVariants);

    }

    protected static PCollection<KV<MutableGATKRead, ReadContextData>> Join(PCollection<MutableGATKRead> pReads, PCollection<KV<MutableGATKRead, ReferenceBases>> kvReadRefBases, PCollection<KV<MutableGATKRead, Iterable<Variant>>> kvReadVariants) {
        // We could add a check that all of the reads in kvReadRefBases, pVariants, and pReads are the same.
        PCollection<KV<UUID, Iterable<Variant>>> UUIDVariants = kvReadVariants.apply(ParDo.of(new DoFn<KV<MutableGATKRead, Iterable<Variant>>, KV<UUID, Iterable<Variant>>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                MutableGATKRead r = c.element().getKey();
                Iterable<Variant> variants = c.element().getValue();
                List<Variant> lVariants = Lists.newArrayList(variants);
                //if (lVariants.isEmpty()) {
                //    lVariants = l
                //}
                c.output(KV.of(r.getUUID(), variants));
            }
        })).setName("KvUUIDiVariants");

        PCollection<KV<UUID, ReferenceBases>> UUIDRefBases = kvReadRefBases.apply(ParDo.of(new DoFn<KV<MutableGATKRead, ReferenceBases>, KV<UUID, ReferenceBases>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                MutableGATKRead r = c.element().getKey();
                c.output(KV.of(r.getUUID(), c.element().getValue()));
            }
        })).setName("KvUUIDRefBases");

        final TupleTag<Iterable<Variant>> variantTag = new TupleTag<>();
        final TupleTag<ReferenceBases> referenceTag = new TupleTag<>();

        PCollection<KV<UUID, CoGbkResult>> coGbkInput = KeyedPCollectionTuple
                .of(variantTag, UUIDVariants)
                .and(referenceTag, UUIDRefBases).apply(CoGroupByKey.<UUID>create());

        PCollection<KV<UUID, ReadContextData>> UUIDcontext = coGbkInput.apply(ParDo.of(new DoFn<KV<UUID, CoGbkResult>, KV<UUID, ReadContextData>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                Iterable<Iterable<Variant>> variants = c.element().getValue().getAll(variantTag);
                Iterable<ReferenceBases> referenceBases = c.element().getValue().getAll(referenceTag);

                // TODO: We should be enforcing that this is a singleton with getOnly somehow...
                List<Iterable<Variant>> vList = Lists.newArrayList(variants);
                if (vList.isEmpty()) {
                    vList.add(Lists.newArrayList());
                }
                if (vList.size() > 1) {
                    throw new GATKException("expected one list of varints, got " + vList.size());
                }
                List<ReferenceBases> bList = Lists.newArrayList(referenceBases);
                if (bList.size() != 1) {
                    throw new GATKException("expected one ReferenceBases, got " + bList.size());
                }

                c.output(KV.of(c.element().getKey(), new ReadContextData(bList.get(0), vList.get(0))));
            }
        })).setName("kVUUIDReadContextData");

        // Now add the reads back in.
        PCollection<KV<UUID, MutableGATKRead>> UUIDRead = pReads.apply(ParDo.of(new DoFn<MutableGATKRead, KV<UUID, MutableGATKRead>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                c.output(KV.of(c.element().getUUID(), c.element()));
            }
        })).setName("AddReadsBackIn");
        final TupleTag<MutableGATKRead> readTag = new TupleTag<>();
        final TupleTag<ReadContextData> contextDataTag = new TupleTag<>();

        PCollection<KV<UUID, CoGbkResult>> coGbkfull = KeyedPCollectionTuple
                .of(readTag, UUIDRead)
                .and(contextDataTag, UUIDcontext).apply(CoGroupByKey.<UUID>create());

        return coGbkfull.apply(ParDo.of(new DoFn<KV<UUID, CoGbkResult>, KV<MutableGATKRead, ReadContextData>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                Iterable<MutableGATKRead> reads = c.element().getValue().getAll(readTag);
                Iterable<ReadContextData> contextDatas = c.element().getValue().getAll(contextDataTag);

                // TODO: We should be enforcing that this is a singleton with getOnly somehow...
                List<MutableGATKRead> rList = Lists.newArrayList();
                for (MutableGATKRead r : reads) {
                    rList.add(r);
                }
                if (rList.size() != 1) {
                    throw new GATKException("expected one Read, got " + rList.size());
                }
                List<ReadContextData> cList = Lists.newArrayList();
                for (ReadContextData cd : contextDatas) {
                    cList.add(cd);
                }
                if (cList.size() != 1) {
                    throw new GATKException("expected one ReadContextData, got " + cList.size());
                }

                c.output(KV.of(rList.get(0), cList.get(0)));
            }
        })).setName("lastAddContextDataToRead");
    }

}