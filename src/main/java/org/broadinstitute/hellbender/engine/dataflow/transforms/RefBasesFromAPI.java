package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPIMetadata;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

public class RefBasesFromAPI {
    public static PCollection<KV<ReferenceBases, Iterable<MutableGATKRead>>> GetBases(PCollection<KV<ReferenceShard, Iterable<MutableGATKRead>>> reads,
                                                                                      RefAPISource refAPISource,
                                                                                      RefAPIMetadata refAPIMetadata) {
        PCollectionView<RefAPISource> sourceView = reads.getPipeline().apply(Create.of(refAPISource)).apply(View.<RefAPISource>asSingleton());
        PCollectionView<RefAPIMetadata> dataView = reads.getPipeline().apply(Create.of(refAPIMetadata)).apply(View.<RefAPIMetadata>asSingleton());
        return reads.apply(ParDo.withSideInputs(sourceView, dataView).of(
                new DoFn<KV<ReferenceShard, Iterable<MutableGATKRead>>, KV<ReferenceBases, Iterable<MutableGATKRead>>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                final ReferenceShard shard = c.element().getKey();
                final Iterable<MutableGATKRead> reads = c.element().getValue();
                int min = Integer.MAX_VALUE;
                int max = 1;
                for (MutableGATKRead r : reads) {
                    if (r.getStart() < min) {
                        min = r.getStart();
                    }
                    if (r.getEnd() > max) {
                        max = r.getEnd();
                    }
                }
                SimpleInterval interval = new SimpleInterval(shard.getContig(), min, max);

                ReferenceBases bases = c.sideInput(sourceView).getReferenceBases(c.getPipelineOptions(), c.sideInput(dataView), interval);
                c.output(KV.of(bases, reads));
            }
        })).setName("RefBasesFromAPI");
    }
}