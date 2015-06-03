package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

public class GroupReadWithRefBases extends PTransform<PCollection<KV<ReferenceBases, Iterable<MutableGATKRead>>>, PCollection<KV<MutableGATKRead, ReferenceBases>>> {
    private static final long serialVersionUID = 1L;

    @Override
    public PCollection<KV<MutableGATKRead, ReferenceBases>> apply(PCollection<KV<ReferenceBases, Iterable<MutableGATKRead>>> input) {
        return input.apply(ParDo.of(new DoFn<KV<ReferenceBases, Iterable<MutableGATKRead>>, KV<MutableGATKRead, ReferenceBases>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                // Each element of the PCollection is a set of reads keyed by a reference shard
                // The shard MUST have all of the reference bases for ALL of the reads. If not
                // it's an error.
                final ReferenceBases shard = c.element().getKey();
                final Iterable<MutableGATKRead> reads = c.element().getValue();
                // For every read, find the subset of the reference that matches it.
                for (MutableGATKRead r : reads) {
                    final ReferenceBases subset = shard.getSubset(new SimpleInterval(r));
                    c.output(KV.of(r, subset));
                }
            }
        })).setName("GroupReadWithRefBases");
    }
}
