package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.GroupByKey;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;

public class GroupReadsForRef extends PTransform<PCollection<MutableGATKRead>, PCollection<KV<ReferenceShard, Iterable<MutableGATKRead>>>> {
    private static final long serialVersionUID = 1L;

    @Override
    public PCollection<KV<ReferenceShard, Iterable<MutableGATKRead>>> apply(PCollection<MutableGATKRead> input) {
        PCollection<KV<ReferenceShard, MutableGATKRead>> keyReadByReferenceShard = input.apply(ParDo.of(new DoFn<MutableGATKRead, KV<ReferenceShard, MutableGATKRead>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                ReferenceShard shard = ReferenceShard.getShardNumberFromInterval(c.element());
                c.output(KV.of(shard, c.element()));
            }
        }).named("KeyReadByReferenceShard"));
        return keyReadByReferenceShard.apply(GroupByKey.<ReferenceShard, MutableGATKRead>create());
    }
}

