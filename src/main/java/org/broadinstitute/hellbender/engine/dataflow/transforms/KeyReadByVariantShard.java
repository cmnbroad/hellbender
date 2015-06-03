package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantShard;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;

import java.util.List;

public class KeyReadByVariantShard extends PTransform<PCollection<MutableGATKRead>, PCollection<KV<VariantShard, MutableGATKRead>>> {
    private static final long serialVersionUID = 1L;

    @Override
    public PCollection<KV<VariantShard, MutableGATKRead>> apply( PCollection<MutableGATKRead> input ) {
        return input.apply(ParDo.of(new DoFn<MutableGATKRead, KV<VariantShard, MutableGATKRead>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                List<VariantShard> shards = VariantShard.getVariantShardsFromInterval(c.element());
                for (VariantShard shard : shards) {
                    c.output(KV.of(shard, c.element()));
                }
            }
        }).named("KeyReadByVariantShard"));
    }
}
