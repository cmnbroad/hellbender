package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantShard;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.util.List;

public class KeyVariantByVariantShard extends PTransform<PCollection<Variant>, PCollection<KV<VariantShard, Variant>>> {
    private static final long serialVersionUID = 1L;

    @Override
    public PCollection<KV<VariantShard, Variant>> apply( PCollection<Variant> input ) {
        return input.apply(ParDo.of(new DoFn<Variant, KV<VariantShard, Variant>>() {
            @Override
            public void processElement( ProcessContext c ) throws Exception {
                List<VariantShard> shards = VariantShard.getVariantShardsFromInterval(c.element());
                for (VariantShard shard : shards) {
                    c.output(KV.of(shard, c.element()));
                }

            }
        }).named("KeyVariantByVariantShard"));
    }
}
