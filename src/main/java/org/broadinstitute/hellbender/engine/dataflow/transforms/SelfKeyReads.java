package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;

import java.util.UUID;

public class SelfKeyReads extends PTransform<PCollection<MutableGATKRead>, PCollection<KV<UUID, MutableGATKRead>>> {
    private static final long serialVersionUID = 1L;

    @Override
    public PCollection<KV<UUID, MutableGATKRead>> apply(PCollection<MutableGATKRead> input) {
        return input.apply(ParDo.of(new DoFn<MutableGATKRead, KV<UUID, MutableGATKRead>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                c.output(KV.of(c.element().getUUID(), c.element()));
            }
        })).setName("SelfKeyReads");
    }
}
