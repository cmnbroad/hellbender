package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.transforms.PairReadsAndVariants;
import org.broadinstitute.hellbender.engine.dataflow.transforms.RemoveDuplicatePairedReadVariants;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;
import org.broadinstitute.hellbender.utils.variant.Variant;


public class KeyVariantsByRead {
    public static PCollection<KV<MutableGATKRead, Iterable<Variant>>> Key(PCollection<Variant> pVariants, PCollection<MutableGATKRead> pReads) {
        PCollection<KV<MutableGATKRead, Variant>> readVariants = PairReadsAndVariants.Pair(pReads, pVariants);

        // At this point, we ALMOST have what we want, but we need to remove duplicates KV<Read, Variant> pairs.
        // And we need to group by Read. Both of these require having deterministic coding, so we need to switch to
        // UUIDS.
        return readVariants.apply(new RemoveDuplicatePairedReadVariants());
    }
}
