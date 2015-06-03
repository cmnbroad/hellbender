package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPIMetadata;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.engine.dataflow.transforms.GroupReadsForRef;
import org.broadinstitute.hellbender.engine.dataflow.transforms.RefBasesFromAPI;
import org.broadinstitute.hellbender.engine.dataflow.transforms.GroupReadWithRefBases;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.read.MutableGATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;


public class APIToRefBasesKeyedByRead {
    public static PCollection<KV<MutableGATKRead, ReferenceBases>> addBases(PCollection<MutableGATKRead> pReads, RefAPISource refAPISource, RefAPIMetadata refAPIMetadata) {
        PCollection<KV<ReferenceShard, Iterable<MutableGATKRead>>> shardAndRead = pReads.apply(new GroupReadsForRef());
        PCollection<KV<ReferenceBases, Iterable<MutableGATKRead>>> groupedReads =
                RefBasesFromAPI.GetBases(shardAndRead, refAPISource, refAPIMetadata);
        return groupedReads.apply(new GroupReadWithRefBases());
    }
}
