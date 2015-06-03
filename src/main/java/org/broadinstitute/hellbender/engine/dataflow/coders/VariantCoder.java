package org.broadinstitute.hellbender.engine.dataflow.coders;


import com.google.cloud.dataflow.sdk.coders.CustomCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.variant.SkeletonVariant;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class VariantCoder extends CustomCoder<Variant> {

    @Override
    public void encode( Variant value, OutputStream outStream, Context context ) throws IOException {
        final Boolean isSkeletonVariant = value.getClass() == SkeletonVariant.class;
        SerializableCoder.of(Boolean.class).encode(isSkeletonVariant, outStream, context);

        if ( isSkeletonVariant ) {
            SerializableCoder.of(SkeletonVariant.class).encode(((SkeletonVariant) value), outStream, context);
        }
        else {
            SerializableCoder.of(VariantContextVariantAdapter.class).encode(((VariantContextVariantAdapter) value), outStream, context);
        }
    }

    @SuppressWarnings("unchecked")
    @Override
    public Variant decode( InputStream inStream, Context context ) throws IOException {
        final Boolean isSkeletonVariant = SerializableCoder.of(Boolean.class).decode(inStream, context);

        if ( isSkeletonVariant ) {
            return SerializableCoder.of(SkeletonVariant.class).decode(inStream, context);
        }
        else {
            return SerializableCoder.of(VariantContextVariantAdapter.class).decode(inStream, context);
        }
    }
}
