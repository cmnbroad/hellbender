package org.broadinstitute.hellbender.utils.reference;

import com.google.cloud.dataflow.sdk.coders.AvroCoder;
import com.google.cloud.dataflow.sdk.coders.DefaultCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;
import java.util.Arrays;

public class ReferenceBases implements Serializable {

    private final byte[] bases;
    private final SimpleInterval interval;

    public ReferenceBases( final byte[] bases, final SimpleInterval interval ) {
        this.bases = bases;
        this.interval = interval;
    }

    @Override
    public String toString() {
        return "ReferenceBases{" +
                "bases=" + Arrays.toString(bases) +
                ", interval=" + interval +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ReferenceBases that = (ReferenceBases) o;

        if (!Arrays.equals(getBases(), that.getBases())) return false;
        return getInterval().equals(that.getInterval());

    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(getBases());
        result = 31 * result + getInterval().hashCode();
        return result;
    }

    public byte[] getBases() {
        return bases;
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    public ReferenceBases getSubset(SimpleInterval interval) {
        // I don't think this is quite right...

        // First check bounds, then if it's a proper subset, then return that subset.

        // Correct check?
        if (!this.interval.contains(interval)) {
            throw new GATKException("Reference doesn't match read...");
        }
        int start = interval.getStart() - this.interval.getStart();
        int end = interval.getEnd() - this.interval.getStart();
        return new ReferenceBases(Arrays.copyOfRange(this.bases, start, end), interval);
    }
}
