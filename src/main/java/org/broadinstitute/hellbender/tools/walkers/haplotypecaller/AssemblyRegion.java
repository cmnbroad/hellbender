package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.GenomeLoc;

import java.util.List;

/**
 * Interface that describes the region of the genome that gets assembled by the local assembly engine.
 */
public interface AssemblyRegion {

    /**
     * Get the span of this region including the extension value
     */
    GenomeLoc getExtendedSpan();

    /**
     * Get the raw span of this region (excluding the extension)
     */
    GenomeLoc getSpan();

    /**
     * Get an unmodifiable list of reads currently in this active region.
     *
     * The reads are sorted by their coordinate position
    */
     List<SAMRecord> getReads();
}
