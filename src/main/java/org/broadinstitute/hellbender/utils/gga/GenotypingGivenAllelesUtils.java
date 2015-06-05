package org.broadinstitute.hellbender.utils.gga;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Compendium of utils to work in GENOTYPE_GIVEN_ALLELES mode.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class GenotypingGivenAllelesUtils {

    /**
     * Create the list of artificial GGA-mode haplotypes by injecting each of the provided alternate alleles into the reference haplotype
     *
     * @param refHaplotype the reference haplotype
     * @param givenHaplotypes the list of alternate alleles in VariantContexts
     * @param activeRegionWindow the window containing the reference haplotype
     *
     * @return a non-null list of haplotypes
     */
    public static List<Haplotype> composeGivenHaplotypes(final Haplotype refHaplotype, final List<VariantContext> givenHaplotypes, final GenomeLoc activeRegionWindow) {
        if (refHaplotype == null) throw new IllegalArgumentException("the reference haplotype cannot be null");
        if (givenHaplotypes == null) throw new IllegalArgumentException("given haplotypes cannot be null");
        if (activeRegionWindow == null) throw new IllegalArgumentException("active region window cannot be null");
        if (activeRegionWindow.size() != refHaplotype.length()) throw new IllegalArgumentException("inconsistent reference haplotype and active region window");

        final Set<Haplotype> returnHaplotypes = new LinkedHashSet<>();
        final int activeRegionStart = refHaplotype.getAlignmentStartHapwrtRef();

        for( final VariantContext compVC : givenHaplotypes ) {
            if (!GATKVariantContextUtils.overlapsRegion(compVC, activeRegionWindow))
                throw new IllegalArgumentException(" some variant provided does not overlap with active region window");
            for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                final Haplotype insertedRefHaplotype = refHaplotype.insertAllele(compVC.getReference(), compAltAllele, activeRegionStart + compVC.getStart() - activeRegionWindow.getStart(), compVC.getStart());
                if( insertedRefHaplotype != null ) { // can be null if the requested allele can't be inserted into the haplotype
                    returnHaplotypes.add(insertedRefHaplotype);
                }
            }
        }

        return new ArrayList<>(returnHaplotypes);
    }
}
