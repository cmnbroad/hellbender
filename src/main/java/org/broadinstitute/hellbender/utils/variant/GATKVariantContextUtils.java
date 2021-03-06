package org.broadinstitute.hellbender.utils.variant;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

public final class GATKVariantContextUtils {

    private static Logger logger = LogManager.getLogger(GATKVariantContextUtils.class);

    public final static String MERGE_FILTER_PREFIX = "filterIn";
    public final static String MERGE_REF_IN_ALL = "ReferenceInAll";
    public final static String MERGE_FILTER_IN_ALL = "FilteredInAll";
    public final static String MERGE_INTERSECTION = "Intersection";

    /**
     * Checks whether a variant-context overlaps with a region.
     *
     * <p>
     *     No event overlaps an unmapped region.
     * </p>
     *
     * @param variantContext variant-context to test the overlap with.
     * @param region region to test the overlap with.
     *
     * @throws IllegalArgumentException if either region or event is {@code null}.
     *
     * @return {@code true} if there is an overlap between the event described and the active region provided.
     */
    public static boolean overlapsRegion(final VariantContext variantContext, final GenomeLoc region) {
        if (region == null) throw new IllegalArgumentException("the active region provided cannot be null");
        if (variantContext == null) throw new IllegalArgumentException("the variant context provided cannot be null");
        if (region.isUnmapped())
            return false;
        if (variantContext.getEnd() < region.getStart())
            return false;
        if (variantContext.getStart() > region.getStop())
            return false;
        if (!variantContext.getContig().equals(region.getContig()))
            return false;
        return true;
    }

    private static boolean hasPLIncompatibleAlleles(final Collection<Allele> alleleSet1, final Collection<Allele> alleleSet2) {
        final Iterator<Allele> it1 = alleleSet1.iterator();
        final Iterator<Allele> it2 = alleleSet2.iterator();

        while ( it1.hasNext() && it2.hasNext() ) {
            final Allele a1 = it1.next();
            final Allele a2 = it2.next();
            if ( ! a1.equals(a2) )
                return true;
        }

        // by this point, at least one of the iterators is empty.  All of the elements
        // we've compared are equal up until this point.  But it's possible that the
        // sets aren't the same size, which is indicated by the test below.  If they
        // are of the same size, though, the sets are compatible
        return it1.hasNext() || it2.hasNext();
    }

    /**
     * Determines the common reference allele
     *
     * @param VCs    the list of VariantContexts
     * @param loc    if not null, ignore records that do not begin at this start location
     * @return possibly null Allele
     */
    public static Allele determineReferenceAllele(final List<VariantContext> VCs, final GenomeLoc loc) {
        Allele ref = null;

        for ( final VariantContext vc : VCs ) {
            if ( contextMatchesLoc(vc, loc) ) {
                final Allele myRef = vc.getReference();
                if ( ref == null || ref.length() < myRef.length() )
                    ref = myRef;
                else if ( ref.length() == myRef.length() && ! ref.equals(myRef) )
                    throw new TribbleException(String.format("The provided variant file(s) have inconsistent references for the same position(s) at %s:%d, %s vs. %s", vc.getContig(), vc.getStart(), ref, myRef));
            }
        }

        return ref;
    }

    /**
     * Calculates the total ploidy of a variant context as the sum of all plodies across genotypes.
     * @param vc the target variant context.
     * @param defaultPloidy the default ploidy to be assume when there is no ploidy information for a genotype.
     * @return never {@code null}.
     */
    public static int totalPloidy(final VariantContext vc, final int defaultPloidy) {
        if (vc == null)
            throw new IllegalArgumentException("the vc provided cannot be null");
        if (defaultPloidy < 0)
            throw new IllegalArgumentException("the default ploidy must 0 or greater");
        int result = 0;
        for (final Genotype genotype : vc.getGenotypes()) {
            final int declaredPloidy = genotype.getPloidy();
            result += declaredPloidy <= 0 ? defaultPloidy : declaredPloidy;
        }

        return result;
    }

    public enum GenotypeMergeType {
        /**
         * Make all sample genotypes unique by file. Each sample shared across RODs gets named sample.ROD.
         */
        UNIQUIFY,
        /**
         * Take genotypes in priority order (see the priority argument).
         */
        PRIORITIZE,
        /**
         * Take the genotypes in any order.
         */
        UNSORTED,
        /**
         * Require that all samples/genotypes be unique between all inputs.
         */
        REQUIRE_UNIQUE
    }

    public enum FilteredRecordMergeType {
        /**
         * Union - leaves the record if any record is unfiltered.
         */
        KEEP_IF_ANY_UNFILTERED,
        /**
         * Requires all records present at site to be unfiltered. VCF files that don't contain the record don't influence this.
         */
        KEEP_IF_ALL_UNFILTERED,
        /**
         * If any record is present at this site (regardless of possibly being filtered), then all such records are kept and the filters are reset.
         */
        KEEP_UNCONDITIONAL
    }

    /**
     * Returns true iff VC is an non-complex indel where every allele represents an expansion or
     * contraction of a series of identical bases in the reference.
     *
     * For example, suppose the ref bases are CTCTCTGA, which includes a 3x repeat of CTCTCT
     *
     * If VC = -/CT, then this function returns true because the CT insertion matches exactly the
     * upcoming reference.
     * If VC = -/CTA then this function returns false because the CTA isn't a perfect match
     *
     * Now consider deletions:
     *
     * If VC = CT/- then again the same logic applies and this returns true
     * The case of CTA/- makes no sense because it doesn't actually match the reference bases.
     *
     * The logic of this function is pretty simple.  Take all of the non-null alleles in VC.  For
     * each insertion allele of n bases, check if that allele matches the next n reference bases.
     * For each deletion allele of n bases, check if this matches the reference bases at n - 2 n,
     * as it must necessarily match the first n bases.  If this test returns true for all
     * alleles you are a tandem repeat, otherwise you are not.
     *
     * @param vc
     * @param refBasesStartingAtVCWithPad not this is assumed to include the PADDED reference
     * @return
     */
    public static boolean isTandemRepeat(final VariantContext vc, final byte[] refBasesStartingAtVCWithPad) {
        final String refBasesStartingAtVCWithoutPad = new String(refBasesStartingAtVCWithPad).substring(1);
        if ( ! vc.isIndel() ) // only indels are tandem repeats
            return false;

        final Allele ref = vc.getReference();

        for ( final Allele allele : vc.getAlternateAlleles() ) {
            if ( ! isRepeatAllele(ref, allele, refBasesStartingAtVCWithoutPad) )
                return false;
        }

        // we've passed all of the tests, so we are a repeat
        return true;
    }

    /**
     *
     * @param vc
     * @param refBasesStartingAtVCWithPad
     * @return
     */
    public static Pair<List<Integer>,byte[]> getNumTandemRepeatUnits(final VariantContext vc, final byte[] refBasesStartingAtVCWithPad) {
        final boolean VERBOSE = false;
        final String refBasesStartingAtVCWithoutPad = new String(refBasesStartingAtVCWithPad).substring(1);
        if ( ! vc.isIndel() ) // only indels are tandem repeats
            return null;

        final Allele refAllele = vc.getReference();
        final byte[] refAlleleBases = Arrays.copyOfRange(refAllele.getBases(), 1, refAllele.length());

        byte[] repeatUnit = null;
        final ArrayList<Integer> lengths = new ArrayList<>();

        for ( final Allele allele : vc.getAlternateAlleles() ) {
            Pair<int[],byte[]> result = getNumTandemRepeatUnits(refAlleleBases, Arrays.copyOfRange(allele.getBases(), 1, allele.length()), refBasesStartingAtVCWithoutPad.getBytes());

            final int[] repetitionCount = result.getLeft();
            // repetition count = 0 means allele is not a tandem expansion of context
            if (repetitionCount[0] == 0 || repetitionCount[1] == 0)
                return null;

            if (lengths.size() == 0) {
                lengths.add(repetitionCount[0]); // add ref allele length only once
            }
            lengths.add(repetitionCount[1]);  // add this alt allele's length

            repeatUnit = result.getRight();
            if (VERBOSE) {
                System.out.println("RefContext:"+refBasesStartingAtVCWithoutPad);
                System.out.println("Ref:"+refAllele.toString()+" Count:" + String.valueOf(repetitionCount[0]));
                System.out.println("Allele:"+allele.toString()+" Count:" + String.valueOf(repetitionCount[1]));
                System.out.println("RU:"+new String(repeatUnit));
            }
        }

        return new MutablePair<>(lengths,repeatUnit);
    }

    public static Pair<int[],byte[]> getNumTandemRepeatUnits(final byte[] refBases, final byte[] altBases, final byte[] remainingRefContext) {
         /* we can't exactly apply same logic as in basesAreRepeated() to compute tandem unit and number of repeated units.
           Consider case where ref =ATATAT and we have an insertion of ATAT. Natural description is (AT)3 -> (AT)2.
         */

        byte[] longB;
        // find first repeat unit based on either ref or alt, whichever is longer
        if (altBases.length > refBases.length)
            longB = altBases;
        else
            longB = refBases;

        // see if non-null allele (either ref or alt, whichever is longer) can be decomposed into several identical tandem units
        // for example, -*,CACA needs to first be decomposed into (CA)2
        final int repeatUnitLength = findRepeatedSubstring(longB);
        final byte[] repeatUnit = Arrays.copyOf(longB, repeatUnitLength);

        final int[] repetitionCount = new int[2];
        // look for repetitions forward on the ref bases (i.e. starting at beginning of ref bases)
        int repetitionsInRef = findNumberOfRepetitions(repeatUnit, refBases, true);
        repetitionCount[0] = findNumberOfRepetitions(repeatUnit, ArrayUtils.addAll(refBases, remainingRefContext), true)-repetitionsInRef;
        repetitionCount[1] = findNumberOfRepetitions(repeatUnit, ArrayUtils.addAll(altBases, remainingRefContext), true)-repetitionsInRef;

        return new MutablePair<>(repetitionCount, repeatUnit);

    }

    /**
     * Find out if a string can be represented as a tandem number of substrings.
     * For example ACTACT is a 2-tandem of ACT,
     * but ACTACA is not.
     *
     * @param bases                 String to be tested
     * @return                      Length of repeat unit, if string can be represented as tandem of substring (if it can't
     *                              be represented as one, it will be just the length of the input string)
     */
    public static int findRepeatedSubstring(byte[] bases) {

        int repLength;
        for (repLength=1; repLength <=bases.length; repLength++) {
            final byte[] candidateRepeatUnit = Arrays.copyOf(bases,repLength);
            boolean allBasesMatch = true;
            for (int start = repLength; start < bases.length; start += repLength ) {
                // check that remaining of string is exactly equal to repeat unit
                final byte[] basePiece = Arrays.copyOfRange(bases,start,start+candidateRepeatUnit.length);
                if (!Arrays.equals(candidateRepeatUnit, basePiece)) {
                    allBasesMatch = false;
                    break;
                }
            }
            if (allBasesMatch)
                return repLength;
        }

        return repLength;
    }

    /**
     * Helper routine that finds number of repetitions a string consists of.
     * For example, for string ATAT and repeat unit AT, number of repetitions = 2
     * @param repeatUnit             Substring
     * @param testString             String to test
     * @oaram lookForward            Look for repetitions forward (at beginning of string) or backward (at end of string)
     * @return                       Number of repetitions (0 if testString is not a concatenation of n repeatUnit's
     */
    public static int findNumberOfRepetitions(byte[] repeatUnit, byte[] testString, boolean lookForward) {
        int numRepeats = 0;
        if (lookForward) {
            // look forward on the test string
            for (int start = 0; start < testString.length; start += repeatUnit.length) {
                int end = start + repeatUnit.length;
                byte[] unit = Arrays.copyOfRange(testString,start, end);
                if(Arrays.equals(unit,repeatUnit))
                    numRepeats++;
                else
                    break;
            }
            return numRepeats;
        }

        // look backward. For example, if repeatUnit = AT and testString = GATAT, number of repeat units is still 2
        // look forward on the test string
        for (int start = testString.length - repeatUnit.length; start >= 0; start -= repeatUnit.length) {
            int end = start + repeatUnit.length;
            byte[] unit = Arrays.copyOfRange(testString,start, end);
            if(Arrays.equals(unit,repeatUnit))
                numRepeats++;
            else
                break;
        }
        return numRepeats;
    }

    /**
     * Helper function for isTandemRepeat that checks that allele matches somewhere on the reference
     * @param ref
     * @param alt
     * @param refBasesStartingAtVCWithoutPad
     * @return
     */
    protected static boolean isRepeatAllele(final Allele ref, final Allele alt, final String refBasesStartingAtVCWithoutPad) {
        if ( ! Allele.oneIsPrefixOfOther(ref, alt) )
            return false; // we require one allele be a prefix of another

        if ( ref.length() > alt.length() ) { // we are a deletion
            return basesAreRepeated(ref.getBaseString(), alt.getBaseString(), refBasesStartingAtVCWithoutPad, 2);
        } else { // we are an insertion
            return basesAreRepeated(alt.getBaseString(), ref.getBaseString(), refBasesStartingAtVCWithoutPad, 1);
        }
    }

    protected static boolean basesAreRepeated(final String l, final String s, final String ref, final int minNumberOfMatches) {
        final String potentialRepeat = l.substring(s.length()); // skip s bases

        for ( int i = 0; i < minNumberOfMatches; i++) {
            final int start = i * potentialRepeat.length();
            final int end = (i+1) * potentialRepeat.length();
            if ( ref.length() < end )
                return false; // we ran out of bases to test
            final String refSub = ref.substring(start, end);
            if ( ! refSub.equals(potentialRepeat) )
                return false; // repeat didn't match, fail
        }

        return true; // we passed all tests, we matched
    }

    public enum GenotypeAssignmentMethod {
        /**
         * set all of the genotype GT values to NO_CALL
         */
        SET_TO_NO_CALL,

        /**
         * Use the subsetted PLs to greedily assigned genotypes
         */
        USE_PLS_TO_ASSIGN,

        /**
         * Try to match the original GT calls, if at all possible
         *
         * Suppose I have 3 alleles: A/B/C and the following samples:
         *
         *       original_GT best_match to A/B best_match to A/C
         * S1 => A/A A/A A/A
         * S2 => A/B A/B A/A
         * S3 => B/B B/B A/A
         * S4 => B/C A/B A/C
         * S5 => C/C A/A C/C
         *
         * Basically, all alleles not in the subset map to ref.  It means that het-alt genotypes
         * when split into 2 bi-allelic variants will be het in each, which is good in some cases,
         * rather than the undetermined behavior when using the PLs to assign, which could result
         * in hom-var or hom-ref for each, depending on the exact PL values.
         */
        BEST_MATCH_TO_ORIGINAL,

        /**
         * do not even bother changing the GTs
         */
        DO_NOT_ASSIGN_GENOTYPES
    }

    /**
     * Subset the samples in VC to reference only information with ref call alleles
     *
     * Preserves DP if present
     *
     * @param vc the variant context to subset down to
     * @param ploidy ploidy to use if a genotype doesn't have any alleles
     * @return a GenotypesContext
     */
    public static GenotypesContext subsetToRefOnly(final VariantContext vc, final int ploidy) {
        if ( vc == null ) throw new IllegalArgumentException("vc cannot be null");
        if ( ploidy < 1 ) throw new IllegalArgumentException("ploidy must be >= 1 but got " + ploidy);

        // the genotypes with PLs
        final GenotypesContext oldGTs = vc.getGenotypes();

        // optimization: if no input genotypes, just exit
        if (oldGTs.isEmpty()) return oldGTs;

        // the new genotypes to create
        final GenotypesContext newGTs = GenotypesContext.create(oldGTs.size());

        final Allele ref = vc.getReference();
        final List<Allele> diploidRefAlleles = Arrays.asList(ref, ref);

        // create the new genotypes
        for ( final Genotype g : vc.getGenotypes() ) {
            final int gPloidy = g.getPloidy() == 0 ? ploidy : g.getPloidy();
            final List<Allele> refAlleles = gPloidy == 2 ? diploidRefAlleles : Collections.nCopies(gPloidy, ref);
            final GenotypeBuilder gb = new GenotypeBuilder(g.getSampleName(), refAlleles);
            if ( g.hasDP() ) gb.DP(g.getDP());
            if ( g.hasGQ() ) gb.GQ(g.getGQ());
            newGTs.add(gb.make());
        }

        return newGTs;
    }

    public static Genotype removePLsAndAD(final Genotype g) {
        return ( g.hasLikelihoods() || g.hasAD() ) ? new GenotypeBuilder(g).noPL().noAD().make() : g;
    }

    //TODO consider refactor variant-context merging code so that we share as much as possible between
    //TODO simpleMerge and referenceConfidenceMerge
    //TODO likely using a separate helper class or hierarchy.
    /**
     * Merges VariantContexts into a single hybrid.  Takes genotypes for common samples in priority order, if provided.
     * If uniquifySamples is true, the priority order is ignored and names are created by concatenating the VC name with
     * the sample name
     *
     * @param unsortedVCs               collection of unsorted VCs
     * @param priorityListOfVCs         priority list detailing the order in which we should grab the VCs
     * @param filteredRecordMergeType   merge type for filtered records
     * @param genotypeMergeOptions      merge option for genotypes
     * @param annotateOrigin            should we annotate the set it came from?
     * @param printMessages             should we print messages?
     * @param setKey                    the key name of the set
     * @param filteredAreUncalled       are filtered records uncalled?
     * @param mergeInfoWithMaxAC        should we merge in info from the VC with maximum allele count?
     * @return new VariantContext       representing the merge of unsortedVCs
     */
    public static VariantContext simpleMerge(final Collection<VariantContext> unsortedVCs,
                                             final List<String> priorityListOfVCs,
                                             final FilteredRecordMergeType filteredRecordMergeType,
                                             final GenotypeMergeType genotypeMergeOptions,
                                             final boolean annotateOrigin,
                                             final boolean printMessages,
                                             final String setKey,
                                             final boolean filteredAreUncalled,
                                             final boolean mergeInfoWithMaxAC ) {
        int originalNumOfVCs = priorityListOfVCs == null ? 0 : priorityListOfVCs.size();
        return simpleMerge(unsortedVCs, priorityListOfVCs, originalNumOfVCs, filteredRecordMergeType, genotypeMergeOptions, annotateOrigin, printMessages, setKey, filteredAreUncalled, mergeInfoWithMaxAC);
    }

    /**
     * Merges VariantContexts into a single hybrid.  Takes genotypes for common samples in priority order, if provided.
     * If uniquifySamples is true, the priority order is ignored and names are created by concatenating the VC name with
     * the sample name.
     * simpleMerge does not verify any more unique sample names EVEN if genotypeMergeOptions == GenotypeMergeType.REQUIRE_UNIQUE. One should use
     * SampleUtils.verifyUniqueSamplesNames to check that before using simpleMerge.
     *
     * For more information on this method see: http://www.thedistractionnetwork.com/programmer-problem/
     *
     * @param unsortedVCs               collection of unsorted VCs
     * @param priorityListOfVCs         priority list detailing the order in which we should grab the VCs
     * @param filteredRecordMergeType   merge type for filtered records
     * @param genotypeMergeOptions      merge option for genotypes
     * @param annotateOrigin            should we annotate the set it came from?
     * @param printMessages             should we print messages?
     * @param setKey                    the key name of the set
     * @param filteredAreUncalled       are filtered records uncalled?
     * @param mergeInfoWithMaxAC        should we merge in info from the VC with maximum allele count?
     * @return new VariantContext       representing the merge of unsortedVCs
     */
    public static VariantContext simpleMerge(final Collection<VariantContext> unsortedVCs,
                                             final List<String> priorityListOfVCs,
                                             final int originalNumOfVCs,
                                             final FilteredRecordMergeType filteredRecordMergeType,
                                             final GenotypeMergeType genotypeMergeOptions,
                                             final boolean annotateOrigin,
                                             final boolean printMessages,
                                             final String setKey,
                                             final boolean filteredAreUncalled,
                                             final boolean mergeInfoWithMaxAC ) {
        if ( unsortedVCs == null || unsortedVCs.size() == 0 )
            return null;

        if (priorityListOfVCs != null && originalNumOfVCs != priorityListOfVCs.size())
            throw new IllegalArgumentException("the number of the original VariantContexts must be the same as the number of VariantContexts in the priority list");

        if ( annotateOrigin && priorityListOfVCs == null && originalNumOfVCs == 0)
            throw new IllegalArgumentException("Cannot merge calls and annotate their origins without a complete priority list of VariantContexts or the number of original VariantContexts");

        final List<VariantContext> preFilteredVCs = sortVariantContextsByPriority(unsortedVCs, priorityListOfVCs, genotypeMergeOptions);
        // Make sure all variant contexts are padded with reference base in case of indels if necessary
        List<VariantContext> VCs = preFilteredVCs.stream()
                .filter(vc -> !filteredAreUncalled || vc.isNotFiltered())
                .collect(Collectors.toList());

        if ( VCs.size() == 0 ) // everything is filtered out and we're filteredAreUncalled
            return null;

        // establish the baseline info from the first VC
        final VariantContext first = VCs.get(0);
        final String name = first.getSource();
        final Allele refAllele = determineReferenceAllele(VCs);

        final LinkedHashSet<Allele> alleles = new LinkedHashSet<>();
        final Set<String> filters = new HashSet<>();
        final Map<String, Object> attributes = new LinkedHashMap<>();
        final Set<String> inconsistentAttributes = new HashSet<>();
        final Set<String> variantSources = new HashSet<>(); // contains the set of sources we found in our set of VCs that are variant
        final Set<String> rsIDs = new LinkedHashSet<>(1); // most of the time there's one id

        VariantContext longestVC = first;
        int depth = 0;
        int maxAC = -1;
        final Map<String, Object> attributesWithMaxAC = new LinkedHashMap<>();
        double log10PError = CommonInfo.NO_LOG10_PERROR;
        boolean anyVCHadFiltersApplied = false;
        VariantContext vcWithMaxAC = null;
        GenotypesContext genotypes = GenotypesContext.create();

        // counting the number of filtered and variant VCs
        int nFiltered = 0;

        boolean remapped = false;

        // cycle through and add info from the other VCs, making sure the loc/reference matches
        for ( final VariantContext vc : VCs ) {
            if ( longestVC.getStart() != vc.getStart() )
                throw new IllegalStateException("BUG: attempting to merge VariantContexts with different start sites: first="+ first.toString() + " second=" + vc.toString());

            if ( VariantContextUtils.getSize(vc) > VariantContextUtils.getSize(longestVC) )
                longestVC = vc; // get the longest location

            nFiltered += vc.isFiltered() ? 1 : 0;
            if ( vc.isVariant() ) variantSources.add(vc.getSource());

            AlleleMapper alleleMapping = resolveIncompatibleAlleles(refAllele, vc, alleles);
            remapped = remapped || alleleMapping.needsRemapping();

            alleles.addAll(alleleMapping.values());

            mergeGenotypes(genotypes, vc, alleleMapping, genotypeMergeOptions == GenotypeMergeType.UNIQUIFY);

            // We always take the QUAL of the first VC with a non-MISSING qual for the combined value
            if ( log10PError == CommonInfo.NO_LOG10_PERROR )
                log10PError =  vc.getLog10PError();

            filters.addAll(vc.getFilters());
            anyVCHadFiltersApplied |= vc.filtersWereApplied();

            //
            // add attributes
            //
            // special case DP (add it up) and ID (just preserve it)
            //
            if (vc.hasAttribute(VCFConstants.DEPTH_KEY))
                depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            if ( vc.hasID() ) rsIDs.add(vc.getID());
            if (mergeInfoWithMaxAC && vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY)) {
                String rawAlleleCounts = vc.getAttributeAsString(VCFConstants.ALLELE_COUNT_KEY, null);
                // lets see if the string contains a "," separator
                if (rawAlleleCounts.contains(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)) {
                    final List<String> alleleCountArray = Arrays.asList(rawAlleleCounts.substring(1, rawAlleleCounts.length() - 1).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR));
                    for (final String alleleCount : alleleCountArray) {
                        final int ac = Integer.valueOf(alleleCount.trim());
                        if (ac > maxAC) {
                            maxAC = ac;
                            vcWithMaxAC = vc;
                        }
                    }
                } else {
                    final int ac = Integer.valueOf(rawAlleleCounts);
                    if (ac > maxAC) {
                        maxAC = ac;
                        vcWithMaxAC = vc;
                    }
                }
            }

            for (final Map.Entry<String, Object> p : vc.getAttributes().entrySet()) {
                final String key = p.getKey();
                final Object value = p.getValue();
                // only output annotations that have the same value in every input VC
                // if we don't like the key already, don't go anywhere
                if ( ! inconsistentAttributes.contains(key) ) {
                    final boolean alreadyFound = attributes.containsKey(key);
                    final Object boundValue = attributes.get(key);
                    final boolean boundIsMissingValue = alreadyFound && boundValue.equals(VCFConstants.MISSING_VALUE_v4);

                    if ( alreadyFound && ! boundValue.equals(value) && ! boundIsMissingValue ) {
                        // we found the value but we're inconsistent, put it in the exclude list
                        inconsistentAttributes.add(key);
                        attributes.remove(key);
                    } else if ( ! alreadyFound || boundIsMissingValue )  { // no value
                        attributes.put(key, value);
                    }
                }
            }
        }

        // if we have more alternate alleles in the merged VC than in one or more of the
        // original VCs, we need to strip out the GL/PLs (because they are no longer accurate), as well as allele-dependent attributes like AC,AF, and AD
        for ( final VariantContext vc : VCs ) {
            if (vc.getAlleles().size() == 1)
                continue;
            if ( hasPLIncompatibleAlleles(alleles, vc.getAlleles())) {
                if ( ! genotypes.isEmpty() ) {
                    logger.debug(String.format("Stripping PLs at %s:%d-%d due to incompatible alleles merged=%s vs. single=%s",
                            vc.getContig(), vc.getStart(), vc.getEnd(), alleles, vc.getAlleles()));
                }
                genotypes = stripPLsAndAD(genotypes);
                // this will remove stale AC,AF attributed from vc
                VariantContextUtils.calculateChromosomeCounts(vc, attributes, true);
                break;
            }
        }

        // take the VC with the maxAC and pull the attributes into a modifiable map
        if ( mergeInfoWithMaxAC && vcWithMaxAC != null ) {
            attributesWithMaxAC.putAll(vcWithMaxAC.getAttributes());
        }

        // if at least one record was unfiltered and we want a union, clear all of the filters
        if ( (filteredRecordMergeType == FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED && nFiltered != VCs.size()) || filteredRecordMergeType == FilteredRecordMergeType.KEEP_UNCONDITIONAL )
            filters.clear();


        if ( annotateOrigin ) { // we care about where the call came from
            String setValue;
            if ( nFiltered == 0 && variantSources.size() == originalNumOfVCs ) // nothing was unfiltered
                setValue = MERGE_INTERSECTION;
            else if ( nFiltered == VCs.size() )     // everything was filtered out
                setValue = MERGE_FILTER_IN_ALL;
            else if ( variantSources.isEmpty() )    // everyone was reference
                setValue = MERGE_REF_IN_ALL;
            else {
                final LinkedHashSet<String> s = new LinkedHashSet<>();
                for ( final VariantContext vc : VCs )
                    if ( vc.isVariant() )
                        s.add( vc.isFiltered() ? MERGE_FILTER_PREFIX + vc.getSource() : vc.getSource() );
                setValue = Utils.join("-", s);
            }

            if ( setKey != null ) {
                attributes.put(setKey, setValue);
                if( mergeInfoWithMaxAC && vcWithMaxAC != null ) {
                    attributesWithMaxAC.put(setKey, setValue);
                }
            }
        }

        if ( depth > 0 )
            attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));

        final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(",", rsIDs);

        final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(ID);
        builder.loc(longestVC.getContig(), longestVC.getStart(), longestVC.getEnd());
        builder.alleles(alleles);
        builder.genotypes(genotypes);
        builder.log10PError(log10PError);
        if ( anyVCHadFiltersApplied ) {
            builder.filters(filters.isEmpty() ? filters : new TreeSet<>(filters));
        }
        builder.attributes(new TreeMap<>(mergeInfoWithMaxAC ? attributesWithMaxAC : attributes));

        // Trim the padded bases of all alleles if necessary
        final VariantContext merged = builder.make();
        if ( printMessages && remapped ) System.out.printf("Remapped => %s%n", merged);
        return merged;
    }

    //TODO as part of a larger refactoring effort remapAlleles can be merged with createAlleleMapping.

    public static GenotypesContext stripPLsAndAD(final GenotypesContext genotypes) {
        final GenotypesContext newGs = GenotypesContext.create(genotypes.size());

        for ( final Genotype g : genotypes ) {
            newGs.add(removePLsAndAD(g));
        }

        return newGs;
    }

    private static Allele determineReferenceAllele(final List<VariantContext> VCs) {
        return determineReferenceAllele(VCs, null);
    }

    public static boolean contextMatchesLoc(final VariantContext vc, final GenomeLoc loc) {
        return loc == null || loc.getStart() == vc.getStart();
    }

    static private AlleleMapper resolveIncompatibleAlleles(final Allele refAllele, final VariantContext vc, final LinkedHashSet<Allele> allAlleles) {
        if ( refAllele.equals(vc.getReference()) )
            return new AlleleMapper(vc);
        else {
            final Map<Allele, Allele> map = createAlleleMapping(refAllele, vc, allAlleles);
            map.put(vc.getReference(), refAllele);
            return new AlleleMapper(map);
        }
    }

    //TODO as part of a larger refactoring effort {@link #createAlleleMapping} can be merged with {@link ReferenceConfidenceVariantContextMerger#remapAlleles}.
    /**
     * Create an allele mapping for the given context where its reference allele must (potentially) be extended to the given allele
     *
     * The refAllele is the longest reference allele seen at this start site.
     * So imagine it is:
     * refAllele: ACGTGA
     * myRef:     ACGT
     * myAlt:     A
     *
     * We need to remap all of the alleles in vc to include the extra GA so that
     * myRef => refAllele and myAlt => AGA
     *
     * @param refAllele          the new (extended) reference allele
     * @param oneVC              the Variant Context to extend
     * @param currentAlleles     the list of alleles already created
     * @return a non-null mapping of original alleles to new (extended) ones
     */
    private static Map<Allele, Allele> createAlleleMapping(final Allele refAllele,
                                                           final VariantContext oneVC,
                                                           final Collection<Allele> currentAlleles) {
        final Allele myRef = oneVC.getReference();
        if ( refAllele.length() <= myRef.length() ) throw new IllegalStateException("BUG: myRef="+myRef+" is longer than refAllele="+refAllele);

        final byte[] extraBases = Arrays.copyOfRange(refAllele.getBases(), myRef.length(), refAllele.length());

        final Map<Allele, Allele> map = new HashMap<>();
        for ( final Allele a : oneVC.getAlternateAlleles() ) {
            if ( isUsableAlternateAllele(a) ) {
                Allele extended = Allele.extend(a, extraBases);
                for ( final Allele b : currentAlleles )
                    if ( extended.equals(b) )
                        extended = b;
                map.put(a, extended);
            }
        }

        return map;
    }

    static private boolean isUsableAlternateAllele(final Allele allele) {
        return ! (allele.isReference() || allele.isSymbolic() );
    }

    public static List<VariantContext> sortVariantContextsByPriority(Collection<VariantContext> unsortedVCs, List<String> priorityListOfVCs, GenotypeMergeType mergeOption ) {
        if ( mergeOption == GenotypeMergeType.PRIORITIZE && priorityListOfVCs == null )
            throw new IllegalArgumentException("Cannot merge calls by priority with a null priority list");

        if ( priorityListOfVCs == null || mergeOption == GenotypeMergeType.UNSORTED )
            return new ArrayList<>(unsortedVCs);
        else {
            ArrayList<VariantContext> sorted = new ArrayList<>(unsortedVCs);
            Collections.sort(sorted, new CompareByPriority(priorityListOfVCs));
            return sorted;
        }
    }

    private static void mergeGenotypes(GenotypesContext mergedGenotypes, VariantContext oneVC, AlleleMapper alleleMapping, boolean uniquifySamples) {
        //TODO: should we add a check for cases when the genotypeMergeOption is REQUIRE_UNIQUE
        for ( final Genotype g : oneVC.getGenotypes() ) {
            final String name = mergedSampleName(oneVC.getSource(), g.getSampleName(), uniquifySamples);
            if ( ! mergedGenotypes.containsSample(name) ) {
                // only add if the name is new
                Genotype newG = g;

                if ( uniquifySamples || alleleMapping.needsRemapping() ) {
                    final List<Allele> alleles = alleleMapping.needsRemapping() ? alleleMapping.remap(g.getAlleles()) : g.getAlleles();
                    newG = new GenotypeBuilder(g).name(name).alleles(alleles).make();
                }

                mergedGenotypes.add(newG);
            }
        }
    }

    /**
     * Cached NO_CALL immutable lists where the position ith contains the list with i elements.
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    private static List<Allele>[] NOCALL_LISTS = new List[] {
            Collections.emptyList(),
            Collections.singletonList(Allele.NO_CALL),
            Collections.nCopies(2,Allele.NO_CALL)
    };

    /**
     * Code to ensure that {@link #NOCALL_LISTS} has enough entries beyond the requested ploidy
     * @param capacity the requested ploidy.
     */
    private static void ensureNoCallListsCapacity(final int capacity) {
        final int currentCapacity = NOCALL_LISTS.length - 1;
        if (currentCapacity >= capacity)
            return;
        NOCALL_LISTS = Arrays.copyOf(NOCALL_LISTS,Math.max(capacity,currentCapacity << 1) + 1);
        for (int i = currentCapacity + 1; i < NOCALL_LISTS.length; i++)
            NOCALL_LISTS[i] = Collections.nCopies(i,Allele.NO_CALL);
    }

    /**
     * Returns a {@link Allele#NO_CALL NO_CALL} allele list provided the ploidy.
     *
     * @param ploidy the required ploidy.
     *
     * @return never {@code null}, but an empty list if {@code ploidy} is equal or less than 0. The returned list
     *   might or might not be mutable.
     */
    public static List<Allele> noCallAlleles(final int ploidy) {
        if (NOCALL_LISTS.length <= ploidy)
            ensureNoCallListsCapacity(ploidy);
        return NOCALL_LISTS[ploidy];
    }


    public static String mergedSampleName(String trackName, String sampleName, boolean uniquify ) {
        return uniquify ? sampleName + "." + trackName : sampleName;
    }

    /**
     * Trim the alleles in inputVC forward and reverse, as requested
     *
     * @param inputVC a non-null input VC whose alleles might need a haircut
     * @param trimForward should we trim up the alleles from the forward direction?
     * @param trimReverse should we trim up the alleles from the reverse direction?
     * @return a non-null VariantContext (may be == to inputVC) with trimmed up alleles
     */
    public static VariantContext trimAlleles(final VariantContext inputVC, final boolean trimForward, final boolean trimReverse) {
        if ( inputVC == null ) throw new IllegalArgumentException("inputVC cannot be null");

        if ( inputVC.getNAlleles() <= 1 || inputVC.isSNP() )
            return inputVC;

        // see whether we need to trim common reference base from all alleles
        final int revTrim = trimReverse ? computeReverseClipping(inputVC.getAlleles(), inputVC.getReference().getDisplayString().getBytes()) : 0;
        final VariantContext revTrimVC = trimAlleles(inputVC, -1, revTrim);
        final int fwdTrim = trimForward ? computeForwardClipping(revTrimVC.getAlleles()) : -1;
        return trimAlleles(revTrimVC, fwdTrim, 0);
    }

    /**
     * Trim up alleles in inputVC, cutting out all bases up to fwdTrimEnd inclusive and
     * the last revTrim bases from the end
     *
     * @param inputVC a non-null input VC
     * @param fwdTrimEnd bases up to this index (can be -1) will be removed from the start of all alleles
     * @param revTrim the last revTrim bases of each allele will be clipped off as well
     * @return a non-null VariantContext (may be == to inputVC) with trimmed up alleles
     */
    protected static VariantContext trimAlleles(final VariantContext inputVC,
                                                final int fwdTrimEnd,
                                                final int revTrim) {
        if( fwdTrimEnd == -1 && revTrim == 0 ) // nothing to do, so just return inputVC unmodified
            return inputVC;

        final List<Allele> alleles = new LinkedList<>();
        final Map<Allele, Allele> originalToTrimmedAlleleMap = new HashMap<>();

        for (final Allele a : inputVC.getAlleles()) {
            if (a.isSymbolic()) {
                alleles.add(a);
                originalToTrimmedAlleleMap.put(a, a);
            } else {
                // get bases for current allele and create a new one with trimmed bases
                final byte[] newBases = Arrays.copyOfRange(a.getBases(), fwdTrimEnd+1, a.length()-revTrim);
                final Allele trimmedAllele = Allele.create(newBases, a.isReference());
                alleles.add(trimmedAllele);
                originalToTrimmedAlleleMap.put(a, trimmedAllele);
            }
        }

        // now we can recreate new genotypes with trimmed alleles
        final AlleleMapper alleleMapper = new AlleleMapper(originalToTrimmedAlleleMap);
        final GenotypesContext genotypes = updateGenotypesWithMappedAlleles(inputVC.getGenotypes(), alleleMapper);

        final int start = inputVC.getStart() + (fwdTrimEnd + 1);
        final VariantContextBuilder builder = new VariantContextBuilder(inputVC);
        builder.start(start);
        builder.stop(start + alleles.get(0).length() - 1);
        builder.alleles(alleles);
        builder.genotypes(genotypes);
        return builder.make();
    }

    protected static GenotypesContext updateGenotypesWithMappedAlleles(final GenotypesContext originalGenotypes, final AlleleMapper alleleMapper) {
        final GenotypesContext updatedGenotypes = GenotypesContext.create(originalGenotypes.size());

        for ( final Genotype genotype : originalGenotypes ) {
            final List<Allele> updatedAlleles = alleleMapper.remap(genotype.getAlleles());
            updatedGenotypes.add(new GenotypeBuilder(genotype).alleles(updatedAlleles).make());
        }

        return updatedGenotypes;
    }

    public static int computeReverseClipping(final List<Allele> unclippedAlleles, final byte[] ref) {
        int clipping = 0;
        boolean stillClipping = true;

        while ( stillClipping ) {
            for ( final Allele a : unclippedAlleles ) {
                if ( a.isSymbolic() )
                    continue;

                // we need to ensure that we don't reverse clip out all of the bases from an allele because we then will have the wrong
                // position set for the VariantContext (although it's okay to forward clip it all out, because the position will be fine).
                if ( a.length() - clipping == 0 )
                    return clipping - 1;

                if ( a.length() - clipping <= 0 || a.length() == 0 ) {
                    stillClipping = false;
                }
                else if ( ref.length == clipping ) {
                    return -1;
                }
                else if ( a.getBases()[a.length()-clipping-1] != ref[ref.length-clipping-1] ) {
                    stillClipping = false;
                }
            }
            if ( stillClipping )
                clipping++;
        }

        return clipping;
    }

    /**
     * Clip out any unnecessary bases off the front of the alleles
     *
     * The VCF spec represents alleles as block substitutions, replacing AC with A for a
     * 1 bp deletion of the C.  However, it's possible that we'd end up with alleles that
     * contain extra bases on the left, such as GAC/GA to represent the same 1 bp deletion.
     * This routine finds an offset among all alleles that can be safely trimmed
     * off the left of each allele and still represent the same block substitution.
     *
     * A/C => A/C
     * AC/A => AC/A
     * ACC/AC => CC/C
     * AGT/CAT => AGT/CAT
     * <DEL>/C => <DEL>/C
     *
     * @param unclippedAlleles a non-null list of alleles that we want to clip
     * @return the offset into the alleles where we can safely clip, inclusive, or
     *   -1 if no clipping is tolerated.  So, if the result is 0, then we can remove
     *   the first base of every allele.  If the result is 1, we can remove the
     *   second base.
     */
    public static int computeForwardClipping(final List<Allele> unclippedAlleles) {
        // cannot clip unless there's at least 1 alt allele
        if ( unclippedAlleles.size() <= 1 )
            return -1;

        // we cannot forward clip any set of alleles containing a symbolic allele
        int minAlleleLength = Integer.MAX_VALUE;
        for ( final Allele a : unclippedAlleles ) {
            if ( a.isSymbolic() )
                return -1;
            minAlleleLength = Math.min(minAlleleLength, a.length());
        }

        final byte[] firstAlleleBases = unclippedAlleles.get(0).getBases();
        int indexOflastSharedBase = -1;

        // the -1 to the stop is that we can never clip off the right most base
        for ( int i = 0; i < minAlleleLength - 1; i++) {
            final byte base = firstAlleleBases[i];

            for ( final Allele allele : unclippedAlleles ) {
                if ( allele.getBases()[i] != base )
                    return indexOflastSharedBase;
            }

            indexOflastSharedBase = i;
        }

        return indexOflastSharedBase;
    }

    private static Map<String, Object> subsetAttributes(final CommonInfo igc, final Collection<String> keysToPreserve) {
        Map<String, Object> attributes = new HashMap<>(keysToPreserve.size());
        for ( final String key : keysToPreserve  ) {
            if ( igc.hasAttribute(key) )
                attributes.put(key, igc.getAttribute(key));
        }
        return attributes;
    }

    protected static class AlleleMapper {
        private VariantContext vc = null;
        private Map<Allele, Allele> map = null;
        public AlleleMapper(VariantContext vc)          { this.vc = vc; }
        public AlleleMapper(Map<Allele, Allele> map)    { this.map = map; }
        public boolean needsRemapping()                 { return this.map != null; }
        public Collection<Allele> values()              { return map != null ? map.values() : vc.getAlleles(); }
        public Allele remap(Allele a)                   { return map != null && map.containsKey(a) ? map.get(a) : a; }

        public List<Allele> remap(List<Allele> as) {
            List<Allele> newAs = as.stream()
                    .map(this::remap)
                    .collect(Collectors.toList());
            //System.out.printf("  Remapping %s => %s%n", a, remap(a));
            return newAs;
        }

    }

    private static class CompareByPriority implements Comparator<VariantContext>, Serializable {
        private static final long serialVersionUID = 0L;

        List<String> priorityListOfVCs;
        public CompareByPriority(List<String> priorityListOfVCs) {
            this.priorityListOfVCs = priorityListOfVCs;
        }

        private int getIndex(VariantContext vc) {
            int i = priorityListOfVCs.indexOf(vc.getSource());
            if ( i == -1 ) throw new IllegalArgumentException("Priority list " + priorityListOfVCs + " doesn't contain variant context " + vc.getSource());
            return i;
        }

        public int compare(VariantContext vc1, VariantContext vc2) {
            return Integer.valueOf(getIndex(vc1)).compareTo(getIndex(vc2));
        }
    }

    /**
     * For testing purposes only.  Create a site-only VariantContext at contig:start containing alleles
     *
     * @param name the name of the VC
     * @param contig the contig for the VC
     * @param start the start of the VC
     * @param alleleStrings a non-null, non-empty list of strings for the alleles.  The first will be the ref allele, and others the
     *                      alt.  Will compute the stop of the VC from the length of the reference allele
     * @return a non-null VariantContext
     */
    public static VariantContext makeFromAlleles(final String name, final String contig, final int start, final List<String> alleleStrings) {
        if ( alleleStrings == null || alleleStrings.isEmpty() )
            throw new IllegalArgumentException("alleleStrings must be non-empty, non-null list");

        final List<Allele> alleles = new LinkedList<>();
        final int length = alleleStrings.get(0).length();

        boolean first = true;
        for ( final String alleleString : alleleStrings ) {
            alleles.add(Allele.create(alleleString, first));
            first = false;
        }
        return new VariantContextBuilder(name, contig, start, start+length-1, alleles).make();
    }

    /**
     * Splits the alleles for the provided variant context into its primitive parts.
     * Requires that the input VC be bi-allelic, so calling methods should first call splitVariantContextToBiallelics() if needed.
     * Currently works only for MNPs.
     *
     * @param vc  the non-null VC to split
     * @return a non-empty list of VCs split into primitive parts or the original VC otherwise
     */
    public static List<VariantContext> splitIntoPrimitiveAlleles(final VariantContext vc) {
        if ( vc == null )
            throw new IllegalArgumentException("Trying to break a null Variant Context into primitive parts");

        if ( !vc.isBiallelic() )
            throw new IllegalArgumentException("Trying to break a multi-allelic Variant Context into primitive parts");

        // currently only works for MNPs
        if ( !vc.isMNP() )
            return Arrays.asList(vc);

        final byte[] ref = vc.getReference().getBases();
        final byte[] alt = vc.getAlternateAllele(0).getBases();

        if ( ref.length != alt.length )
            throw new IllegalStateException("ref and alt alleles for MNP have different lengths");

        final List<VariantContext> result = new ArrayList<>(ref.length);

        for ( int i = 0; i < ref.length; i++ ) {

            // if the ref and alt bases are different at a given position, create a new SNP record (otherwise do nothing)
            if ( ref[i] != alt[i] ) {

                // create the ref and alt SNP alleles
                final Allele newRefAllele = Allele.create(ref[i], true);
                final Allele newAltAllele = Allele.create(alt[i], false);

                // create a new VariantContext with the new SNP alleles
                final VariantContextBuilder newVC = new VariantContextBuilder(vc).start(vc.getStart() + i).stop(vc.getStart() + i).alleles(Arrays.asList(newRefAllele, newAltAllele));

                // create new genotypes with updated alleles
                final Map<Allele, Allele> alleleMap = new HashMap<>();
                alleleMap.put(vc.getReference(), newRefAllele);
                alleleMap.put(vc.getAlternateAllele(0), newAltAllele);
                final GenotypesContext newGenotypes = updateGenotypesWithMappedAlleles(vc.getGenotypes(), new AlleleMapper(alleleMap));

                result.add(newVC.genotypes(newGenotypes).make());
            }
        }

        if ( result.isEmpty() )
            result.add(vc);

        return result;
    }

    /**
     * Are vc1 and 2 equal including their position and alleles?
     * @param vc1 non-null VariantContext
     * @param vc2 non-null VariantContext
     * @return true if vc1 and vc2 are equal, false otherwise
     */
    public static boolean equalSites(final VariantContext vc1, final VariantContext vc2) {
        if ( vc1 == null ) throw new IllegalArgumentException("vc1 cannot be null");
        if ( vc2 == null ) throw new IllegalArgumentException("vc2 cannot be null");

        if ( vc1.getStart() != vc2.getStart() ) return false;
        if ( vc1.getEnd() != vc2.getEnd() ) return false;
        if ( !vc1.getContig().equals(vc2.getContig())) return false;
        if ( ! vc1.getAlleles().equals(vc2.getAlleles()) ) return false;
        return true;
    }

    /**
     * Returns the absolute 0-based index of an allele.
     *
     * <p/>
     * If the allele is equal to the reference, the result is 0, if it equal to the first alternative the result is 1
     * and so forth.
     * <p/>
     * Therefore if you want the 0-based index within the alternative alleles you need to do the following:
     *
     * <p/>
     * You can indicate whether the Java object reference comparator {@code ==} can be safelly used by setting {@code useEquals} to {@code false}.
     *
     * @param vc the target variant context.
     * @param allele the target allele.
     * @param ignoreRefState whether the reference states of the allele is important at all. Has no effect if {@code useEquals} is {@code false}.
     * @param considerRefAllele whether the reference allele should be considered. You should set it to {@code false} if you are only interested in alternative alleles.
     * @param useEquals whether equal method should be used in the search: {@link Allele#equals(Allele,boolean)}.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null}.
     * @return {@code -1} if there is no such allele that satify those criteria, a value between 0 and {@link VariantContext#getNAlleles()} {@code -1} otherwise.
     */
    public static int indexOfAllele(final VariantContext vc, final Allele allele, final boolean ignoreRefState, final boolean considerRefAllele, final boolean useEquals) {
        if (allele == null) throw new IllegalArgumentException();
        return useEquals ? indexOfEqualAllele(vc,allele,ignoreRefState,considerRefAllele) : indexOfSameAllele(vc,allele,considerRefAllele);
    }

    /**
     * Returns the relative 0-based index of an alternative allele.
     * <p/>
     * The the query allele is the same as the first alternative allele, the result is 0,
     * if it is equal to the second 1 and so forth.
     *
     *
     * <p/>
     * Notice that the ref-status of the query {@code allele} is ignored.
     *
     * @param vc the target variant context.
     * @param allele the query allele.
     * @param useEquals  whether equal method should be used in the search: {@link Allele#equals(Allele,boolean)}.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null}.
     *
     * @return {@code -1} if there is no such allele that satify those criteria, a value between 0 and the number
     *  of alternative alleles - 1.
     */
    public static int indexOfAltAllele(final VariantContext vc, final Allele allele, final boolean useEquals) {
        final int absoluteIndex = indexOfAllele(vc,allele,true,false,useEquals);
        return absoluteIndex == -1 ? -1 : absoluteIndex - 1;
    }

    // Impements index search using equals.
    private static int indexOfEqualAllele(final VariantContext vc, final Allele allele, final boolean ignoreRefState,
                                          final boolean considerRefAllele) {
        int i = 0;
        for (final Allele a : vc.getAlleles())
            if (a.equals(allele,ignoreRefState))
                return i == 0 ? (considerRefAllele ? 0 : -1) : i;
            else
                i++;
        return -1;
    }

    // Implements index search using ==.
    private static int indexOfSameAllele(final VariantContext vc, final Allele allele, final boolean considerRefAllele) {
        int i = 0;

        for (final Allele a : vc.getAlleles())
            if (a == allele)
                return i == 0 ? (considerRefAllele ? 0 : -1) : i;
            else
                i++;

        return -1;
    }
}
