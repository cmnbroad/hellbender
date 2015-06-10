package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class LocalAssemblyEngineUnitTest extends BaseTest {
    private GenomeLocParser genomeLocParser;
    private IndexedFastaSequenceFile seq;
    private SAMFileHeader header;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        seq = new CachingIndexedFastaSequenceFile(new File(hg19_chr1_1M_Reference));
        genomeLocParser = new GenomeLocParser(seq);
        header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
    }

    @DataProvider(name = "AssembleIntervalsData")
    public Object[][] makeAssembleIntervalsData() {
        List<Object[]> tests = new ArrayList<>();

        final String contig = "20";
        final int start = 10000000;
        final int end   = 10100000;
        final int windowSize = 100;
        final int stepSize = 200;
        final int nReadsToUse = 5;

        for ( int startI = start; startI < end; startI += stepSize) {
            final int endI = startI + windowSize;
            final GenomeLoc refLoc = genomeLocParser.createGenomeLoc(contig, startI, endI);
            tests.add(new Object[]{new ReadThreadingAssembler(), refLoc, nReadsToUse});
        }

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "AssembleIntervalsWithVariantData")
    public Object[][] makeAssembleIntervalsWithVariantData() {
        List<Object[]> tests = new ArrayList<>();

        final String contig = "20";
        final int start = 10000000;
        final int end   = 10001000;
        final int windowSize = 100;
        final int stepSize = 200;
        final int variantStepSize = 1;
        final int nReadsToUse = 5;

        for ( int startI = start; startI < end; startI += stepSize) {
            final int endI = startI + windowSize;
            final GenomeLoc refLoc = genomeLocParser.createGenomeLoc(contig, startI, endI);
            for ( int variantStart = windowSize / 2 - 10; variantStart < windowSize / 2 + 10; variantStart += variantStepSize ) {
                tests.add(new Object[]{new ReadThreadingAssembler(), refLoc, nReadsToUse, variantStart});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AssembleIntervalsData")
    public void testAssembleRef(final ReadThreadingAssembler assembler, final GenomeLoc loc, final int nReadsToUse) {
        final byte[] refBases = seq.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();

        final List<SAMRecord> reads = new LinkedList<>();
        for ( int i = 0; i < nReadsToUse; i++ ) {
            final byte[] bases = refBases.clone();
            final byte[] quals = Utils.dupBytes((byte) 30, refBases.length);
            final String cigar = refBases.length + "M";
            final SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, loc.getContig(), loc.getContigIndex(), loc.getStart(), bases, quals, cigar);
            reads.add(read);
        }

        // TODO -- generalize to all assemblers
        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final List<Haplotype> haplotypes = assemble(assembler, refBases, loc, reads);
        Assert.assertEquals(haplotypes, Collections.singletonList(refHaplotype));
    }

    @Test(dataProvider = "AssembleIntervalsWithVariantData")
    public void testAssembleRefAndSNP(final ReadThreadingAssembler assembler, final GenomeLoc loc, final int nReadsToUse, final int variantSite) {
        final byte[] refBases = seq.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
        final Allele refBase = Allele.create(refBases[variantSite], true);
        final Allele altBase = Allele.create((byte) (refBase.getBases()[0] == 'A' ? 'C' : 'A'), false);
        final VariantContextBuilder vcb = new VariantContextBuilder("x", loc.getContig(), variantSite, variantSite, Arrays.asList(refBase, altBase));
        testAssemblyWithVariant(assembler, refBases, loc, nReadsToUse, vcb.make());
    }

    @Test(dataProvider = "AssembleIntervalsWithVariantData")
    public void testAssembleRefAndDeletion(final ReadThreadingAssembler assembler, final GenomeLoc loc, final int nReadsToUse, final int variantSite) {
        final byte[] refBases = seq.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
        for ( int deletionLength = 1; deletionLength < 10; deletionLength++ ) {
            final Allele refBase = Allele.create(new String(refBases).substring(variantSite, variantSite + deletionLength + 1), true);
            final Allele altBase = Allele.create(refBase.getBases()[0], false);
            final VariantContextBuilder vcb = new VariantContextBuilder("x", loc.getContig(), variantSite, variantSite + deletionLength, Arrays.asList(refBase, altBase));
            testAssemblyWithVariant(assembler, refBases, loc, nReadsToUse, vcb.make());
        }
    }

    @Test(dataProvider = "AssembleIntervalsWithVariantData")
    public void testAssembleRefAndInsertion(final ReadThreadingAssembler assembler, final GenomeLoc loc, final int nReadsToUse, final int variantSite) {
        final byte[] refBases = seq.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
        for ( int insertionLength = 1; insertionLength < 10; insertionLength++ ) {
            final Allele refBase = Allele.create(refBases[variantSite], false);
            final Allele altBase = Allele.create(new String(refBases).substring(variantSite, variantSite + insertionLength + 1), true);
            final VariantContextBuilder vcb = new VariantContextBuilder("x", loc.getContig(), variantSite, variantSite + insertionLength, Arrays.asList(refBase, altBase));
            testAssemblyWithVariant(assembler, refBases, loc, nReadsToUse, vcb.make());
        }
    }

    private void testAssemblyWithVariant(final ReadThreadingAssembler assembler, final byte[] refBases, final GenomeLoc loc, final int nReadsToUse, final VariantContext site) {
        final String preRef = new String(refBases).substring(0, site.getStart());
        final String postRef = new String(refBases).substring(site.getEnd() + 1, refBases.length);
        final byte[] altBases = (preRef + site.getAlternateAllele(0).getBaseString() + postRef).getBytes();

//        logger.warn("ref " + new String(refBases));
//        logger.warn("alt " + new String(altBases));

        final List<SAMRecord> reads = new LinkedList<>();
        for ( int i = 0; i < nReadsToUse; i++ ) {
            final byte[] bases = altBases.clone();
            final byte[] quals = Utils.dupBytes((byte) 30, altBases.length);
            final String cigar = altBases.length + "M";
            final SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, loc.getContig(), loc.getContigIndex(), loc.getStart(), bases, quals, cigar);
            reads.add(read);
        }

        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final Haplotype altHaplotype = new Haplotype(altBases, false);
        final List<Haplotype> haplotypes = assemble(assembler, refBases, loc, reads);
        Assert.assertEquals(haplotypes, Arrays.asList(refHaplotype, altHaplotype));
    }


    private List<Haplotype> assemble(final ReadThreadingAssembler assembler, final byte[] refBases, final GenomeLoc loc, final List<SAMRecord> reads) {
        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final Cigar c = new Cigar();
        c.add(new CigarElement(refHaplotype.getBases().length, CigarOperator.M));
        refHaplotype.setCigar(c);

        final AssemblyRegion activeRegion = new AssemblyRegion(loc, null, true, genomeLocParser, 0);
        activeRegion.addAll(reads);
//        logger.warn("Assembling " + activeRegion + " with " + engine);
        final AssemblyResultSet assemblyResultSet =  assembler.runLocalAssembly(activeRegion, refHaplotype, refBases, loc, Collections.<VariantContext>emptyList(), null);
        return assemblyResultSet.getHaplotypeList();
    }

    @DataProvider(name = "SimpleAssemblyTestData")
    public Object[][] makeSimpleAssemblyTestData() {
        List<Object[]> tests = new ArrayList<>();

        final String contig  = "20";
        final int start      = 10000000;
        final int windowSize = 200;
        final int end        = start + windowSize;

        final int excludeVariantsWithinXbp = 25; // TODO -- decrease to zero when the edge calling problem is fixed

        final String ref = new String(seq.getSubsequenceAt(contig, start, end).getBases());
        final GenomeLoc refLoc = genomeLocParser.createGenomeLoc(contig, start, end);

            for ( int snpPos = 0; snpPos < windowSize; snpPos++) {
                if ( snpPos > excludeVariantsWithinXbp && (windowSize - snpPos) >= excludeVariantsWithinXbp ) {
                    final byte[] altBases = ref.getBytes();
                    altBases[snpPos] = altBases[snpPos] == 'A' ? (byte)'C' : (byte)'A';
                    final String alt = new String(altBases);
                    tests.add(new Object[]{"SNP at " + snpPos, new ReadThreadingAssembler(), refLoc, ref, alt});
                }
            }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "SimpleAssemblyTestData")
    public void testSimpleAssembly(final String name, final ReadThreadingAssembler assembler, final GenomeLoc loc, final String ref, final String alt) {
        final byte[] refBases = ref.getBytes();
        final byte[] altBases = alt.getBytes();

        final List<SAMRecord> reads = new LinkedList<>();
        for ( int i = 0; i < 20; i++ ) {
            final byte[] bases = altBases.clone();
            final byte[] quals = Utils.dupBytes((byte) 30, altBases.length);
            final String cigar = altBases.length + "M";
            final SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, loc.getContig(), loc.getContigIndex(), loc.getStart(), bases, quals, cigar);
            reads.add(read);
        }

        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final Haplotype altHaplotype = new Haplotype(altBases, false);
        final List<Haplotype> haplotypes = assemble(assembler, refBases, loc, reads);
        Assert.assertTrue(haplotypes.size() > 0, "Failed to find ref haplotype");
        Assert.assertEquals(haplotypes.get(0), refHaplotype);

        Assert.assertEquals(haplotypes.size(), 2, "Failed to find single alt haplotype");
        Assert.assertEquals(haplotypes.get(1), altHaplotype);
    }
}
