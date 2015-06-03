package org.broadinstitute.hellbender.utils.read;

import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import htsjdk.samtools.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

public final class ArtificialReadUtils {
    public static final int DEFAULT_READ_LENGTH = 50;

    private static final String DEFAULT_READ_GROUP_PREFIX = "ReadGroup";
    private static final String DEFAULT_PLATFORM_UNIT_PREFIX = "Lane";
    private static final String DEFAULT_PLATFORM_PREFIX = "Platform";
    private static final String DEFAULT_SAMPLE_NAME = "SampleX";

    /**
     * Creates an artificial sam header, matching the parameters, chromosomes which will be labeled chr1, chr2, etc
     *
     * @param numberOfChromosomes the number of chromosomes to create
     * @param startingChromosome  the starting number for the chromosome (most likely set to 1)
     * @param chromosomeSize      the length of each chromosome
     * @return
     */
    public static SAMFileHeader createArtificialSamHeader(int numberOfChromosomes, int startingChromosome, int chromosomeSize) {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(htsjdk.samtools.SAMFileHeader.SortOrder.coordinate);
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        // make up some sequence records
        for (int x = startingChromosome; x < startingChromosome + numberOfChromosomes; x++) {
            SAMSequenceRecord rec = new SAMSequenceRecord( Integer.toString(x), chromosomeSize /* size */);
            rec.setSequenceLength(chromosomeSize);
            dict.addSequence(rec);
        }
        header.setSequenceDictionary(dict);
        return header;
    }

    /**
     * Creates an artificial sam header, matching the parameters, chromosomes which will be labeled chr1, chr2, etc
     * It also adds read groups.
     *
     * @param numberOfChromosomes the number of chromosomes to create
     * @param startingChromosome  the starting number for the chromosome (most likely set to 1)
     * @param chromosomeSize      the length of each chromosome
     * @param groupCount          the number of groups to make
     */
    public static SAMFileHeader createArtificialSamHeaderWithGroups(int numberOfChromosomes, int startingChromosome, int chromosomeSize, int groupCount) {
        final SAMFileHeader header = createArtificialSamHeader(numberOfChromosomes, startingChromosome, chromosomeSize);

        final List<SAMReadGroupRecord> readGroups = new ArrayList<>();
        for (int i = 0; i < groupCount; i++) {
            SAMReadGroupRecord rec = new SAMReadGroupRecord(DEFAULT_READ_GROUP_PREFIX + i);
            rec.setSample(DEFAULT_SAMPLE_NAME);
            readGroups.add(rec);
        }
        header.setReadGroups(readGroups);

        for (int i = 0; i < groupCount; i++) {
            final SAMReadGroupRecord groupRecord = header.getReadGroup(readGroups.get(i).getId());
            groupRecord.setPlatform(DEFAULT_PLATFORM_PREFIX + ((i % 2) + 1));
            groupRecord.setPlatformUnit(DEFAULT_PLATFORM_UNIT_PREFIX + ((i % 3) + 1));
        }
        return header;
    }

    /**
     * Creates an artificial sam header with standard test parameters
     *
     * @return the sam header
     */
    public static SAMFileHeader createArtificialSamHeader() {
        return createArtificialSamHeader(1, 1, 1000000);
    }

    public static SAMFileHeader createArtificialSamHeaderWithReadGroup( final SAMReadGroupRecord readGroup ) {
        final SAMFileHeader header = createArtificialSamHeader();
        header.addReadGroup(readGroup);
        return header;
    }

    private static SAMRecord createArtificialSAMRecord(SAMFileHeader header, String name, int refIndex, int alignmentStart, int length) {
        if ((refIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart != SAMRecord.NO_ALIGNMENT_START) ||
                (refIndex != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart == SAMRecord.NO_ALIGNMENT_START))
            throw new IllegalArgumentException("Invalid alignment start for artificial read, start = " + alignmentStart);
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceIndex(refIndex);
        record.setAlignmentStart(alignmentStart);
        List<CigarElement> elements = new ArrayList<>();
        elements.add(new CigarElement(length, CigarOperator.characterToEnum('M')));
        record.setCigar(new Cigar(elements));
        record.setProperPairFlag(false);

        // our reads and quals are all 'A's by default
        byte[] c = new byte[length];
        byte[] q = new byte[length];
        for (int x = 0; x < length; x++)
            c[x] = q[x] = 'A';
        record.setReadBases(c);
        record.setBaseQualities(q);

        if (refIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
            record.setReadUnmappedFlag(true);
        }

        return record;
    }
    /**
     * Create an artificial read based on the parameters.  The cigar string will be *M, where * is the length of the read
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param length         the length of the read
     * @return the artificial read
     */
    public static MutableGATKRead createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, int length) {
        if ((refIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart != SAMRecord.NO_ALIGNMENT_START) ||
                (refIndex != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart == SAMRecord.NO_ALIGNMENT_START))
            throw new IllegalArgumentException("Invalid alignment start for artificial read, start = " + alignmentStart);
        SAMRecord record = createArtificialSAMRecord(header, name, refIndex, alignmentStart, length);
        return new SAMRecordToGATKReadAdapter(record);
    }

    /**
     * Create an artificial read based on the parameters.  The cigar string will be *M, where * is the length of the read
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param bases          the sequence of the read
     * @param qual           the qualities of the read
     * @return the artificial read
     */
    public static MutableGATKRead createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, byte[] bases, byte[] qual) {
        if (bases.length != qual.length) {
            throw new IllegalArgumentException("Passed in read string is different length then the quality array");
        }
        MutableGATKRead rec = createArtificialRead(header, name, refIndex, alignmentStart, bases.length);
        rec.setBases(Arrays.copyOf(bases, bases.length));
        rec.setBaseQualities(Arrays.copyOf(qual, qual.length));
        rec.setAttribute(SAMTag.PG.name(), new SAMReadGroupRecord("x").getId());
        if (refIndex == -1) {
            rec.setIsUnmapped();
        }

        return rec;
    }

    /**
     * Create an artificial read based on the parameters
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param bases          the sequence of the read
     * @param qual           the qualities of the read
     * @param cigar          the cigar string of the read
     * @return the artificial read
     */
    public static MutableGATKRead createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, byte[] bases, byte[] qual, String cigar) {
        MutableGATKRead rec = createArtificialRead(header, name, refIndex, alignmentStart, bases, qual);
        rec.setCigar(cigar);
        return rec;
    }

    /**
     * Create an artificial read with the following default parameters :
     * header:
     * numberOfChromosomes = 1
     * startingChromosome = 1
     * chromosomeSize = 1000000
     * read:
     * name = "default_read"
     * refIndex = 0
     * alignmentStart = 10000
     *
     * @param header SAM header for the read
     * @param bases the sequence of the read
     * @param qual  the qualities of the read
     * @param cigar the cigar string of the read
     * @return the artificial read
     */
    public static MutableGATKRead createArtificialRead(final SAMFileHeader header, final byte[] bases, final byte[] qual, final String cigar) {
        return createArtificialRead(header, "default_read", 0, 10000, bases, qual, cigar);
    }

    public static MutableGATKRead createArtificialRead(final byte[] bases, final byte[] qual, final String cigar) {
        SAMFileHeader header = createArtificialSamHeader();
        return createArtificialRead(header, "default_read", 0, 10000, bases, qual, cigar);
    }

    public static MutableGATKRead createArtificialRead(final SAMFileHeader header, final Cigar cigar) {
        int length = cigar.getReadLength();
        byte [] base = {'A'};
        byte [] qual = {30};
        byte [] bases = Utils.arrayFromArrayWithLength(base, length);
        byte [] quals = Utils.arrayFromArrayWithLength(qual, length);
        return createArtificialRead(header, "default_read", 0, 10000, bases, quals, cigar.toString());
    }

    public static MutableGATKRead createArtificialRead(final Cigar cigar) {
        final SAMFileHeader header = createArtificialSamHeader();
        return createArtificialRead(header, cigar);
    }

    public static List<MutableGATKRead> createPair(SAMFileHeader header, String name, int readLen, int leftStart, int rightStart, boolean leftIsFirst, boolean leftIsNegative) {
        MutableGATKRead left = createArtificialRead(header, name, 0, leftStart, readLen);
        MutableGATKRead right = createArtificialRead(header, name, 0, rightStart, readLen);

        left.setIsPaired(true);
        right.setIsPaired(true);

        left.setIsProperlyPaired(true);
        right.setIsProperlyPaired(true);

        left.setIsFirstOfPair(leftIsFirst);
        right.setIsFirstOfPair(!leftIsFirst);

        left.setIsReverseStrand(leftIsNegative);
        left.setMateIsReverseStrand(!leftIsNegative);
        right.setIsReverseStrand(!leftIsNegative);
        right.setMateIsReverseStrand(leftIsNegative);

        left.setMatePosition(header.getSequence(0).getSequenceName(), right.getStart());
        right.setMatePosition(header.getSequence(0).getSequenceName(), left.getStart());

        int isize = rightStart + readLen - leftStart;
        left.setFragmentLength(isize);
        right.setFragmentLength(-isize);

        return Arrays.asList(left, right);
    }

    public static MutableGATKRead createRandomRead(SAMFileHeader header, int length) {
        List<CigarElement> cigarElements = new LinkedList<>();
        cigarElements.add(new CigarElement(length, CigarOperator.M));
        Cigar cigar = new Cigar(cigarElements);
        return createArtificialRead(header, cigar);
    }

    public static MutableGATKRead createRandomRead(int length) {
        SAMFileHeader header = createArtificialSamHeader();
        return createRandomRead(header, length);
    }

    public static MutableGATKRead createRandomRead(int start, int length) {
        List<CigarElement> cigarElements = new LinkedList<>();
        cigarElements.add(new CigarElement(length, CigarOperator.M));
        Cigar cigar = new Cigar(cigarElements);
        MutableGATKRead artificialRead = ArtificialReadUtils.createArtificialRead(cigar);
        artificialRead.setPosition(artificialRead.getContig(), start);
        return artificialRead;
    }

    public static MutableGATKRead createRandomRead(int start, int length, int UUIDSeed) {
        List<CigarElement> cigarElements = new LinkedList<>();
        cigarElements.add(new CigarElement(length, CigarOperator.M));
        Cigar cigar = new Cigar(cigarElements);
        MutableGATKRead artificialRead = ArtificialReadUtils.createArtificialRead(cigar);
        artificialRead.setPosition(artificialRead.getContig(), start);
        // We know that a SAMRecordToGATKReadAdapter is backing the MutableGATKRead at this point.
        SAMRecordToGATKReadAdapter sam = (SAMRecordToGATKReadAdapter) artificialRead;
        sam.setUUID(new UUID(UUIDSeed, UUIDSeed));
        return artificialRead;
    }

    public static MutableGATKRead createRandomGoogleRead(int start, int length, int UUIDSeed) {
        List<CigarElement> cigarElements = new LinkedList<>();
        cigarElements.add(new CigarElement(length, CigarOperator.M));
        Cigar cigar = new Cigar(cigarElements);
        final SAMFileHeader header = createArtificialSamHeader();

        //int length = cigar.getReadLength();

        byte [] base = {'A'};
        byte [] qual = {30};
        byte [] bases = Utils.arrayFromArrayWithLength(base, length);
        byte [] quals = Utils.arrayFromArrayWithLength(qual, length);
        String name = "default_read";
        int refIndex = 0;
        int alignmentStart = 10000;
        //createArtificialRead(header, 0, 10000, bases, quals, cigar.toString());
        //MutableGATKRead rec = createArtificialRead(header, name, refIndex, alignmentStart, bases, qual);
        if (bases.length != qual.length) {
            throw new IllegalArgumentException("Passed in read string is different length then the quality array");
        }

        //MutableGATKRead rec = createArtificialRead(header, name, refIndex, alignmentStart, bases.length);
        //if (refIndex == -1) {
        //    rec.setIsUnmapped();
        //}
        //if ((refIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart != SAMRecord.NO_ALIGNMENT_START) ||
        //        (refIndex != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart == SAMRecord.NO_ALIGNMENT_START))
        //    throw new IllegalArgumentException("Invalid alignment start for artificial read, start = " + alignmentStart);

        SAMRecord rec = createArtificialSAMRecord(header, name, refIndex, alignmentStart, length);
        rec.setCigar(cigar);

        GoogleGenomicsReadToGATKReadAdapter record = new GoogleGenomicsReadToGATKReadAdapter(ReadConverter.makeRead(rec));
        record.setBases(Arrays.copyOf(bases, bases.length));
        record.setBaseQualities(Arrays.copyOf(qual, qual.length));
        record.setAttribute(SAMTag.PG.name(), new SAMReadGroupRecord("x").getId());
        return record;
    }

    public static MutableGATKRead createRandomRead(SAMFileHeader header, int length, boolean allowNs) {
        byte[] quals = createRandomReadQuals(length);
        byte[] bbases = createRandomReadBases(length, allowNs);
        return createArtificialRead(bbases, quals, bbases.length + "M");
    }

    public static MutableGATKRead createRandomRead(int length, boolean allowNs) {
        SAMFileHeader header = createArtificialSamHeader();
        return createRandomRead(header, length, allowNs);
    }

    /**
     * Create random read qualities
     *
     * @param length the length of the read
     * @return an array with randomized base qualities between 0 and 50
     */
    public static byte[] createRandomReadQuals(int length) {
        Random random = Utils.getRandomGenerator();
        byte[] quals = new byte[length];
        for (int i = 0; i < length; i++)
            quals[i] = (byte) random.nextInt(50);
        return quals;
    }

    /**
     * Create random read qualities
     *
     * @param length  the length of the read
     * @param allowNs whether or not to allow N's in the read
     * @return an array with randomized bases (A-N) with equal probability
     */
    public static byte[] createRandomReadBases(int length, boolean allowNs) {
        Random random = Utils.getRandomGenerator();
        int numberOfBases = allowNs ? 5 : 4;
        byte[] bases = new byte[length];
        for (int i = 0; i < length; i++) {
            switch (random.nextInt(numberOfBases)) {
                case 0:
                    bases[i] = 'A';
                    break;
                case 1:
                    bases[i] = 'C';
                    break;
                case 2:
                    bases[i] = 'G';
                    break;
                case 3:
                    bases[i] = 'T';
                    break;
                case 4:
                    bases[i] = 'N';
                    break;
                default:
                    throw new GATKException("Something went wrong, this is just impossible");
            }
        }
        return bases;
    }
    /**
     * create an iterator containing the specified read piles
     *
     * @param startingChr the chromosome (reference ID) to start from
     * @param endingChr   the id to end with
     * @param readCount   the number of reads per chromosome
     * @return iterator representing the specified amount of fake data
     */
    public static ArtificialReadQueryIterator mappedReadIterator(int startingChr, int endingChr, int readCount) {
        SAMFileHeader header = createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialReadQueryIterator(startingChr, endingChr, readCount, 0, header);
    }

    /**
     * create an iterator containing the specified read piles
     *
     * @param startingChr       the chromosome (reference ID) to start from
     * @param endingChr         the id to end with
     * @param readCount         the number of reads per chromosome
     * @param unmappedReadCount the count of unmapped reads to place at the end of the iterator, like in a sorted bam file
     * @return iterator representing the specified amount of fake data
     */
    public static ArtificialReadQueryIterator mappedAndUnmappedReadIterator(int startingChr, int endingChr, int readCount, int unmappedReadCount) {
        SAMFileHeader header = createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialReadQueryIterator(startingChr, endingChr, readCount, unmappedReadCount, header);
    }

    /**
     * Creates an artificial sam header based on the sequence dictionary dict
     *
     * @return a new sam header
     */
    public static SAMFileHeader createArtificialSamHeader(final SAMSequenceDictionary dict) {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(htsjdk.samtools.SAMFileHeader.SortOrder.coordinate);
        header.setSequenceDictionary(dict);
        return header;
    }
}
