package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import java.util.Arrays;

/**
 * Fast wrapper for byte[] kmers
 *
 * This objects has several important features that make it better than using a raw byte[] for a kmer:
 *
 * -- Can create kmer from a range of a larger byte[], allowing us to avoid Array.copyOfRange
 * -- Fast equals and hashcode methods
 * -- can get actual byte[] of the kmer, even if it's from a larger byte[], and this operation
 *    only does the work of that operation once, updating its internal state
 */
public class Kmer {
    // this values may be updated in the course of interacting with this kmer
    byte[] bases;
    int start;

    // two constants
    final int length;
    final int hash;

    /**
     * Create a new kmer using all bases in kmer
     * @param kmer a non-null byte[]
     */
    public Kmer(final byte[] kmer) {
        this(kmer, 0, kmer.length);
    }

    /**
     * Create a new kmer based on the string kmer
     *
     * This is not a good method to use for performance
     *
     * @param kmer the bases as a string
     */
    public Kmer(final String kmer) {
        this(kmer.getBytes());
    }

    /**
     * Create a new kmer backed by the bases in bases, spanning start -> start + length
     *
     * Under no circumstances can bases be modified anywhere in the client code.  This does not make a copy
     * of bases for performance reasons
     *
     * @param bases an array of bases
     * @param start the start of the kmer in bases, must be >= 0 and < bases.length
     * @param length the length of the kmer.  Must be >= 0 and start + length < bases.length
     */
    public Kmer(final byte[] bases, final int start, final int length) {
        if ( bases == null ) throw new IllegalArgumentException("bases cannot be null");
        if ( start < 0 ) throw new IllegalArgumentException("start must be >= 0 but got " + start);
        if ( length < 0 ) throw new IllegalArgumentException("length must be >= 0 but got " + length);
        if ( (start + length) > bases.length ) throw new IllegalArgumentException("start + length " + (start + length) + " must be <= bases.length " + bases.length + " but got " + start + " with length " + length);
        this.bases = bases;
        this.start = start;
        this.length = length;
        this.hash = myHashCode(bases, start, length);
    }

    /**
     * Create a new kmer that's a shallow copy of kmer
     * @param kmer the kmer to shallow copy
     */
    public Kmer(final Kmer kmer) {
        this.bases = kmer.bases;
        this.start = kmer.start;
        this.length = kmer.length;
        this.hash = kmer.hash;
    }

    public Kmer(final Kmer kmer, final byte nextChar) {
        final byte[] sequence = new byte[kmer.length];
        System.arraycopy(kmer.bases, kmer.start + 1, sequence, 0, kmer.length - 1);
        sequence[kmer.length - 1] = nextChar;
        bases = sequence;
        start = 0;
        length = kmer.length;
        hash = myHashCode(bases,start,length);
    }

    /**
     * Create a derived shallow kmer that starts at newStart and has newLength bases
     * @param newStart the new start of kmer, where 0 means that start of the kmer, 1 means skip the first base
     * @param newLength the new length
     * @return a new kmer based on the data in this kmer.  Does not make a copy, so shares most of the data
     */
    public Kmer subKmer(final int newStart, final int newLength) {
        return new Kmer(bases, start + newStart, newLength);
    }

    /**
     * Get the bases of this kmer.  May create a copy of the bases, depending on how this kmer was constructed.
     *
     * Note that this function is efficient in that if it needs to copy the bases this only occurs once.
     *
     * @return a non-null byte[] containing length() bases of this kmer, regardless of how this kmer was created
     */
    public byte[] bases() {
        if ( start != 0 || bases.length != length ) {
            // update operation.  Rip out the exact byte[] and update start so we don't ever do this again
            bases = Arrays.copyOfRange(bases, start, start + length);
            start = 0;
        }

        return bases;
    }

    /**
     * Get a string representation of the bases of this kmer
     * @return a non-null string
     */
    public String baseString() {
        return new String(bases());
    }

    /**
     * The length of this kmer
     * @return an integer >= 0
     */
    public int length() {
        return length;
    }

    /**
     * Gets a set of differing positions and bases from another k-mer, limiting up to a max distance.
     * For example, if this = "ACATT" and other = "ACGGT":
     * - if maxDistance < 2 then -1 will be returned, since distance between kmers is 2.
     * - If maxDistance >=2, then 2 will be returned, and arrays will be filled as follows:
     * differingIndeces = {2,3}
     * differingBases = {'G','G'}
     * @param other                 Other k-mer to test
     * @param maxDistance           Maximum distance to search. If this and other k-mers are beyond this Hamming distance,
     *                              search is aborted and a null is returned
     * @param differingIndeces      Array with indices of differing bytes in array
     * @param differingBases        Actual differing bases
     * @return                      Set of mappings of form (int->byte), where each elements represents index
     *                              of k-mer array where bases mismatch, and the byte is the base from other kmer.
     *                              If both k-mers differ by more than maxDistance, returns null
     */
    public int getDifferingPositions(final Kmer other,
                                                   final int maxDistance,
                                                   final int[] differingIndeces,
                                                   final byte[] differingBases) {


        int dist = 0;
        if (length == other.length()) {
            final byte[] f2 = other.bases;
            for (int i=0; i < length; i++) {
                if (bases[start + i] != f2[i]) {
                    differingIndeces[dist] = i;
                    differingBases[dist++] = f2[i];
                    if (dist > maxDistance) {
                        return -1;
                    }
                }
            }

        }
        return dist;
    }

    @Override
    public String toString() {
        return "Kmer{" + new String(bases,start,length) + "}";
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || !Kmer.class.isAssignableFrom(o.getClass())) return false;

        final Kmer kmer = (Kmer) o;

        // very fast test.  If hash aren't equal you are done, otherwise compare the bases
        if ( hash != kmer.hash ) return false;
        if ( length != kmer.length ) return false;

        for ( int i = 0; i < length; i++ ) {
            if (bases[start + i] != kmer.bases[kmer.start + i]) {
                return false;
            }
        }

        return true;
    }

    @Override
    public int hashCode() {
        return hash;
    }

    /**
     * Helper method that computes the hashcode for this kmer based only the bases in
     * a[], starting at start and running length bases
     *
     * @param a a non-null bases array
     * @param start where to start in bases
     * @param length the length of the bases
     * @return a hashcode value appropriate for a[start] -> a[start + length]
     */
    private static int myHashCode(final byte a[], final int start, final int length) {
        if (a == null) {
            return 0;
        }

        int result = 1;
        for (int i = 0; i < length; i++) {
            result = 31 * result + a[start + i];
        }

        return result;
    }

    public byte base(final int i) {
        return bases[start + i];
    }

}
