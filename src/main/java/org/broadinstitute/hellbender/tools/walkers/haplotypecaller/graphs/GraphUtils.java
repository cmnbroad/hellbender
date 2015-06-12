package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.collections.PrimitivePair;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Utility functions used in the graphs package
 */
public final class GraphUtils {
    private GraphUtils() {}

    /**
     * Compute the maximum shared prefix length of list of bytes.
     *
     * @param listOfBytes a list of bytes with at least one element
     * @param minLength the min. length among all byte[] in listOfBytes
     * @return the number of shared bytes common at the start of all bytes
     */
    public static int compPrefixLen(final List<byte[]> listOfBytes, final int minLength) {
        for ( int i = 0; i < minLength; i++ ) {
            final byte b = listOfBytes.get(0)[i];
            for ( int j = 1; j < listOfBytes.size(); j++ ) {
                if ( b != listOfBytes.get(j)[i] ) {
                    return i;
                }
            }
        }

        return minLength;
    }

    /**
     * Compute the maximum shared suffix length of list of bytes.
     *
     * @param listOfBytes a list of bytes with at least one element
     * @param minLength the min. length among all byte[] in listOfBytes
     * @return the number of shared bytes common at the end of all bytes
     */
    public static int compSuffixLen(final List<byte[]> listOfBytes, final int minLength) {
        for ( int suffixLen = 0; suffixLen < minLength; suffixLen++ ) {
            final byte b = listOfBytes.get(0)[listOfBytes.get(0).length - suffixLen - 1];
            for ( int j = 1; j < listOfBytes.size(); j++ ) {
                if ( b != listOfBytes.get(j)[listOfBytes.get(j).length - suffixLen - 1] )
                    return suffixLen;
            }
        }
        return minLength;
    }

    /**
     * Get the list of kmers as byte[] from the vertices in the graph
     *
     * @param vertices a collection of vertices
     * @return a list of their kmers in order of the iterator on vertices
     */
    public static List<byte[]> getKmers(final Collection<SeqVertex> vertices) {
        final List<byte[]> kmers = new ArrayList<>(vertices.size());
        for ( final SeqVertex v : vertices ) {
            kmers.add(v.getSequence());
        }
        return kmers;
    }

    /**
     * Get the minimum length of a collection of byte[]
     *
     * @param kmers a list of kmers whose .length min we want
     * @return the min of the kmers, if kmers is empty the result is 0
     */
    public static int minKmerLength(final Collection<byte[]> kmers) {
        if ( kmers == null ) throw new IllegalArgumentException("kmers cannot be null");

        if ( kmers.isEmpty() ) return 0;
        int min = Integer.MAX_VALUE;
        for ( final byte[] kmer : kmers ) {
            min = Math.min(min, kmer.length);
        }
        return min;
    }

    /**
     * Find the ending position of the longest uniquely matching
     * run of bases of kmer in seq.
     *
     * for example, if seq = ACGT and kmer is NAC, this function returns 1,2 as we have the following
     * match:
     *
     *  0123
     * .ACGT
     * NAC..
     *
     * @param seq a non-null sequence of bytes
     * @param kmer a non-null kmer
     * @return the ending position and length where kmer matches uniquely in sequence, or null if no
     *         unique longest match can be found
     */
    public static PrimitivePair.Int findLongestUniqueSuffixMatch(final byte[] seq, final byte[] kmer) {
        int longestPos = -1;
        int length = 0;
        boolean foundDup = false;

        for ( int i = 0; i < seq.length; i++ ) {
            final int matchSize = longestSuffixMatch(seq, kmer, i);
            if ( matchSize > length ) {
                longestPos = i;
                length = matchSize;
                foundDup = false;
            } else if ( matchSize == length ) {
                foundDup = true;
            }
        }

        return foundDup ? null : new PrimitivePair.Int(longestPos, length);
    }

    /**
     * calculates the longest suffix match between a sequence and a smaller kmer
     *
     * @param seq         the (reference) sequence
     * @param kmer        the smaller kmer sequence
     * @param seqStart    the index (inclusive) on seq to start looking backwards from
     * @return the longest matching suffix
     */
    public static int longestSuffixMatch(final byte[] seq, final byte[] kmer, final int seqStart) {
        for ( int len = 1; len <= kmer.length; len++ ) {
            final int seqI = seqStart - len + 1;
            final int kmerI = kmer.length - len;
            if ( seqI < 0 || seq[seqI] != kmer[kmerI] ) {
                return len - 1;
            }
        }
        return kmer.length;
    }
}
