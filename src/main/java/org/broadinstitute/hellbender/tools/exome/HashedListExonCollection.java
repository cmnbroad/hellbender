package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Exon collection supported by a list of intervals sorted by location and with quick
 * look-up by name using a hash.
 * <p>
 *     Intervals are sorted using {@link IntervalUtils#LEXICOGRAPHICAL_ORDER_COMPARATOR}.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class HashedListExonCollection<T extends Locatable> implements ExonCollection<T> {

    /**
     * Map from interval name to interval object.
     */
    private final Map<String,T> intervalsByName;

    /**
     * Sorted list of intervals.
     */
    private final List<T> sortedIntervals;

    /**
     * Cached index of the last overlapping interval found by
     * {@link #cachedBinarySearch(SimpleInterval)}.
     *
     * <p>
     *     This results in a noticeable performance gain when processing data in
     *     genomic order.
     * </p>
     */
    private int lastBinarySearchResult = -1;

    /**
     * Creates a exon data-base give a sorted list of intervals.
     *
     * <p>
     *     Intervals will be sorted by their contig id lexicographical order.
     * </p>
     *
     * @param intervals the input interval list. Is assumed to be sorted
     * @throws IllegalArgumentException if {@code intervals} is {@code null}.
     */
    public HashedListExonCollection(final List<T> intervals) {
        if (intervals == null) {
            throw new IllegalArgumentException("the input intervals cannot be null");
        }
        if (intervals.contains(null)) {
            throw new IllegalArgumentException("the input cannot contain null");
        }
        sortedIntervals = intervals.stream().sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList());
        checkForOverlaps(sortedIntervals);
        this.intervalsByName = composeIntervalsByName(sortedIntervals);
    }

    /**
     * Fails with an exception if the intervals collection has overlapping intervals.
     * @param sortedIntervals the intervals sorted.
     */
    private static <T extends Locatable> void checkForOverlaps(final List<T> sortedIntervals) {
        final OptionalInt failureIndex = IntStream.range(1, sortedIntervals.size())
                .filter(i -> IntervalUtils.overlaps(sortedIntervals.get(i-1),sortedIntervals.get(i)))
                .findFirst();

        if (failureIndex.isPresent()) {
            final int index = failureIndex.getAsInt();
            throw new IllegalArgumentException(
                    String.format("input intervals contain at least two overlapping intervals: %s and %s",
                            sortedIntervals.get(index-1),sortedIntervals.get(index)));
        }
    }

    /**
     * Composes a map from name to interval.
     *
     * <p>
     *     Also check whether there are more than one interval with the same name and in that case
     *     it fails with an exception.
     * </p>
     *
     * @param intervals the input intervals.
     * @throws IllegalArgumentException if {@code intervals} is {@code null} or it contains {@code nulls},
     *                  or is not a coherent interval list as per the criteria above.
     */
    private Map<String,T> composeIntervalsByName(final List<T> intervals) {
        final Map<String,T> result = new HashMap<>(intervals.size());
        for (final T location : intervals) {
            final String name = name(location);
            if (name != null) {
                final T previous = result.put(name,location);
                if (previous != null) {
                    throw new IllegalStateException(
                            String.format("more than one interval in the input list results in the same name (%s); " +
                                            "perhaps repeated: '%s' and '%s'.",
                                    name,previous,location));
                }
            }
        }
        return result;
    }

    @Override
    public String name(final T exon) {
        Utils.nonNull(exon,"null exon not allowed");
        final String contig = exon.getContig();
        if (contig == null) {
            return null;
        } else {
            return String.format("%s:%d-%d", contig, exon.getStart(), exon.getEnd());
        }
    }

    @Override
    public List<T> exons() {
        return sortedIntervals;
    }

    @Override
    public int exonCount() {
        return sortedIntervals.size();
    }

    @Override
    public T exon(final int index) {
        Utils.validIndex(index, sortedIntervals.size());
        return sortedIntervals.get(index);
    }

    @Override
    public T exon(final String name) {
        Utils.nonNull(name,"the input name cannot be null");
        return intervalsByName.get(name);
    }

    @Override
    public int index(final String name) {
        final T exon = intervalsByName.get(name);
        if (exon == null) {
            return -1;
        } else {
            final int searchIndex = uncachedBinarySearch(location(exon));
            if (searchIndex < 0) { // checking just in case.
                throw new IllegalStateException("could not found named interval amongst sorted intervals, impossible");
            }
            return searchIndex;
        }
    }

    @Override
    public IndexRange indexRange(final SimpleInterval location) {
        Utils.nonNull(location, "the input location cannot be null");
        final int searchIndex = cachedBinarySearch(location);
        if (searchIndex < 0) {
            return new IndexRange(-searchIndex - 1, -searchIndex - 1);
        } else {
            final int firstOverlappingIndex = extendSearchIndexBackwards(location, searchIndex);
            final int lastOverlappingIndex = extendSearchIndexForward(location, searchIndex);
            return new IndexRange(firstOverlappingIndex, lastOverlappingIndex + 1);
        }
    }

    @Override
    public SimpleInterval location(final T exon) {
        return new SimpleInterval(Utils.nonNull(exon,"the exon cannot be null"));
    }

    @Override
    public SimpleInterval location(final int index) {
        Utils.validIndex(index,sortedIntervals.size());
        return new SimpleInterval(sortedIntervals.get(index));
    }

    @Override
    public T exon(final SimpleInterval overlapRegion) {
        final int searchIndex = index(overlapRegion);
        return searchIndex < 0 ? null : sortedIntervals.get(searchIndex);
    }

    @Override
    public int index(final SimpleInterval location) {
        final IndexRange range = indexRange(location);
        switch (range.size()) {
            case 1:
                return range.from;
            case 0:
                return - (range.from + 1);
            default:
                throw new AmbiguousExonException(
                        String.format("location '%s' overlaps with %d exons: from '%s' to '%s'.",
                  location,range.size(),exon(range.from),exon(range.to - 1)));
        }
    }

    @Override
    public List<T> exons(final SimpleInterval overlapRegion) {
        final IndexRange range = indexRange(overlapRegion);
        return sortedIntervals.subList(range.from, range.to);
    }

    /**
     * Implements a cached binary search of the overlapping intervals.
     *
     * <p>
     *     This was found to improve performance significantly when analyzing empirical
     *     exome data in sequential order.
     * </p>
     *
     * <p>
     *     First it checks
     *     whether the query interval overlaps with the result of the
     *     last search (whether a hit or a miss (taking the insertion position of the miss)).
     * </p>
     * <p>
     *     If the last search was a hit but the new query location does not overlap,
     *     it checks on the next interval and if not fails over to the regular binary search.
     * </p>
     *
     * <p>
     *     A positive (0 or greater) returned value indicates a hit where the corresponding interval
     *     is guaranteed to overlap the query {@code location}. There might be more intervals that
     *     overlap this location and they must all be contiguous to the returned index.
     * <p>
     *     One can use {@link #extendSearchIndexBackwards(SimpleInterval, int)}
     *     and {@link #extendSearchIndexForward(SimpleInterval, int)}
     *     to find the actual overlapping index range.
     * </p>
     *
     * <p>
     *     In contrast, a negative result indicates that there is no overlapping interval.
     *     This value encodes the insertion position for the query location; that is, where would
     *     the query location be inserted were it to be added to the sorted list of intervals:
     * </p>
     * <p>
     *     <code>insertion index == -(result + 1)</code>
     * </p>
     *
     * @param location the query location.
     * @return any integer between <code>-{@link #exonCount()}-1</code> and <code>{@link #exonCount()} - 1</code>.
     */
    private int cachedBinarySearch(final SimpleInterval location) {

        if (lastBinarySearchResult < 0) {
            final int candidate = -(lastBinarySearchResult + 1);
            if (candidate >= sortedIntervals.size()) {
                return lastBinarySearchResult = uncachedBinarySearch(location);
            } else if (IntervalUtils.overlaps(sortedIntervals.get(candidate),location)) {
                return lastBinarySearchResult = candidate;
            } else {
                return lastBinarySearchResult = uncachedBinarySearch(location);
            }
        } else {
            if (IntervalUtils.overlaps(sortedIntervals.get(lastBinarySearchResult),location)) {
                return lastBinarySearchResult;
            } else {
                final int candidate = lastBinarySearchResult + 1;
                if (candidate == sortedIntervals.size()) {
                    return lastBinarySearchResult = uncachedBinarySearch(location);
                } else if (IntervalUtils.overlaps(sortedIntervals.get(candidate),location)) {
                    return lastBinarySearchResult = candidate;
                } else {
                    return lastBinarySearchResult = uncachedBinarySearch(location);
                }
            }
        }
    }

    /**
     * Implements a binary search of the overlapping intervals.
     *
     * <p>
     *     A positive (0 or greater) returned value indicates a hit where the corresponding interval
     *     is guaranteed to overlap the query {@code location}. There might be more intervals that
     *     overlap this location and they must all be contiguous to the returned index.
     * <p>
     *     One can use {@link #extendSearchIndexBackwards(SimpleInterval, int)}
     *     and {@link #extendSearchIndexForward(SimpleInterval, int)}
     *     to find the actual overlapping index range.
     * </p>
     *
     * <p>
     *     In contrast, a negative result indicates that there is no overlapping interval.
     *     This value encodes the insertion position for the query location; that is, where would
     *     the query location be inserted were it to be added to the sorted list of intervals:
     * </p>
     * <p>
     *     <code>insertion index == -(result + 1)</code>
     * </p>
     *
     * @param location the query location.
     * @return any integer between <code>-{@link #exonCount()}-1</code> and <code>{@link #exonCount()} - 1</code>.
     */
    private int uncachedBinarySearch(final SimpleInterval location) {
        if (sortedIntervals.size() == 0) {
            return -1;
        }

        final int searchResult = Collections.binarySearch(sortedIntervals, location, IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        if (searchResult >= 0) {
            return searchResult;
        } else {
            final int insertIndex = - (searchResult + 1);
            if (insertIndex < sortedIntervals.size() && IntervalUtils.overlaps(sortedIntervals.get(insertIndex),location)) {
                return insertIndex;
            } if (insertIndex > 0 && IntervalUtils.overlaps(sortedIntervals.get(insertIndex - 1),location)) {
                return insertIndex - 1;
            } else {
                return searchResult;
            }
        }
    }

    /**
     * Looks for the last index in {@link #sortedIntervals} that has an overlap with the input {@code location}.
     * starting at {@code startIndex} and assuming that the element at that index has an
     * overlap with {@code location}.
     */
    private int extendSearchIndexForward(final SimpleInterval location, final int startIndex) {
        final ListIterator<T> it = sortedIntervals.listIterator(startIndex + 1);
        while (it.hasNext()) {
            final T next = it.next();
            if (!IntervalUtils.overlaps(location,next)) {
                return it.previousIndex() - 1;
            }
        }
        return it.previousIndex();
    }

    /**
     * Looks for the first index in {@link #sortedIntervals} that has an overlap with the input {@code location}
     * starting at {@code startIndex} and assuming that the element at that index has an overlap with {@code location}.
     */
    private int extendSearchIndexBackwards(final SimpleInterval location, final int startIndex) {
        final ListIterator<T> it = sortedIntervals.listIterator(startIndex);
        while (it.hasPrevious()) {
            final T previous = it.previous();
            if (!IntervalUtils.overlaps(location,previous)) {
                return it.nextIndex() + 1;
            }
        }
        return it.nextIndex();
    }

}
