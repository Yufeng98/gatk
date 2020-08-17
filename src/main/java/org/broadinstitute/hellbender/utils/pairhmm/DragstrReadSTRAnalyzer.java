package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Arrays;

public class DragstrReadSTRAnalyzer {
    private final int[][] repeatsByPeriodAndPosition;
    private final int[] periodWithMostRepeats;
    private final int maxPeriod;
    private int seqLength;

    DragstrReadSTRAnalyzer(final int maxSequenceLength, final int maxPeriod) {
        repeatsByPeriodAndPosition = new int[maxPeriod][maxSequenceLength];
        this.maxPeriod = maxPeriod;
        this.periodWithMostRepeats = new int[maxSequenceLength];
    }

    public int numberOfRepeats(final int position, final int period) {
        if (period <= 0 || period > maxPeriod) {
            return 0;
        } else if (position < 0 || position >= seqLength) {
            throw new IllegalArgumentException("cannot query outside requested boundaries");
        } else {
            return repeatsByPeriodAndPosition[period - 1][position];
        }
    }

    public int mostRepeatedPeriod(final int position) {
        if (position >= 0 && position < seqLength) {
            return periodWithMostRepeats[position];
        } else {
            throw new IllegalArgumentException("cannot query outside requested boundaries");
        }
    }

    public int numberOfMostRepeats(final int position) {
        if (position >= 0 && position < seqLength) {
            return repeatsByPeriodAndPosition[periodWithMostRepeats[position] - 1][position];
        } else {
            throw new IllegalArgumentException("cannot query outside requested boundaries");
        }
    }

    public void load(final byte[] sequence) {
        if (sequence.length > repeatsByPeriodAndPosition[0].length) {
            throw new IllegalArgumentException("input sequence is too long");
        }
        final int seqLength = this.seqLength = sequence.length;
        // periodIndex == periodLength - 1 since in Java indexes start with 0.

        loadPeriodOne(sequence); // simpler and faster code for period-length == 0.
        for (int periodIndex = 1; periodIndex < maxPeriod; periodIndex++) {
            final int periodLength = periodIndex + 1;
            final int[] runLength = repeatsByPeriodAndPosition[periodIndex];
            if (sequence.length < periodLength) {
                Arrays.fill(runLength, 0, sequence.length, 0);
                continue;
            }
            int position, matchedCycles, prevValue, cycleIndex, rightMargin, positionPlusPeriod;

            // we calculate the margins, first right-margin, by looking for
            // the first mismatch after then sequence of interest.

            // Calculate backward the repeat run lengths from a position forward.
            // so that runLength[i] would contain the number of repetitions of seq[i .. i + period) in seq[i .. end_of_seq]
            // First We set period - 1 last positions to 0 as the units will fall off the end of read.
            for (position = seqLength - 1, cycleIndex = 1; cycleIndex < periodLength; position--, cycleIndex++) {
                runLength[position] = 0;
            }

            // Then the main loop.
            // We basically go backward checking what was the number of reported repeats downstream
            // and increasing as we get to a number exact matched cycles that equal to the period (length).
            // When we found a mismatch we reset that count to 0 and the run length to 1 (the minimum run-length).
            prevValue = runLength[position--] = 1; //prevValue holds the num of repeats reported in the previous (+1) position.
            for (positionPlusPeriod = position + periodLength, matchedCycles = 0; position >= 0; position--, positionPlusPeriod--) {
                if (sequence[position] == sequence[positionPlusPeriod]) { // we keep matching repeat unit bases.
                    if (++matchedCycles == periodLength) { // we go a new full repeat matched so the run length increases:
                        prevValue = runLength[position] = prevValue + 1;
                        matchedCycles = 0; // we reset the match-run-length to 0.
                    } else { // we simply copy the run length from the +1 position base.
                        runLength[position] = prevValue;
                    }
                } else { // we bump into a mismatch that end the run.
                    prevValue = runLength[position] = 1;
                    matchedCycles = 0; // we reset the match-run-length to 0.
                }
            }
            // propagate forward the total run-length over to the other repeats in the run.
            // we do it per cycle in the period so that we first deal with repeat units whose offset
            // is zero respect the beginning of the sequence, the2 1 base, then 2 etc.
            for (cycleIndex = 0; cycleIndex < periodLength; cycleIndex++) {
                // The left most repeated unit runLength[i] contains the actual run length for all
                // the units to the right. We copy that value forward to the other run-length units.
                // We do this by iterating over consecutive repeat runs.
                for (position = cycleIndex; position < seqLength; position += periodLength) {
                    final int totalRunLength = runLength[position];
                    for (int repeatInRun = 1; repeatInRun < totalRunLength; repeatInRun++) {
                        runLength[position += periodLength] = totalRunLength;
                    }
                }
            }

            // Now we calculate the max repeat length that overlaps any given position.

            // we skip period == 1 (periodIndex == 0) since is already resolved.

            // for period == 2 the code can be simplified a bit since it only requires the combination of the
            // current value and the previous one

                // for period 3 or above we could use a special heap to get in k log k the longest length in a window
                // but since in practice the max perior is something like 8 or 10 bases it seems a bit of an overkill.
                final int[] periodValues = repeatsByPeriodAndPosition[periodIndex];
                // we fill the expected trailing 0 with the count for the first non-zero which is exactly
                // at length - periodLength (e.q. length - periodIndex - 1).
                int windowEnd, windowStart, maxInWindow;
                windowStart = seqLength - periodLength;
                windowEnd = seqLength;
                maxInWindow = periodValues[windowStart]; // guaranteed to be the the best as other are always 0.
                while (windowEnd > 1) { // not > 0 since the very first value does not need to be changed.
                    int valueOut = periodValues[--windowEnd];
                    int valueIn = --windowStart < 0 ? -1 : periodValues[windowStart];
                    periodValues[windowEnd] = Math.max(maxInWindow, valueIn);
                    if (valueIn < valueOut && valueOut >= maxInWindow) {
                        // here we pay a bit of a penalty to make sure that we get the maxInWindow value correct since
                        // we are removing a value that is the same (or larger) than the current maxInWindow ant the value
                        // coming in is not larger than that.
                        maxInWindow = MathUtils.arrayMax(periodValues, Math.max(0, windowStart), windowEnd, -1);
                    } else if (valueIn > maxInWindow) {
                        // the value coming in is larger than the current max so the maxInWindow update is trivial.
                        maxInWindow = valueIn;
                    }
                }
        }

        // finally we update the periodWithMostRepeats array:

        Arrays.fill(periodWithMostRepeats, 0, seqLength, 1);
        final int[] mostRepeats = repeatsByPeriodAndPosition[0].clone();
        for (int periodIndex = 1; periodIndex < maxPeriod; periodIndex++) {
            final int periodLength = periodIndex + 1;
            final int[] periodValues = repeatsByPeriodAndPosition[periodIndex];
            for (int position = 0; position < seqLength; position++) {
                final int repeats = periodValues[position];
                if (repeats > mostRepeats[position]) {
                    mostRepeats[position] = repeats;
                    periodWithMostRepeats[position] = periodLength;
                }
            }
        }
    }

    private void loadPeriodOne(final byte[] sequence) {
        final int[] runLengths = repeatsByPeriodAndPosition[0];
        final int rightMargin = sequence.length - 1;
        byte last = sequence[rightMargin];
        runLengths[rightMargin] = 1;
        int carryBack = 1;
        for (int position = rightMargin - 1; position >= 0; position--) {
            final byte next = sequence[position];
            runLengths[position] = next == last ? ++carryBack : (carryBack = 1);
            last = next;
        }
        // forward phase:
        last = sequence[0];
        int prevRunLength = 0;
        int carryForward = 0;
        for (int position = 0; position <= rightMargin; position++) {
            final byte next = sequence[position];
            if (next == last) {
                runLengths[position] += carryForward++;
            } else {
                carryForward = 0;
            }
            last = next;
            final int thisRunLength = runLengths[position];
            runLengths[position] = Math.max(thisRunLength, prevRunLength);
            prevRunLength = thisRunLength;
        }
    }
}