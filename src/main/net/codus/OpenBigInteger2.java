package net.codus;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.Random;

/**
 * Immutable arbitrary-precision integers.  All operations behave as if
 * OpenBigIntegers were represented in two's-complement notation (like Java's
 * primitive integer types).  OpenBigInteger provides analogues to all of Java's
 * primitive integer operators, and all relevant methods from java.lang.Math.
 * Additionally, OpenBigInteger provides operations for modular arithmetic, GCD
 * calculation, primality testing, prime generation, bit manipulation,
 * and a few other miscellaneous operations.
 *
 * <p>Semantics of arithmetic operations exactly mimic those of Java's integer
 * arithmetic operators, as defined in <i>The Java Language Specification</i>.
 * For example, division by zero throws an {@code ArithmeticException}, and
 * division of a negative by a positive yields a negative (or zero) remainder.
 * All of the details in the Spec concerning overflow are ignored, as
 * OpenBigIntegers are made as large as necessary to accommodate the results of an
 * operation.
 *
 * <p>Semantics of shift operations extend those of Java's shift operators
 * to allow for negative shift distances.  A right-shift with a negative
 * shift distance results in a left shift, and vice-versa.  The unsigned
 * right shift operator ({@code >>>}) is omitted, as this operation makes
 * little sense in combination with the "infinite word size" abstraction
 * provided by this class.
 *
 * <p>Semantics of bitwise logical operations exactly mimic those of Java's
 * bitwise integer operators.  The binary operators ({@code and},
 * {@code or}, {@code xor}) implicitly perform sign extension on the shorter
 * of the two operands prior to performing the operation.
 *
 * <p>Comparison operations perform signed integer comparisons, analogous to
 * those performed by Java's relational and equality operators.
 *
 * <p>Modular arithmetic operations are provided to compute residues, perform
 * exponentiation, and compute multiplicative inverses.  These methods always
 * return a non-negative result, between {@code 0} and {@code (modulus - 1)},
 * inclusive.
 *
 * <p>Bit operations operate on a single bit of the two's-complement
 * representation of their operand.  If necessary, the operand is sign-
 * extended so that it contains the designated bit.  None of the single-bit
 * operations can produce a OpenBigInteger with a different sign from the
 * OpenBigInteger being operated on, as they affect only a single bit, and the
 * "infinite word size" abstraction provided by this class ensures that there
 * are infinitely many "virtual sign bits" preceding each OpenBigInteger.
 *
 * <p>For the sake of brevity and clarity, pseudo-code is used throughout the
 * descriptions of OpenBigInteger methods.  The pseudo-code expression
 * {@code (i + j)} is shorthand for "a OpenBigInteger whose value is
 * that of the OpenBigInteger {@code i} plus that of the OpenBigInteger {@code j}."
 * The pseudo-code expression {@code (i == j)} is shorthand for
 * "{@code true} if and only if the OpenBigInteger {@code i} represents the same
 * value as the OpenBigInteger {@code j}."  Other pseudo-code expressions are
 * interpreted similarly.
 *
 * <p>All methods and constructors in this class throw
 * {@code NullPointerException} when passed
 * a null object reference for any input parameter.
 *
 * @see     BigDecimal
 * @author  Josh Bloch
 * @author  Michael McCloskey
 * @since JDK1.1
 */

public class OpenBigInteger2 extends Number implements Comparable<OpenBigInteger2> {
    /**
     * The signum of this OpenBigInteger: -1 for negative, 0 for zero, or
     * 1 for positive.  Note that the OpenBigInteger zero <i>must</i> have
     * a signum of 0.  This is necessary to ensures that there is exactly one
     * representation for each OpenBigInteger value.
     *
     * @serial
     */
    int signum;

    /**
     * The magnitude of this OpenBigInteger, in <i>big-endian</i> order: the
     * zeroth element of this array is the most-significant int of the
     * magnitude.  The magnitude must be "minimal" in that the most-significant
     * int ({@code mag[0]}) must be non-zero.  This is necessary to
     * ensure that there is exactly one representation for each OpenBigInteger
     * value.  Note that this implies that the OpenBigInteger zero has a
     * zero-length mag array.
     */
    int[] mag;

    // These "redundant fields" are initialized with recognizable nonsense
    // values, and cached the first time they are needed (or never, if they
    // aren't needed).

     /**
     * One plus the bitCount of this OpenBigInteger. Zeros means unitialized.
     *
     * @serial
     * @see #bitCount
     * @deprecated Deprecated since logical value is offset from stored
     * value and correction factor is applied in accessor method.
     */
    @Deprecated
    private int bitCount;

    /**
     * One plus the bitLength of this OpenBigInteger. Zeros means unitialized.
     * (either value is acceptable).
     *
     * @serial
     * @see #bitLength()
     * @deprecated Deprecated since logical value is offset from stored
     * value and correction factor is applied in accessor method.
     */
    @Deprecated
    private int bitLength;

    /**
     * Two plus the lowest set bit of this OpenBigInteger, as returned by
     * getLowestSetBit().
     *
     * @serial
     * @see #getLowestSetBit
     * @deprecated Deprecated since logical value is offset from stored
     * value and correction factor is applied in accessor method.
     */
    @Deprecated
    private int lowestSetBit;

    /**
     * Two plus the index of the lowest-order int in the magnitude of this
     * OpenBigInteger that contains a nonzero int, or -2 (either value is acceptable).
     * The least significant int has int-number 0, the next int in order of
     * increasing significance has int-number 1, and so forth.
     * @deprecated Deprecated since logical value is offset from stored
     * value and correction factor is applied in accessor method.
     */
    @Deprecated
    private int firstNonzeroIntNum;

    /**
     * This mask is used to obtain the value of an int as if it were unsigned.
     */
    final static long LONG_MASK = 0xffffffffL;

    //Constructors

    /**
     * Translates a byte array containing the two's-complement binary
     * representation of a OpenBigInteger into a OpenBigInteger.  The input array is
     * assumed to be in <i>big-endian</i> byte-order: the most significant
     * byte is in the zeroth element.
     *
     * @param  val big-endian two's-complement binary representation of
     *         OpenBigInteger.
     * @throws NumberFormatException {@code val} is zero bytes long.
     */
    public OpenBigInteger2(byte[] val) {
        if (val.length == 0)
            throw new NumberFormatException("Zero length OpenBigInteger");

        if (val[0] < 0) {
            mag = makePositive(val);
            signum = -1;
        } else {
            mag = stripLeadingZeroBytes(val);
            signum = (mag.length == 0 ? 0 : 1);
        }
    }

    /**
     * This private constructor translates an int array containing the
     * two's-complement binary representation of a OpenBigInteger into a
     * OpenBigInteger. The input array is assumed to be in <i>big-endian</i>
     * int-order: the most significant int is in the zeroth element.
     */
    private OpenBigInteger2(int[] val) {
        if (val.length == 0)
            throw new NumberFormatException("Zero length OpenBigInteger");

        if (val[0] < 0) {
            mag = makePositive(val);
            signum = -1;
        } else {
            mag = trustedStripLeadingZeroInts(val);
            signum = (mag.length == 0 ? 0 : 1);
        }
    }

    /**
     * Translates the sign-magnitude representation of a OpenBigInteger into a
     * OpenBigInteger.  The sign is represented as an integer signum value: -1 for
     * negative, 0 for zero, or 1 for positive.  The magnitude is a byte array
     * in <i>big-endian</i> byte-order: the most significant byte is in the
     * zeroth element.  A zero-length magnitude array is permissible, and will
     * result in a OpenBigInteger value of 0, whether signum is -1, 0 or 1.
     *
     * @param  signum signum of the number (-1 for negative, 0 for zero, 1
     *         for positive).
     * @param  magnitude big-endian binary representation of the magnitude of
     *         the number.
     * @throws NumberFormatException {@code signum} is not one of the three
     *         legal values (-1, 0, and 1), or {@code signum} is 0 and
     *         {@code magnitude} contains one or more non-zero bytes.
     */
    public OpenBigInteger2(int signum, byte[] magnitude) {
        this.mag = stripLeadingZeroBytes(magnitude);

        if (signum < -1 || signum > 1)
            throw(new NumberFormatException("Invalid signum value"));

        if (this.mag.length==0) {
            this.signum = 0;
        } else {
            if (signum == 0)
                throw(new NumberFormatException("signum-magnitude mismatch"));
            this.signum = signum;
        }
    }

    /**
     * A constructor for internal use that translates the sign-magnitude
     * representation of a OpenBigInteger into a OpenBigInteger. It checks the
     * arguments and copies the magnitude so this constructor would be
     * safe for external use.
     */
    private OpenBigInteger2(int signum, int[] magnitude) {
        this.mag = stripLeadingZeroInts(magnitude);

        if (signum < -1 || signum > 1)
            throw(new NumberFormatException("Invalid signum value"));

        if (this.mag.length==0) {
            this.signum = 0;
        } else {
            if (signum == 0)
                throw(new NumberFormatException("signum-magnitude mismatch"));
            this.signum = signum;
        }
    }

    /**
     * Translates the String representation of a OpenBigInteger in the
     * specified radix into a OpenBigInteger.  The String representation
     * consists of an optional minus or plus sign followed by a
     * sequence of one or more digits in the specified radix.  The
     * character-to-digit mapping is provided by {@code
     * Character.digit}.  The String may not contain any extraneous
     * characters (whitespace, for example).
     *
     * @param val String representation of OpenBigInteger.
     * @param radix radix to be used in interpreting {@code val}.
     * @throws NumberFormatException {@code val} is not a valid representation
     *         of a OpenBigInteger in the specified radix, or {@code radix} is
     *         outside the range from {@link Character#MIN_RADIX} to
     *         {@link Character#MAX_RADIX}, inclusive.
     * @see    Character#digit
     */
    public OpenBigInteger2(String val, int radix) {
        int cursor = 0, numDigits;
        final int len = val.length();

        if (radix < Character.MIN_RADIX || radix > Character.MAX_RADIX)
            throw new NumberFormatException("Radix out of range");
        if (len == 0)
            throw new NumberFormatException("Zero length OpenBigInteger");

        // Check for at most one leading sign
        int sign = 1;
        int index1 = val.lastIndexOf('-');
        int index2 = val.lastIndexOf('+');
        if ((index1 + index2) <= -1) {
            // No leading sign character or at most one leading sign character
            if (index1 == 0 || index2 == 0) {
                cursor = 1;
                if (len == 1)
                    throw new NumberFormatException("Zero length OpenBigInteger");
            }
            if (index1 == 0)
                sign = -1;
        } else
            throw new NumberFormatException("Illegal embedded sign character");

        // Skip leading zeros and compute number of digits in magnitude
        while (cursor < len &&
               Character.digit(val.charAt(cursor), radix) == 0)
            cursor++;
        if (cursor == len) {
            signum = 0;
            mag = ZERO.mag;
            return;
        }

        numDigits = len - cursor;
        signum = sign;

        // Pre-allocate array of expected size. May be too large but can
        // never be too small. Typically exact.
        int numBits = (int)(((numDigits * bitsPerDigit[radix]) >>> 10) + 1);
        int numWords = (numBits + 31) >>> 5;
        int[] magnitude = new int[numWords];

        // Process first (potentially short) digit group
        int firstGroupLen = numDigits % digitsPerInt[radix];
        if (firstGroupLen == 0)
            firstGroupLen = digitsPerInt[radix];
        String group = val.substring(cursor, cursor += firstGroupLen);
        magnitude[numWords - 1] = Integer.parseInt(group, radix);
        if (magnitude[numWords - 1] < 0)
            throw new NumberFormatException("Illegal digit");

        // Process remaining digit groups
        int superRadix = intRadix[radix];
        int groupVal = 0;
        while (cursor < len) {
            group = val.substring(cursor, cursor += digitsPerInt[radix]);
            groupVal = Integer.parseInt(group, radix);
            if (groupVal < 0)
                throw new NumberFormatException("Illegal digit");
            destructiveMulAdd(magnitude, superRadix, groupVal);
        }
        // Required for cases where the array was overallocated.
        mag = trustedStripLeadingZeroInts(magnitude);
    }

    // Constructs a new OpenBigInteger using a char array with radix=10
    OpenBigInteger2(char[] val) {
        int cursor = 0, numDigits;
        int len = val.length;

        // Check for leading minus sign
        int sign = 1;
        if (val[0] == '-') {
            if (len == 1)
                throw new NumberFormatException("Zero length OpenBigInteger");
            sign = -1;
            cursor = 1;
        } else if (val[0] == '+') {
            if (len == 1)
                throw new NumberFormatException("Zero length OpenBigInteger");
            cursor = 1;
        }

        // Skip leading zeros and compute number of digits in magnitude
        while (cursor < len && Character.digit(val[cursor], 10) == 0)
            cursor++;
        if (cursor == len) {
            signum = 0;
            mag = ZERO.mag;
            return;
        }

        numDigits = len - cursor;
        signum = sign;

        // Pre-allocate array of expected size
        int numWords;
        if (len < 10) {
            numWords = 1;
        } else {
            int numBits = (int)(((numDigits * bitsPerDigit[10]) >>> 10) + 1);
            numWords = (numBits + 31) >>> 5;
        }
        int[] magnitude = new int[numWords];

        // Process first (potentially short) digit group
        int firstGroupLen = numDigits % digitsPerInt[10];
        if (firstGroupLen == 0)
            firstGroupLen = digitsPerInt[10];
        magnitude[numWords - 1] = parseInt(val, cursor,  cursor += firstGroupLen);

        // Process remaining digit groups
        while (cursor < len) {
            int groupVal = parseInt(val, cursor, cursor += digitsPerInt[10]);
            destructiveMulAdd(magnitude, intRadix[10], groupVal);
        }
        mag = trustedStripLeadingZeroInts(magnitude);
    }

    // Create an integer with the digits between the two indexes
    // Assumes start < end. The result may be negative, but it
    // is to be treated as an unsigned value.
    private int parseInt(char[] source, int start, int end) {
        int result = Character.digit(source[start++], 10);
        if (result == -1)
            throw new NumberFormatException(new String(source));

        for (int index = start; index<end; index++) {
            int nextVal = Character.digit(source[index], 10);
            if (nextVal == -1)
                throw new NumberFormatException(new String(source));
            result = 10*result + nextVal;
        }

        return result;
    }

    // bitsPerDigit in the given radix times 1024
    // Rounded up to avoid underallocation.
    private static long bitsPerDigit[] = { 0, 0,
        1024, 1624, 2048, 2378, 2648, 2875, 3072, 3247, 3402, 3543, 3672,
        3790, 3899, 4001, 4096, 4186, 4271, 4350, 4426, 4498, 4567, 4633,
        4696, 4756, 4814, 4870, 4923, 4975, 5025, 5074, 5120, 5166, 5210,
                                           5253, 5295};

    // Multiply x array times word y in place, and add word z
    private static void destructiveMulAdd(int[] x, int y, int z) {
        // Perform the multiplication word by word
        long ylong = y & LONG_MASK;
        long zlong = z & LONG_MASK;
        int len = x.length;

        long product = 0;
        long carry = 0;
        for (int i = len-1; i >= 0; i--) {
            product = ylong * (x[i] & LONG_MASK) + carry;
            x[i] = (int)product;
            carry = product >>> 32;
        }

        // Perform the addition
        long sum = (x[len-1] & LONG_MASK) + zlong;
        x[len-1] = (int)sum;
        carry = sum >>> 32;
        for (int i = len-2; i >= 0; i--) {
            sum = (x[i] & LONG_MASK) + carry;
            x[i] = (int)sum;
            carry = sum >>> 32;
        }
    }

    /**
     * Translates the decimal String representation of a OpenBigInteger into a
     * OpenBigInteger.  The String representation consists of an optional minus
     * sign followed by a sequence of one or more decimal digits.  The
     * character-to-digit mapping is provided by {@code Character.digit}.
     * The String may not contain any extraneous characters (whitespace, for
     * example).
     *
     * @param val decimal String representation of OpenBigInteger.
     * @throws NumberFormatException {@code val} is not a valid representation
     *         of a OpenBigInteger.
     * @see    Character#digit
     */
    public OpenBigInteger2(String val) {
        this(val, 10);
    }

    /**
     * Constructs a randomly generated OpenBigInteger, uniformly distributed over
     * the range 0 to (2<sup>{@code numBits}</sup> - 1), inclusive.
     * The uniformity of the distribution assumes that a fair source of random
     * bits is provided in {@code rnd}.  Note that this constructor always
     * constructs a non-negative OpenBigInteger.
     *
     * @param  numBits maximum bitLength of the new OpenBigInteger.
     * @param  rnd source of randomness to be used in computing the new
     *         OpenBigInteger.
     * @throws IllegalArgumentException {@code numBits} is negative.
     * @see #bitLength()
     */
    public OpenBigInteger2(int numBits, Random rnd) {
        this(1, randomBits(numBits, rnd));
    }

    private static byte[] randomBits(int numBits, Random rnd) {
        if (numBits < 0)
            throw new IllegalArgumentException("numBits must be non-negative");
        int numBytes = (int)(((long)numBits+7)/8); // avoid overflow
        byte[] randomBits = new byte[numBytes];

        // Generate random bytes and mask out any excess bits
        if (numBytes > 0) {
            rnd.nextBytes(randomBits);
            int excessBits = 8*numBytes - numBits;
            randomBits[0] &= (1 << (8-excessBits)) - 1;
        }
        return randomBits;
    }

    // Minimum size in bits that the requested prime number has
    // before we use the large prime number generating algorithms
    private static final int SMALL_PRIME_THRESHOLD = 95;

    // Certainty required to meet the spec of probablePrime
    private static final int DEFAULT_PRIME_CERTAINTY = 100;

    private static final OpenBigInteger2 SMALL_PRIME_PRODUCT
                       = valueOf(3L*5*7*11*13*17*19*23*29*31*37*41);

    private static volatile Random staticRandom;

    private static Random getSecureRandom() {
        if (staticRandom == null) {
            staticRandom = new java.security.SecureRandom();
        }
        return staticRandom;
    }

    /**
     * This internal constructor differs from its public cousin
     * with the arguments reversed in two ways: it assumes that its
     * arguments are correct, and it doesn't copy the magnitude array.
     */
    OpenBigInteger2(int[] magnitude, int signum) {
        this.signum = (magnitude.length==0 ? 0 : signum);
        this.mag = magnitude;
    }

    /**
     * This private constructor is for internal use and assumes that its
     * arguments are correct.
     */
    private OpenBigInteger2(byte[] magnitude, int signum) {
        this.signum = (magnitude.length==0 ? 0 : signum);
        this.mag = stripLeadingZeroBytes(magnitude);
    }

    //Static Factory Methods

    /**
     * Returns a OpenBigInteger whose value is equal to that of the
     * specified {@code long}.  This "static factory method" is
     * provided in preference to a ({@code long}) constructor
     * because it allows for reuse of frequently used OpenBigIntegers.
     *
     * @param  val value of the OpenBigInteger to return.
     * @return a OpenBigInteger with the specified value.
     */
    public static OpenBigInteger2 valueOf(long val) {
        // If -MAX_CONSTANT < val < MAX_CONSTANT, return stashed constant
        if (val == 0)
            return ZERO;
        if (val > 0 && val <= MAX_CONSTANT)
            return posConst[(int) val];
        else if (val < 0 && val >= -MAX_CONSTANT)
            return negConst[(int) -val];

        return new OpenBigInteger2(val);
    }

    /**
     * Constructs a OpenBigInteger with the specified value, which may not be zero.
     */
    private OpenBigInteger2(long val) {
        if (val < 0) {
            val = -val;
            signum = -1;
        } else {
            signum = 1;
        }

        int highWord = (int)(val >>> 32);
        if (highWord==0) {
            mag = new int[1];
            mag[0] = (int)val;
        } else {
            mag = new int[2];
            mag[0] = highWord;
            mag[1] = (int)val;
        }
    }

    /**
     * Returns a OpenBigInteger with the given two's complement representation.
     * Assumes that the input array will not be modified (the returned
     * OpenBigInteger will reference the input array if feasible).
     */
    private static OpenBigInteger2 valueOf(int val[]) {
        return (val[0]>0 ? new OpenBigInteger2(val, 1) : new OpenBigInteger2(val));
    }

    // Constants

    /**
     * Initialize static constant array when class is loaded.
     */
    private final static int MAX_CONSTANT = 16;
    private static OpenBigInteger2 posConst[] = new OpenBigInteger2[MAX_CONSTANT+1];
    private static OpenBigInteger2 negConst[] = new OpenBigInteger2[MAX_CONSTANT+1];
    static {
        for (int i = 1; i <= MAX_CONSTANT; i++) {
            int[] magnitude = new int[1];
            magnitude[0] = i;
            posConst[i] = new OpenBigInteger2(magnitude,  1);
            negConst[i] = new OpenBigInteger2(magnitude, -1);
        }
    }

    /**
     * The OpenBigInteger constant zero.
     *
     * @since   1.2
     */
    public static final OpenBigInteger2 ZERO = new OpenBigInteger2(new int[0], 0);

    /**
     * The OpenBigInteger constant one.
     *
     * @since   1.2
     */
    public static final OpenBigInteger2 ONE = valueOf(1);

    /**
     * The OpenBigInteger constant two.  (Not exported.)
     */
    private static final OpenBigInteger2 TWO = valueOf(2);

    /**
     * The OpenBigInteger constant ten.
     *
     * @since   1.5
     */
    public static final OpenBigInteger2 TEN = valueOf(10);

    // Arithmetic Operations

    /**
     * Returns a OpenBigInteger whose value is {@code (this + val)}.
     *
     * @param  val value to be added to this OpenBigInteger.
     * @return {@code this + val}
     */
    public OpenBigInteger2 add(OpenBigInteger2 val) {
        if (val.signum == 0)
            return this;
        if (signum == 0)
            return val;
        if (val.signum == signum)
            return new OpenBigInteger2(add(mag, val.mag), signum);

        int cmp = compareMagnitude(val);
        if (cmp == 0)
            return ZERO;
        int[] resultMag = (cmp > 0 ? subtract(mag, start, val.mag,val.start)
                           : subtract(val.mag, val.start, mag, start));
        resultMag = trustedStripLeadingZeroInts(resultMag);

        return new OpenBigInteger2(resultMag, cmp == signum ? 1 : -1);
    }

    /**
     * Adds the contents of the int arrays x and y. This method allocates
     * a new int array to hold the answer and returns a reference to that
     * array.
     */
    private static int[] add(int[] x, int[] y) {
        // If x is shorter, swap the two arrays
        if (x.length < y.length) {
            int[] tmp = x;
            x = y;
            y = tmp;
        }

        int xIndex = x.length;
        int yIndex = y.length;
        int result[] = new int[xIndex];
        long sum = 0;

        // Add common parts of both numbers
        while(yIndex > 0) {
            sum = (x[--xIndex] & LONG_MASK) +
                  (y[--yIndex] & LONG_MASK) + (sum >>> 32);
            result[xIndex] = (int)sum;
        }

        // Copy remainder of longer number while carry propagation is required
        boolean carry = (sum >>> 32 != 0);
        while (xIndex > 0 && carry)
            carry = ((result[--xIndex] = x[xIndex] + 1) == 0);

        // Copy remainder of longer number
        while (xIndex > 0)
            result[--xIndex] = x[xIndex];

        // Grow result if necessary
        if (carry) {
            int bigger[] = new int[result.length + 1];
            System.arraycopy(result, 0, bigger, 1, result.length);
            bigger[0] = 0x01;
            return bigger;
        }
        return result;
    }

    /**
     * Returns a OpenBigInteger whose value is {@code (this - val)}.
     *
     * @param  val value to be subtracted from this OpenBigInteger.
     * @return {@code this - val}
     */
    public OpenBigInteger2 subtract(OpenBigInteger2 val) {
        if (val.signum == 0)
            return this;
        if (signum == 0)
            return val.negate();
        if (val.signum != signum)
            return new OpenBigInteger2(add(mag, val.mag), signum);

        int cmp = compareMagnitude(val);
        if (cmp == 0)
            return ZERO;
        int[] resultMag = (cmp > 0 ? subtract(mag, start, val.mag, val.start)
                           : subtract(val.mag, val.start, mag, start));
        resultMag = trustedStripLeadingZeroInts(resultMag);
        return new OpenBigInteger2(resultMag, cmp == signum ? 1 : -1);
    }

    /**
     * Subtracts the contents of the second int arrays (little) from the
     * first (big).  The first int array (big) must represent a larger number
     * than the second.  This method allocates the space necessary to hold the
     * answer.
     */
    private static int[] subtract(int[] big, int bigStart, int[] little, int littleStart) {
        int bigIndex = big.length;
        int result[] = new int[bigIndex];
        int littleIndex = little.length;
        long difference = 0;

        // Subtract common parts of both numbers
        while(littleIndex > littleStart) {
            difference = (big[--bigIndex] & LONG_MASK) -
                         (little[--littleIndex] & LONG_MASK) +
                         (difference >> 32);
            result[bigIndex] = (int)difference;
        }

        // Subtract remainder of longer number while borrow propagates
        boolean borrow = (difference >> 32 != 0);
        while (bigIndex > bigStart && borrow)
            borrow = ((result[--bigIndex] = big[bigIndex] - 1) == -1);

        // Copy remainder of longer number
        while (bigIndex > bigStart)
            result[--bigIndex] = big[bigIndex];

        return result;
    }

    /**
     * Returns a OpenBigInteger whose value is {@code (this * val)}.
     *
     * @param  val value to be multiplied by this OpenBigInteger.
     * @return {@code this * val}
     */
    public OpenBigInteger2 multiply(OpenBigInteger2 val) {
        if (val.signum == 0 || signum == 0)
            return ZERO;

        int[] result = multiplyToLen(mag, mag.length, 0, 
                                     val.mag, val.mag.length, null);
        result = trustedStripLeadingZeroInts(result);
        return new OpenBigInteger2(result, signum == val.signum ? 1 : -1);
    }

    /**
     * Multiplies int arrays x and y to the specified lengths and places
     * the result into z. There will be no leading zeros in the resultant array.
     */
    private int[] multiplyToLen(int[] x, int xlen, int xignore, int[] y, int ylen, int[] z) {
        int xstart = xlen - 1;
        int ystart = ylen - 1;

        if (z == null || z.length < (xlen+ ylen))
            z = new int[xlen+ylen];

        long carry = 0;
        for (int j=ystart, k=ystart+1+xstart; j>=0; j--, k--) {
            long product = (y[j] & LONG_MASK) *
                           (x[xstart] & LONG_MASK) + carry;
            z[k] = (int)product;
            carry = product >>> 32;
        }
        z[xstart] = (int)carry;

        for (int i = xstart-1; i >= xignore; i--) {
            carry = 0;
            for (int j=ystart, k=ystart+1+i; j>=0; j--, k--) {
                long product = (y[j] & LONG_MASK) *
                               (x[i] & LONG_MASK) +
                               (z[k] & LONG_MASK) + carry;
                z[k] = (int)product;
                carry = product >>> 32;
            }
            z[i] = (int)carry;
        }
        return z;
    }

    /**
     * Returns a OpenBigInteger whose value is {@code (this<sup>2</sup>)}.
     *
     * @return {@code this<sup>2</sup>}
     */
    private OpenBigInteger2 square() {
        if (signum == 0)
            return ZERO;
        int[] z = squareToLen(mag, mag.length, null);
        return new OpenBigInteger2(trustedStripLeadingZeroInts(z), 1);
    }

    /**
     * Squares the contents of the int array x. The result is placed into the
     * int array z.  The contents of x are not changed.
     */
    private static final int[] squareToLen(int[] x, int len, int[] z) {
        /*
         * The algorithm used here is adapted from Colin Plumb's C library.
         * Technique: Consider the partial products in the multiplication
         * of "abcde" by itself:
         *
         *               a  b  c  d  e
         *            *  a  b  c  d  e
         *          ==================
         *              ae be ce de ee
         *           ad bd cd dd de
         *        ac bc cc cd ce
         *     ab bb bc bd be
         *  aa ab ac ad ae
         *
         * Note that everything above the main diagonal:
         *              ae be ce de = (abcd) * e
         *           ad bd cd       = (abc) * d
         *        ac bc             = (ab) * c
         *     ab                   = (a) * b
         *
         * is a copy of everything below the main diagonal:
         *                       de
         *                 cd ce
         *           bc bd be
         *     ab ac ad ae
         *
         * Thus, the sum is 2 * (off the diagonal) + diagonal.
         *
         * This is accumulated beginning with the diagonal (which
         * consist of the squares of the digits of the input), which is then
         * divided by two, the off-diagonal added, and multiplied by two
         * again.  The low bit is simply a copy of the low bit of the
         * input, so it doesn't need special care.
         */
        int zlen = len << 1;
        if (z == null || z.length < zlen)
            z = new int[zlen];

        // Store the squares, right shifted one bit (i.e., divided by 2)
        int lastProductLowWord = 0;
        for (int j=0, i=0; j<len; j++) {
            long piece = (x[j] & LONG_MASK);
            long product = piece * piece;
            z[i++] = (lastProductLowWord << 31) | (int)(product >>> 33);
            z[i++] = (int)(product >>> 1);
            lastProductLowWord = (int)product;
        }

        // Add in off-diagonal sums
        for (int i=len, offset=1; i>0; i--, offset+=2) {
            int t = x[i-1];
            t = mulAdd(z, x, offset, i-1, t);
            addOne(z, offset-1, i, t);
        }

        // Shift back up and set low bit
        primitiveLeftShift(z, zlen, 1);
        z[zlen-1] |= x[len-1] & 1;

        return z;
    }

    /**
     * Returns a OpenBigInteger whose value is <tt>(this<sup>exponent</sup>)</tt>.
     * Note that {@code exponent} is an integer rather than a OpenBigInteger.
     *
     * @param  exponent exponent to which this OpenBigInteger is to be raised.
     * @return <tt>this<sup>exponent</sup></tt>
     * @throws ArithmeticException {@code exponent} is negative.  (This would
     *         cause the operation to yield a non-integer value.)
     */
    public OpenBigInteger2 pow(int exponent) {
        if (exponent < 0)
            throw new ArithmeticException("Negative exponent");
        if (signum==0)
            return (exponent==0 ? ONE : this);

        // Perform exponentiation using repeated squaring trick
        int newSign = (signum<0 && (exponent&1)==1 ? -1 : 1);
        int[] baseToPow2 = this.mag;
        int[] result = {1};

        while (exponent != 0) {
            if ((exponent & 1)==1) {
                result = multiplyToLen(result, result.length, 0,
                                       baseToPow2, baseToPow2.length, null);
                result = trustedStripLeadingZeroInts(result);
            }
            if ((exponent >>>= 1) != 0) {
                baseToPow2 = squareToLen(baseToPow2, baseToPow2.length, null);
                baseToPow2 = trustedStripLeadingZeroInts(baseToPow2);
            }
        }
        return new OpenBigInteger2(result, newSign);
    }

    /**
     * Package private method to return bit length for an integer.
     */
    static int bitLengthForInt(int n) {
        return 32 - Integer.numberOfLeadingZeros(n);
    }

    /**
     * Left shift int array a up to len by n bits. Returns the array that
     * results from the shift since space may have to be reallocated.
     */
    private static int[] leftShift(int[] a, int len, int n) {
        int nInts = n >>> 5;
        int nBits = n&0x1F;
        int bitsInHighWord = bitLengthForInt(a[0]);

        // If shift can be done without recopy, do so
        if (n <= (32-bitsInHighWord)) {
            primitiveLeftShift(a, len, nBits);
            return a;
        } else { // Array must be resized
            if (nBits <= (32-bitsInHighWord)) {
                int result[] = new int[nInts+len];
                for (int i=0; i<len; i++)
                    result[i] = a[i];
                primitiveLeftShift(result, result.length, nBits);
                return result;
            } else {
                int result[] = new int[nInts+len+1];
                for (int i=0; i<len; i++)
                    result[i] = a[i];
                primitiveRightShift(result, result.length, 32 - nBits);
                return result;
            }
        }
    }

    // shifts a up to len right n bits assumes no leading zeros, 0<n<32
    static void primitiveRightShift(int[] a, int len, int n) {
        int n2 = 32 - n;
        for (int i=len-1, c=a[i]; i>0; i--) {
            int b = c;
            c = a[i-1];
            a[i] = (c << n2) | (b >>> n);
        }
        a[0] >>>= n;
    }

    // shifts a up to len left n bits assumes no leading zeros, 0<=n<32
    static void primitiveLeftShift(int[] a, int len, int n) {
        if (len == 0 || n == 0)
            return;

        int n2 = 32 - n;
        for (int i=0, c=a[i], m=i+len-1; i<m; i++) {
            int b = c;
            c = a[i+1];
            a[i] = (b << n) | (c >>> n2);
        }
        a[len-1] <<= n;
    }

    /**
     * Calculate bitlength of contents of the first len elements an int array,
     * assuming there are no leading zero ints.
     */
    private static int bitLength(int[] val, int len) {
        if (len == 0)
            return 0;
        return ((len - 1) << 5) + bitLengthForInt(val[0]);
    }

    /**
     * Returns a OpenBigInteger whose value is the absolute value of this
     * OpenBigInteger.
     *
     * @return {@code abs(this)}
     */
    public OpenBigInteger2 abs() {
        return (signum >= 0 ? this : this.negate());
    }

    /**
     * Returns a OpenBigInteger whose value is {@code (-this)}.
     *
     * @return {@code -this}
     */
    public OpenBigInteger2 negate() {
        return new OpenBigInteger2(this.mag, -this.signum);
    }

    /**
     * Returns the signum function of this OpenBigInteger.
     *
     * @return -1, 0 or 1 as the value of this OpenBigInteger is negative, zero or
     *         positive.
     */
    public int signum() {
        return this.signum;
    }

    // Modular Arithmetic Operations

    static int[] bnExpModThreshTable = {7, 25, 81, 241, 673, 1793,
                                                Integer.MAX_VALUE}; // Sentinel

    /**
     * Montgomery reduce n, modulo mod.  This reduces modulo mod and divides
     * by 2^(32*mlen). Adapted from Colin Plumb's C library.
     */
    private static int[] montReduce(int[] n, int[] mod, int mlen, int inv) {
        int c=0;
        int len = mlen;
        int offset=0;

        do {
            int nEnd = n[n.length-1-offset];
            int carry = mulAdd(n, mod, offset, mlen, inv * nEnd);
            c += addOne(n, offset, mlen, carry);
            offset++;
        } while(--len > 0);

        while(c>0)
            c += subN(n, mod, mlen);

        while (intArrayCmpToLen(n, mod, mlen) >= 0)
            subN(n, mod, mlen);

        return n;
    }


    /*
     * Returns -1, 0 or +1 as big-endian unsigned int array arg1 is less than,
     * equal to, or greater than arg2 up to length len.
     */
    private static int intArrayCmpToLen(int[] arg1, int[] arg2, int len) {
        for (int i=0; i<len; i++) {
            long b1 = arg1[i] & LONG_MASK;
            long b2 = arg2[i] & LONG_MASK;
            if (b1 < b2)
                return -1;
            if (b1 > b2)
                return 1;
        }
        return 0;
    }

    /**
     * Subtracts two numbers of same length, returning borrow.
     */
    private static int subN(int[] a, int[] b, int len) {
        long sum = 0;

        while(--len >= 0) {
            sum = (a[len] & LONG_MASK) -
                 (b[len] & LONG_MASK) + (sum >> 32);
            a[len] = (int)sum;
        }

        return (int)(sum >> 32);
    }

    /**
     * Multiply an array by one word k and add to result, return the carry
     */
    static int mulAdd(int[] out, int[] in, int offset, int len, int k) {
        long kLong = k & LONG_MASK;
        long carry = 0;

        offset = out.length-offset - 1;
        for (int j=len-1; j >= 0; j--) {
            long product = (in[j] & LONG_MASK) * kLong +
                           (out[offset] & LONG_MASK) + carry;
            out[offset--] = (int)product;
            carry = product >>> 32;
        }
        return (int)carry;
    }

    /**
     * Add one word to the number a mlen words into a. Return the resulting
     * carry.
     */
    static int addOne(int[] a, int offset, int mlen, int carry) {
        offset = a.length-1-mlen-offset;
        long t = (a[offset] & LONG_MASK) + (carry & LONG_MASK);

        a[offset] = (int)t;
        if ((t >>> 32) == 0)
            return 0;
        while (--mlen >= 0) {
            if (--offset < 0) { // Carry out of number
                return 1;
            } else {
                a[offset]++;
                if (a[offset] != 0)
                    return 0;
            }
        }
        return 1;
    }

    /**
     * Returns a OpenBigInteger whose value is (this ** exponent) mod (2**p)
     */
    private OpenBigInteger2 modPow2(OpenBigInteger2 exponent, int p) {
        /*
         * Perform exponentiation using repeated squaring trick, chopping off
         * high order bits as indicated by modulus.
         */
        OpenBigInteger2 result = valueOf(1);
        OpenBigInteger2 baseToPow2 = this.mod2(p);
        int expOffset = 0;

        int limit = exponent.bitLength();

        if (this.testBit(0))
           limit = (p-1) < limit ? (p-1) : limit;

        while (expOffset < limit) {
            if (exponent.testBit(expOffset))
                result = result.multiply(baseToPow2).mod2(p);
            expOffset++;
            if (expOffset < limit)
                baseToPow2 = baseToPow2.square().mod2(p);
        }

        return result;
    }

    /**
     * Returns a OpenBigInteger whose value is this mod(2**p).
     * Assumes that this {@code OpenBigInteger >= 0} and {@code p > 0}.
     */
    private OpenBigInteger2 mod2(int p) {
        if (bitLength() <= p)
            return this;

        // Copy remaining ints of mag
        int numInts = (p + 31) >>> 5;
        int[] mag = new int[numInts];
        for (int i=0; i<numInts; i++)
            mag[i] = this.mag[i + (this.mag.length - numInts)];

        // Mask out any excess bits
        int excessBits = (numInts << 5) - p;
        mag[0] &= (1L << (32-excessBits)) - 1;

        return (mag[0]==0 ? new OpenBigInteger2(1, mag) : new OpenBigInteger2(mag, 1));
    }

    // Shift Operations

    /**
     * Returns a OpenBigInteger whose value is {@code (this << n)}.
     * The shift distance, {@code n}, may be negative, in which case
     * this method performs a right shift.
     * (Computes <tt>floor(this * 2<sup>n</sup>)</tt>.)
     *
     * @param  n shift distance, in bits.
     * @return {@code this << n}
     * @throws ArithmeticException if the shift distance is {@code
     *         Integer.MIN_VALUE}.
     * @see #shiftRight
     */
    public OpenBigInteger2 shiftLeft(int n) {
        if (signum == 0)
            return ZERO;
        if (n==0)
            return this;
        if (n<0) {
            if (n == Integer.MIN_VALUE) {
                throw new ArithmeticException("Shift distance of Integer.MIN_VALUE not supported.");
            } else {
                return shiftRight(-n);
            }
        }

        int nInts = n >>> 5;
        int nBits = n & 0x1f;
        int magLen = mag.length;
        int newMag[] = null;

        if (nBits == 0) {
            newMag = new int[magLen + nInts];
            for (int i=0; i<magLen; i++)
                newMag[i] = mag[i];
        } else {
            int i = 0;
            int nBits2 = 32 - nBits;
            int highBits = mag[0] >>> nBits2;
            if (highBits != 0) {
                newMag = new int[magLen + nInts + 1];
                newMag[i++] = highBits;
            } else {
                newMag = new int[magLen + nInts];
            }
            int j=0;
            while (j < magLen-1)
                newMag[i++] = mag[j++] << nBits | mag[j] >>> nBits2;
            newMag[i] = mag[j] << nBits;
        }

        return new OpenBigInteger2(newMag, signum);
    }

    /**
     * Returns a OpenBigInteger whose value is {@code (this >> n)}.  Sign
     * extension is performed.  The shift distance, {@code n}, may be
     * negative, in which case this method performs a left shift.
     * (Computes <tt>floor(this / 2<sup>n</sup>)</tt>.)
     *
     * @param  n shift distance, in bits.
     * @return {@code this >> n}
     * @throws ArithmeticException if the shift distance is {@code
     *         Integer.MIN_VALUE}.
     * @see #shiftLeft
     */
    public OpenBigInteger2 shiftRight(int n) {
        if (n==0)
            return this;
        if (n<0) {
            if (n == Integer.MIN_VALUE) {
                throw new ArithmeticException("Shift distance of Integer.MIN_VALUE not supported.");
            } else {
                return shiftLeft(-n);
            }
        }

        int nInts = n >>> 5;
        int nBits = n & 0x1f;
        int magLen = mag.length;
        int newMag[] = null;

        // Special case: entire contents shifted off the end
        if (nInts >= magLen)
            return (signum >= 0 ? ZERO : negConst[1]);

        if (nBits == 0) {
            int newMagLen = magLen - nInts;
            newMag = new int[newMagLen];
            for (int i=0; i<newMagLen; i++)
                newMag[i] = mag[i];
        } else {
            int i = 0;
            int highBits = mag[0] >>> nBits;
            if (highBits != 0) {
                newMag = new int[magLen - nInts];
                newMag[i++] = highBits;
            } else {
                newMag = new int[magLen - nInts -1];
            }

            int nBits2 = 32 - nBits;
            int j=0;
            while (j < magLen - nInts - 1)
                newMag[i++] = (mag[j++] << nBits2) | (mag[j] >>> nBits);
        }

        if (signum < 0) {
            // Find out whether any one-bits were shifted off the end.
            boolean onesLost = false;
            for (int i=magLen-1, j=magLen-nInts; i>=j && !onesLost; i--)
                onesLost = (mag[i] != 0);
            if (!onesLost && nBits != 0)
                onesLost = (mag[magLen - nInts - 1] << (32 - nBits) != 0);

            if (onesLost)
                newMag = javaIncrement(newMag);
        }

        return new OpenBigInteger2(newMag, signum);
    }

    int[] javaIncrement(int[] val) {
        int lastSum = 0;
        for (int i=val.length-1;  i >= 0 && lastSum == 0; i--)
            lastSum = (val[i] += 1);
        if (lastSum == 0) {
            val = new int[val.length+1];
            val[0] = 1;
        }
        return val;
    }

    // Bitwise Operations

    /**
     * Returns a OpenBigInteger whose value is {@code (this & val)}.  (This
     * method returns a negative OpenBigInteger if and only if this and val are
     * both negative.)
     *
     * @param val value to be AND'ed with this OpenBigInteger.
     * @return {@code this & val}
     */
    public OpenBigInteger2 and(OpenBigInteger2 val) {
        int[] result = new int[Math.max(intLength(), val.intLength())];
        for (int i=0; i<result.length; i++)
            result[i] = (getInt(result.length-i-1)
                         & val.getInt(result.length-i-1));

        return valueOf(result);
    }

    /**
     * Returns a OpenBigInteger whose value is {@code (this | val)}.  (This method
     * returns a negative OpenBigInteger if and only if either this or val is
     * negative.)
     *
     * @param val value to be OR'ed with this OpenBigInteger.
     * @return {@code this | val}
     */
    public OpenBigInteger2 or(OpenBigInteger2 val) {
        int[] result = new int[Math.max(intLength(), val.intLength())];
        for (int i=0; i<result.length; i++)
            result[i] = (getInt(result.length-i-1)
                         | val.getInt(result.length-i-1));

        return valueOf(result);
    }

    /**
     * Returns a OpenBigInteger whose value is {@code (this ^ val)}.  (This method
     * returns a negative OpenBigInteger if and only if exactly one of this and
     * val are negative.)
     *
     * @param val value to be XOR'ed with this OpenBigInteger.
     * @return {@code this ^ val}
     */
    public OpenBigInteger2 xor(OpenBigInteger2 val) {
        int[] result = new int[Math.max(intLength(), val.intLength())];
        for (int i=0; i<result.length; i++)
            result[i] = (getInt(result.length-i-1)
                         ^ val.getInt(result.length-i-1));

        return valueOf(result);
    }

    /**
     * Returns a OpenBigInteger whose value is {@code (~this)}.  (This method
     * returns a negative value if and only if this OpenBigInteger is
     * non-negative.)
     *
     * @return {@code ~this}
     */
    public OpenBigInteger2 not() {
        int[] result = new int[intLength()];
        for (int i=0; i<result.length; i++)
            result[i] = ~getInt(result.length-i-1);

        return valueOf(result);
    }

    /**
     * Returns a OpenBigInteger whose value is {@code (this & ~val)}.  This
     * method, which is equivalent to {@code and(val.not())}, is provided as
     * a convenience for masking operations.  (This method returns a negative
     * OpenBigInteger if and only if {@code this} is negative and {@code val} is
     * positive.)
     *
     * @param val value to be complemented and AND'ed with this OpenBigInteger.
     * @return {@code this & ~val}
     */
    public OpenBigInteger2 andNot(OpenBigInteger2 val) {
        int[] result = new int[Math.max(intLength(), val.intLength())];
        for (int i=0; i<result.length; i++)
            result[i] = (getInt(result.length-i-1)
                         & ~val.getInt(result.length-i-1));

        return valueOf(result);
    }


    // Single Bit Operations

    /**
     * Returns {@code true} if and only if the designated bit is set.
     * (Computes {@code ((this & (1<<n)) != 0)}.)
     *
     * @param  n index of bit to test.
     * @return {@code true} if and only if the designated bit is set.
     * @throws ArithmeticException {@code n} is negative.
     */
    public boolean testBit(int n) {
        if (n<0)
            throw new ArithmeticException("Negative bit address");

        return (getInt(n >>> 5) & (1 << (n & 31))) != 0;
    }

    /**
     * Returns a OpenBigInteger whose value is equivalent to this OpenBigInteger
     * with the designated bit set.  (Computes {@code (this | (1<<n))}.)
     *
     * @param  n index of bit to set.
     * @return {@code this | (1<<n)}
     * @throws ArithmeticException {@code n} is negative.
     */
    public OpenBigInteger2 setBit(int n) {
        if (n<0)
            throw new ArithmeticException("Negative bit address");

        int intNum = n >>> 5;
        int[] result = new int[Math.max(intLength(), intNum+2)];

        for (int i=0; i<result.length; i++)
            result[result.length-i-1] = getInt(i);

        result[result.length-intNum-1] |= (1 << (n & 31));

        return valueOf(result);
    }

    /**
     * Returns a OpenBigInteger whose value is equivalent to this OpenBigInteger
     * with the designated bit cleared.
     * (Computes {@code (this & ~(1<<n))}.)
     *
     * @param  n index of bit to clear.
     * @return {@code this & ~(1<<n)}
     * @throws ArithmeticException {@code n} is negative.
     */
    public OpenBigInteger2 clearBit(int n) {
        if (n<0)
            throw new ArithmeticException("Negative bit address");

        int intNum = n >>> 5;
        int[] result = new int[Math.max(intLength(), ((n + 1) >>> 5) + 1)];

        for (int i=0; i<result.length; i++)
            result[result.length-i-1] = getInt(i);

        result[result.length-intNum-1] &= ~(1 << (n & 31));

        return valueOf(result);
    }

    /**
     * Returns a OpenBigInteger whose value is equivalent to this OpenBigInteger
     * with the designated bit flipped.
     * (Computes {@code (this ^ (1<<n))}.)
     *
     * @param  n index of bit to flip.
     * @return {@code this ^ (1<<n)}
     * @throws ArithmeticException {@code n} is negative.
     */
    public OpenBigInteger2 flipBit(int n) {
        if (n<0)
            throw new ArithmeticException("Negative bit address");

        int intNum = n >>> 5;
        int[] result = new int[Math.max(intLength(), intNum+2)];

        for (int i=0; i<result.length; i++)
            result[result.length-i-1] = getInt(i);

        result[result.length-intNum-1] ^= (1 << (n & 31));

        return valueOf(result);
    }

    /**
     * Returns the index of the rightmost (lowest-order) one bit in this
     * OpenBigInteger (the number of zero bits to the right of the rightmost
     * one bit).  Returns -1 if this OpenBigInteger contains no one bits.
     * (Computes {@code (this==0? -1 : log2(this & -this))}.)
     *
     * @return index of the rightmost one bit in this OpenBigInteger.
     */
    public int getLowestSetBit() {
        @SuppressWarnings("deprecation") int lsb = lowestSetBit - 2;
        if (lsb == -2) {  // lowestSetBit not initialized yet
            lsb = 0;
            if (signum == 0) {
                lsb -= 1;
            } else {
                // Search for lowest order nonzero int
                int i,b;
                for (i=0; (b = getInt(i))==0; i++)
                    ;
                lsb += (i << 5) + Integer.numberOfTrailingZeros(b);
            }
            lowestSetBit = lsb + 2;
        }
        return lsb;
    }


    // Miscellaneous Bit Operations

    /**
     * Returns the number of bits in the minimal two's-complement
     * representation of this OpenBigInteger, <i>excluding</i> a sign bit.
     * For positive OpenBigIntegers, this is equivalent to the number of bits in
     * the ordinary binary representation.  (Computes
     * {@code (ceil(log2(this < 0 ? -this : this+1)))}.)
     *
     * @return number of bits in the minimal two's-complement
     *         representation of this OpenBigInteger, <i>excluding</i> a sign bit.
     */
    public int bitLength() {
        @SuppressWarnings("deprecation") int n = bitLength - 1;
        if (n == -1) { // bitLength not initialized yet
            int[] m = mag;
            int len = m.length - start;
            if (len == 0) {
                n = 0; // offset by one to initialize
            }  else {
                // Calculate the bit length of the magnitude
                int magBitLength = ((len - 1) << 5) + bitLengthForInt(mag[start]);
                 if (signum < 0) {
                     // Check if magnitude is a power of two
                     boolean pow2 = (Integer.bitCount(mag[start]) == 1);
                     for(int i=1; i< len && pow2; i++)
                         pow2 = (mag[i+start] == 0);

                     n = (pow2 ? magBitLength -1 : magBitLength);
                 } else {
                     n = magBitLength;
                 }
            }
            bitLength = n + 1;
        }
        return n;
    }

    /**
     * Returns the number of bits in the two's complement representation
     * of this OpenBigInteger that differ from its sign bit.  This method is
     * useful when implementing bit-vector style sets atop OpenBigIntegers.
     *
     * @return number of bits in the two's complement representation
     *         of this OpenBigInteger that differ from its sign bit.
     */
    public int bitCount() {
        @SuppressWarnings("deprecation") int bc = bitCount - 1;
        if (bc == -1) {  // bitCount not initialized yet
            bc = 0;      // offset by one to initialize
            // Count the bits in the magnitude
            for (int i=0; i<mag.length; i++)
                bc += Integer.bitCount(mag[i]);
            if (signum < 0) {
                // Count the trailing zeros in the magnitude
                int magTrailingZeroCount = 0, j;
                for (j=mag.length-1; mag[j]==0; j--)
                    magTrailingZeroCount += 32;
                magTrailingZeroCount += Integer.numberOfTrailingZeros(mag[j]);
                bc += magTrailingZeroCount - 1;
            }
            bitCount = bc + 1;
        }
        return bc;
    }

    // Comparison Operations

    /**
     * Compares this OpenBigInteger with the specified OpenBigInteger.  This
     * method is provided in preference to individual methods for each
     * of the six boolean comparison operators ({@literal <}, ==,
     * {@literal >}, {@literal >=}, !=, {@literal <=}).  The suggested
     * idiom for performing these comparisons is: {@code
     * (x.compareTo(y)} &lt;<i>op</i>&gt; {@code 0)}, where
     * &lt;<i>op</i>&gt; is one of the six comparison operators.
     *
     * @param  val OpenBigInteger to which this OpenBigInteger is to be compared.
     * @return -1, 0 or 1 as this OpenBigInteger is numerically less than, equal
     *         to, or greater than {@code val}.
     */
    public int compareTo(OpenBigInteger2 val) {
        if (signum == val.signum) {
            switch (signum) {
            case 1:
                return compareMagnitude(val);
            case -1:
                return val.compareMagnitude(this);
            default:
                return 0;
            }
        }
        return signum > val.signum ? 1 : -1;
    }

    /**
     * Compares the magnitude array of this OpenBigInteger with the specified
     * OpenBigInteger's. This is the version of compareTo ignoring sign.
     *
     * @param val OpenBigInteger whose magnitude array to be compared.
     * @return -1, 0 or 1 as this magnitude array is less than, equal to or
     *         greater than the magnitude aray for the specified OpenBigInteger's.
     */
    final int compareMagnitude(OpenBigInteger2 val) {
        int[] m1 = mag;
        int len1 = m1.length;
        int[] m2 = val.mag;
        int len2 = m2.length;
        if (len1 < len2)
            return -1;
        if (len1 > len2)
            return 1;
        for (int i = 0; i < len1; i++) {
            int a = m1[i];
            int b = m2[i];
            if (a != b)
                return ((a & LONG_MASK) < (b & LONG_MASK)) ? -1 : 1;
        }
        return 0;
    }

    /**
     * Compares this OpenBigInteger with the specified Object for equality.
     *
     * @param  x Object to which this OpenBigInteger is to be compared.
     * @return {@code true} if and only if the specified Object is a
     *         OpenBigInteger whose value is numerically equal to this OpenBigInteger.
     */
    public boolean equals(Object x) {
        // This test is just an optimization, which may or may not help
        if (x == this)
            return true;

        if (!(x instanceof OpenBigInteger2))
            return false;

        OpenBigInteger2 xInt = (OpenBigInteger2) x;
        if (xInt.signum != signum)
            return false;

        int[] m = mag;
        int len = m.length;
        int[] xm = xInt.mag;
        if (len != xm.length)
            return false;

        for (int i = 0; i < len; i++)
            if (xm[i] != m[i])
                return false;

        return true;
    }

    /**
     * Returns the minimum of this OpenBigInteger and {@code val}.
     *
     * @param  val value with which the minimum is to be computed.
     * @return the OpenBigInteger whose value is the lesser of this OpenBigInteger and
     *         {@code val}.  If they are equal, either may be returned.
     */
    public OpenBigInteger2 min(OpenBigInteger2 val) {
        return (compareTo(val)<0 ? this : val);
    }

    /**
     * Returns the maximum of this OpenBigInteger and {@code val}.
     *
     * @param  val value with which the maximum is to be computed.
     * @return the OpenBigInteger whose value is the greater of this and
     *         {@code val}.  If they are equal, either may be returned.
     */
    public OpenBigInteger2 max(OpenBigInteger2 val) {
        return (compareTo(val)>0 ? this : val);
    }


    // Hash Function

    /**
     * Returns the hash code for this OpenBigInteger.
     *
     * @return hash code for this OpenBigInteger.
     */
    public int hashCode() {
        int hashCode = 0;

        for (int i=0; i<mag.length; i++)
            hashCode = (int)(31*hashCode + (mag[i] & LONG_MASK));

        return hashCode * signum;
    }

//    /**
//     * Returns the String representation of this OpenBigInteger in the
//     * given radix.  If the radix is outside the range from {@link
//     * Character#MIN_RADIX} to {@link Character#MAX_RADIX} inclusive,
//     * it will default to 10 (as is the case for
//     * {@code Integer.toString}).  The digit-to-character mapping
//     * provided by {@code Character.forDigit} is used, and a minus
//     * sign is prepended if appropriate.  (This representation is
//     * compatible with the {@link #OpenBigInteger(String, int) (String,
//     * int)} constructor.)
//     *
//     * @param  radix  radix of the String representation.
//     * @return String representation of this OpenBigInteger in the given radix.
//     * @see    Integer#toString
//     * @see    Character#forDigit
//     * @see    #OpenBigInteger(java.lang.String, int)
//     */
//    public String toString(int radix) {
//        if (signum == 0)
//            return "0";
//        if (radix < Character.MIN_RADIX || radix > Character.MAX_RADIX)
//            radix = 10;
//
//        // Compute upper bound on number of digit groups and allocate space
//        int maxNumDigitGroups = (4*mag.length + 6)/7;
//        String digitGroup[] = new String[maxNumDigitGroups];
//
//        // Translate number to string, a digit group at a time
//        OpenBigInteger tmp = this.abs();
//        int numGroups = 0;
//        while (tmp.signum != 0) {
//            OpenBigInteger d = longRadix[radix];
//
//            MutableOpenBigInteger q = new MutableOpenBigInteger(),
//                              a = new MutableOpenBigInteger(tmp.mag),
//                              b = new MutableOpenBigInteger(d.mag);
//            MutableOpenBigInteger r = a.divide(b, q);
//            OpenBigInteger q2 = q.toOpenBigInteger(tmp.signum * d.signum);
//            OpenBigInteger r2 = r.toOpenBigInteger(tmp.signum * d.signum);
//
//            digitGroup[numGroups++] = Long.toString(r2.longValue(), radix);
//            tmp = q2;
//        }
//
//        // Put sign (if any) and first digit group into result buffer
//        StringBuilder buf = new StringBuilder(numGroups*digitsPerLong[radix]+1);
//        if (signum<0)
//            buf.append('-');
//        buf.append(digitGroup[numGroups-1]);
//
//        // Append remaining digit groups padded with leading zeros
//        for (int i=numGroups-2; i>=0; i--) {
//            // Prepend (any) leading zeros for this digit group
//            int numLeadingZeros = digitsPerLong[radix]-digitGroup[i].length();
//            if (numLeadingZeros != 0)
//                buf.append(zeros[numLeadingZeros]);
//            buf.append(digitGroup[i]);
//        }
//        return buf.toString();
//    }

    /* zero[i] is a string of i consecutive zeros. */
    private static String zeros[] = new String[64];
    static {
        zeros[63] =
            "000000000000000000000000000000000000000000000000000000000000000";
        for (int i=0; i<63; i++)
            zeros[i] = zeros[63].substring(0, i);
    }

//    /**
//     * Returns the decimal String representation of this OpenBigInteger.
//     * The digit-to-character mapping provided by
//     * {@code Character.forDigit} is used, and a minus sign is
//     * prepended if appropriate.  (This representation is compatible
//     * with the {@link #OpenBigInteger(String) (String)} constructor, and
//     * allows for String concatenation with Java's + operator.)
//     *
//     * @return decimal String representation of this OpenBigInteger.
//     * @see    Character#forDigit
//     * @see    #OpenBigInteger(java.lang.String)
//     */
//    public String toString() {
//        return toString(10);
//    }

    /**
     * Returns a byte array containing the two's-complement
     * representation of this OpenBigInteger.  The byte array will be in
     * <i>big-endian</i> byte-order: the most significant byte is in
     * the zeroth element.  The array will contain the minimum number
     * of bytes required to represent this OpenBigInteger, including at
     * least one sign bit, which is {@code (ceil((this.bitLength() +
     * 1)/8))}.  (This representation is compatible with the
     * {@link #OpenBigInteger(byte[]) (byte[])} constructor.)
     *
     * @return a byte array containing the two's-complement representation of
     *         this OpenBigInteger.
     * @see    #OpenBigInteger(byte[])
     */
    public byte[] toByteArray() {
        int byteLen = bitLength()/8 + 1;
        byte[] byteArray = new byte[byteLen];

        for (int i=byteLen-1, bytesCopied=4, nextInt=0, intIndex=0; i>=0; i--) {
            if (bytesCopied == 4) {
                nextInt = getInt(intIndex++);
                bytesCopied = 1;
            } else {
                nextInt >>>= 8;
                bytesCopied++;
            }
            byteArray[i] = (byte)nextInt;
        }
        return byteArray;
    }

    /**
     * Converts this OpenBigInteger to an {@code int}.  This
     * conversion is analogous to a
     * <i>narrowing primitive conversion</i> from {@code long} to
     * {@code int} as defined in section 5.1.3 of
     * <cite>The Java&trade; Language Specification</cite>:
     * if this OpenBigInteger is too big to fit in an
     * {@code int}, only the low-order 32 bits are returned.
     * Note that this conversion can lose information about the
     * overall magnitude of the OpenBigInteger value as well as return a
     * result with the opposite sign.
     *
     * @return this OpenBigInteger converted to an {@code int}.
     */
    public int intValue() {
        int result = 0;
        result = getInt(0);
        return result;
    }

    /**
     * Converts this OpenBigInteger to a {@code long}.  This
     * conversion is analogous to a
     * <i>narrowing primitive conversion</i> from {@code long} to
     * {@code int} as defined in section 5.1.3 of
     * <cite>The Java&trade; Language Specification</cite>:
     * if this OpenBigInteger is too big to fit in a
     * {@code long}, only the low-order 64 bits are returned.
     * Note that this conversion can lose information about the
     * overall magnitude of the OpenBigInteger value as well as return a
     * result with the opposite sign.
     *
     * @return this OpenBigInteger converted to a {@code long}.
     */
    public long longValue() {
        long result = 0;

        for (int i=1; i>=0; i--)
            result = (result << 32) + (getInt(i) & LONG_MASK);
        return result;
    }

    /**
     * Converts this OpenBigInteger to a {@code float}.  This
     * conversion is similar to the
     * <i>narrowing primitive conversion</i> from {@code double} to
     * {@code float} as defined in section 5.1.3 of
     * <cite>The Java&trade; Language Specification</cite>:
     * if this OpenBigInteger has too great a magnitude
     * to represent as a {@code float}, it will be converted to
     * {@link Float#NEGATIVE_INFINITY} or {@link
     * Float#POSITIVE_INFINITY} as appropriate.  Note that even when
     * the return value is finite, this conversion can lose
     * information about the precision of the OpenBigInteger value.
     *
     * @return this OpenBigInteger converted to a {@code float}.
     */
    public float floatValue() {
        // Somewhat inefficient, but guaranteed to work.
        return Float.parseFloat(this.toString());
    }

    /**
     * Converts this OpenBigInteger to a {@code double}.  This
     * conversion is similar to the
     * <i>narrowing primitive conversion</i> from {@code double} to
     * {@code float} as defined in section 5.1.3 of
     * <cite>The Java&trade; Language Specification</cite>:
     * if this OpenBigInteger has too great a magnitude
     * to represent as a {@code double}, it will be converted to
     * {@link Double#NEGATIVE_INFINITY} or {@link
     * Double#POSITIVE_INFINITY} as appropriate.  Note that even when
     * the return value is finite, this conversion can lose
     * information about the precision of the OpenBigInteger value.
     *
     * @return this OpenBigInteger converted to a {@code double}.
     */
    public double doubleValue() {
        // Somewhat inefficient, but guaranteed to work.
        return Double.parseDouble(this.toString());
    }

    /**
     * Returns a copy of the input array stripped of any leading zero bytes.
     */
    private static int[] stripLeadingZeroInts(int val[]) {
        int vlen = val.length;
        int keep;

        // Find first nonzero byte
        for (keep = 0; keep < vlen && val[keep] == 0; keep++)
            ;
        return java.util.Arrays.copyOfRange(val, keep, vlen);
    }

    /**
     * Returns the input array stripped of any leading zero bytes.
     * Since the source is trusted the copying may be skipped.
     */
    private static int[] trustedStripLeadingZeroInts(int val[]) {
        int vlen = val.length;
        int keep;

        // Find first nonzero byte
        for (keep = 0; keep < vlen && val[keep] == 0; keep++)
            ;
        return keep == 0 ? val : java.util.Arrays.copyOfRange(val, keep, vlen);
    }

    /**
     * Returns a copy of the input array stripped of any leading zero bytes.
     */
    private static int[] stripLeadingZeroBytes(byte a[]) {
        int byteLength = a.length;
        int keep;

        // Find first nonzero byte
        for (keep = 0; keep < byteLength && a[keep]==0; keep++)
            ;

        // Allocate new array and copy relevant part of input array
        int intLength = ((byteLength - keep) + 3) >>> 2;
        int[] result = new int[intLength];
        int b = byteLength - 1;
        for (int i = intLength-1; i >= 0; i--) {
            result[i] = a[b--] & 0xff;
            int bytesRemaining = b - keep + 1;
            int bytesToTransfer = Math.min(3, bytesRemaining);
            for (int j=8; j <= (bytesToTransfer << 3); j += 8)
                result[i] |= ((a[b--] & 0xff) << j);
        }
        return result;
    }

    /**
     * Takes an array a representing a negative 2's-complement number and
     * returns the minimal (no leading zero bytes) unsigned whose value is -a.
     */
    private static int[] makePositive(byte a[]) {
        int keep, k;
        int byteLength = a.length;

        // Find first non-sign (0xff) byte of input
        for (keep=0; keep<byteLength && a[keep]==-1; keep++)
            ;


        /* Allocate output array.  If all non-sign bytes are 0x00, we must
         * allocate space for one extra output byte. */
        for (k=keep; k<byteLength && a[k]==0; k++)
            ;

        int extraByte = (k==byteLength) ? 1 : 0;
        int intLength = ((byteLength - keep + extraByte) + 3)/4;
        int result[] = new int[intLength];

        /* Copy one's complement of input into output, leaving extra
         * byte (if it exists) == 0x00 */
        int b = byteLength - 1;
        for (int i = intLength-1; i >= 0; i--) {
            result[i] = a[b--] & 0xff;
            int numBytesToTransfer = Math.min(3, b-keep+1);
            if (numBytesToTransfer < 0)
                numBytesToTransfer = 0;
            for (int j=8; j <= 8*numBytesToTransfer; j += 8)
                result[i] |= ((a[b--] & 0xff) << j);

            // Mask indicates which bits must be complemented
            int mask = -1 >>> (8*(3-numBytesToTransfer));
            result[i] = ~result[i] & mask;
        }

        // Add one to one's complement to generate two's complement
        for (int i=result.length-1; i>=0; i--) {
            result[i] = (int)((result[i] & LONG_MASK) + 1);
            if (result[i] != 0)
                break;
        }

        return result;
    }

    /**
     * Takes an array a representing a negative 2's-complement number and
     * returns the minimal (no leading zero ints) unsigned whose value is -a.
     */
    private static int[] makePositive(int a[]) {
        int keep, j;

        // Find first non-sign (0xffffffff) int of input
        for (keep=0; keep<a.length && a[keep]==-1; keep++)
            ;

        /* Allocate output array.  If all non-sign ints are 0x00, we must
         * allocate space for one extra output int. */
        for (j=keep; j<a.length && a[j]==0; j++)
            ;
        int extraInt = (j==a.length ? 1 : 0);
        int result[] = new int[a.length - keep + extraInt];

        /* Copy one's complement of input into output, leaving extra
         * int (if it exists) == 0x00 */
        for (int i = keep; i<a.length; i++)
            result[i - keep + extraInt] = ~a[i];

        // Add one to one's complement to generate two's complement
        for (int i=result.length-1; ++result[i]==0; i--)
            ;

        return result;
    }

    /*
     * The following two arrays are used for fast String conversions.  Both
     * are indexed by radix.  The first is the number of digits of the given
     * radix that can fit in a Java long without "going negative", i.e., the
     * highest integer n such that radix**n < 2**63.  The second is the
     * "long radix" that tears each number into "long digits", each of which
     * consists of the number of digits in the corresponding element in
     * digitsPerLong (longRadix[i] = i**digitPerLong[i]).  Both arrays have
     * nonsense values in their 0 and 1 elements, as radixes 0 and 1 are not
     * used.
     */
    private static int digitsPerLong[] = {0, 0,
        62, 39, 31, 27, 24, 22, 20, 19, 18, 18, 17, 17, 16, 16, 15, 15, 15, 14,
        14, 14, 14, 13, 13, 13, 13, 13, 13, 12, 12, 12, 12, 12, 12, 12, 12};

    private static OpenBigInteger2 longRadix[] = {null, null,
        valueOf(0x4000000000000000L), valueOf(0x383d9170b85ff80bL),
        valueOf(0x4000000000000000L), valueOf(0x6765c793fa10079dL),
        valueOf(0x41c21cb8e1000000L), valueOf(0x3642798750226111L),
        valueOf(0x1000000000000000L), valueOf(0x12bf307ae81ffd59L),
        valueOf( 0xde0b6b3a7640000L), valueOf(0x4d28cb56c33fa539L),
        valueOf(0x1eca170c00000000L), valueOf(0x780c7372621bd74dL),
        valueOf(0x1e39a5057d810000L), valueOf(0x5b27ac993df97701L),
        valueOf(0x1000000000000000L), valueOf(0x27b95e997e21d9f1L),
        valueOf(0x5da0e1e53c5c8000L), valueOf( 0xb16a458ef403f19L),
        valueOf(0x16bcc41e90000000L), valueOf(0x2d04b7fdd9c0ef49L),
        valueOf(0x5658597bcaa24000L), valueOf( 0x6feb266931a75b7L),
        valueOf( 0xc29e98000000000L), valueOf(0x14adf4b7320334b9L),
        valueOf(0x226ed36478bfa000L), valueOf(0x383d9170b85ff80bL),
        valueOf(0x5a3c23e39c000000L), valueOf( 0x4e900abb53e6b71L),
        valueOf( 0x7600ec618141000L), valueOf( 0xaee5720ee830681L),
        valueOf(0x1000000000000000L), valueOf(0x172588ad4f5f0981L),
        valueOf(0x211e44f7d02c1000L), valueOf(0x2ee56725f06e5c71L),
        valueOf(0x41c21cb8e1000000L)};

    /*
     * These two arrays are the integer analogue of above.
     */
    private static int digitsPerInt[] = {0, 0, 30, 19, 15, 13, 11,
        11, 10, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5};

    private static int intRadix[] = {0, 0,
        0x40000000, 0x4546b3db, 0x40000000, 0x48c27395, 0x159fd800,
        0x75db9c97, 0x40000000, 0x17179149, 0x3b9aca00, 0xcc6db61,
        0x19a10000, 0x309f1021, 0x57f6c100, 0xa2f1b6f,  0x10000000,
        0x18754571, 0x247dbc80, 0x3547667b, 0x4c4b4000, 0x6b5a6e1d,
        0x6c20a40,  0x8d2d931,  0xb640000,  0xe8d4a51,  0x1269ae40,
        0x17179149, 0x1cb91000, 0x23744899, 0x2b73a840, 0x34e63b41,
        0x40000000, 0x4cfa3cc1, 0x5c13d840, 0x6d91b519, 0x39aa400
    };

    /**
     * These routines provide access to the two's complement representation
     * of OpenBigIntegers.
     */

    /**
     * Returns the length of the two's complement representation in ints,
     * including space for at least one sign bit.
     */
    private int intLength() {
        return (bitLength() >>> 5) + 1;
    }

    /* Returns sign bit */
    private int signBit() {
        return signum < 0 ? 1 : 0;
    }

    /* Returns an int of sign bits */
    private int signInt() {
        return signum < 0 ? -1 : 0;
    }

    /**
     * Returns the specified int of the little-endian two's complement
     * representation (int 0 is the least significant).  The int number can
     * be arbitrarily high (values are logically preceded by infinitely many
     * sign ints).
     */
    private int getInt(int n) {
        if (n < 0)
            return 0;
        if (n >= mag.length)
            return signInt();

        int magInt = mag[mag.length-n-1];

        return (signum >= 0 ? magInt :
                (n <= firstNonzeroIntNum() ? -magInt : ~magInt));
    }

    /**
     * Returns the index of the int that contains the first nonzero int in the
     * little-endian binary representation of the magnitude (int 0 is the
     * least significant). If the magnitude is zero, return value is undefined.
     */
     private int firstNonzeroIntNum() {
         int fn = firstNonzeroIntNum - 2;
         if (fn == -2) { // firstNonzeroIntNum not initialized yet
             fn = 0;

             // Search for the first nonzero int
             int i;
             int mlen = mag.length;
             for (i = mlen - 1; i >= 0 && mag[i] == 0; i--)
                 ;
             fn = mlen - i - 1;
             firstNonzeroIntNum = fn + 2; // offset by two to initialize
         }
         return fn;
     }

    /** use serialVersionUID from JDK 1.1. for interoperability */
    private static final long serialVersionUID = -8287574255936472291L;

	private OpenBigInteger2 multiplyer;
	private int[] multiMag = null;

	public BigInteger toBigInteger() {
		return new BigInteger(this.toByteArray());
	}

	public int getLastDigit() {
		int lastInt = getInt(0);
		return lastInt & 0x0000FFFF;
	}

	int start = 0;
	
	public void divideByBase() {
		int last16bytes = 0 ,prevBytes = 0;
		for (int n=start; n<mag.length; n++){
			int magInt = mag[n];
			
			// store the last 16 bytes
			prevBytes = last16bytes;
			last16bytes = (magInt & 0x0000FFFF) << 16;
			mag[n] >>= 16;
			mag[n] &=  0x0000FFFF;
		
			if (n>start) {
				mag[n] |= prevBytes;
			}
		}
		
		// If the most significant byte is zero, then we should ignore the first int
		// by incrementing start
		if (mag[start] == 0 && (mag.length - start) > 1){
			start++;
			this.bitLength = 0;
		}
	}

	int[] lastDigitHolder = new int[1];
	
	// working now
	public void subtractMultipication(OpenBigInteger2 multiplyer, int lastDigit) {
		lastDigitHolder[0] = lastDigit;
		multiplyToLen(multiplyer.mag, multiplyer.mag.length, this.start, lastDigitHolder, 1, multiMag);
		
		int keep;
		for (keep = 0; keep < multiMag.length && multiMag[keep] == 0; keep++)
            ;
		int multimagStart = keep;
		// old code :
//		int[] multiMagToUse;
//		if (multiMag[0] == 0)
//			multiMagToUse = java.util.Arrays.copyOfRange(multiMag, 1, multiMag.length);
//		else
//			multiMagToUse = multiMag;
		
		int cmp = compareMagnitude(multiMag, multimagStart);
		//System.out.println("multiMagToUse:" + Arrays.toString(multiMag) + " multimagStart:"+multimagStart+ " cmp:" + cmp);
        if (cmp == 0){
        	this.mag = new int[0];
        	this.signum = 0;
        	start = 0;
        } else {
	        this.mag = (cmp > 0 ? subtract(mag, start,  multiMag , multimagStart)
	                           : subtract(multiMag, multimagStart, mag, start));
//	        resultMag = trustedStripLeadingZeroInts(resultMag);

	        // Find first nonzero byte
	        for (keep = 0; keep < mag.length && mag[keep] == 0; keep++)
	            ;
	        start = keep;
	        this.signum = (this.mag.length -start) == 0 ? 0 : (cmp == signum ? 1 : -1);
        }
        this.bitLength = 0;
        this.bitCount = 0;
	}
	
	final int compareMagnitude(int[] m2, int m2Ignore) {
        int[] m1 = mag;
        int len1 = m1.length - start;
        int len2 = m2.length - m2Ignore;
        if (len1 < len2)
            return -1;
        if (len1 > len2)
            return 1;
        for (int i = 0; i < len1; i++) {
            int a = m1[start+i];
            int b = m2[m2Ignore+i];
            if (a != b)
                return ((a & LONG_MASK) < (b & LONG_MASK)) ? -1 : 1;
        }
        return 0;
    }

	public void initializeMultiplyer(OpenBigInteger2 multiplyer) {
		this.multiplyer = multiplyer;
		this.multiMag = new int[multiplyer.mag.length + 1];
	}

	public void step() {
		// Get the last digit :
		int lastDigit = this.getLastDigit();
		
		// Divide by the base :
		divideByBase();
//		System.out.println("step :" + this.toBigInteger().toString());
//		System.out.println("Before substraction :" + Arrays.toString(mag) +" length :" + mag.length +" start :" + start);
		subtractMultipication(multiplyer, lastDigit);
//		System.out.println("After substraction :" + Arrays.toString(mag) +" length :" + mag.length +" start :" + start);
	}

}
