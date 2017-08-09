#!/usr/bin/python

'''Run this python script to generate values for the uv_txsize_lookup array.

This array stores four transform sizes for each block size and starting
transform size. These transform sizes are the transform size to use when
partitioning the block: not at all, horizontally, vertically, both ways.

For example, if the block size is BLOCK_8X8 and the transform size is TX_8X4,
the list of sizes is { TX_8X4, TX_8X4 }, { TX_4X4, TX_4X4 }. This is because
the partitioned block sizes are BLOCK_8X8, BLOCK_4X8, BLOCK_4X8, BLOCK_4X4 and
each transform listed is the largest transform that fits in both TX_8X4 and its
corresponding partitioned block size.

The output has #if / #endif pairs to deal with configuration options. The
simplest question is whether a given block size exists (BLOCK_2X2 through
BLOCK_4X2 require CONFIG_CHROMA_2X2 or CONFIG_CHROMA_SUB8X8; BLOCK_64X128
through BLOCK_128X128 require CONFIG_EXT_PARTITION). More complicated, the
correct answer for an entry might depend on the set of config options. For
example, see BLOCK_4X4: if CONFIG_CHROMA_2X2 is set, the blocks will divide
down to 2X2, 2X4 or 4X2, giving TX_2X2 in the output. If not, the transform
sizes will never get smaller than 4X4.

At the moment, there is a constraint that a transform never gets less square
than it started (so a square never gets chopped into a 2:1 rectangle, for
example). There is also a constraint about when 2X2 transforms get generated,
which tries to match the previous handwritten table. These constraints are
guaranteed by the code in greatest_txforms_within, which can easily be changed
if we decide to allow such divisions.

'''

from __future__ import print_function
import sys

# The block sizes from enums.h
__BLOCK_SIZES = [
    (2, 2),
    (2, 4),
    (4, 2),
    (4, 4),
    (4, 8),
    (8, 4),
    (8, 8),
    (8, 16),
    (16, 8),
    (16, 16),
    (16, 32),
    (32, 16),
    (32, 32),
    (32, 64),
    (64, 32),
    (64, 64),
    (64, 128),
    (128, 64),
    (128, 128),
    (4, 16),
    (16, 4),
    (8, 32),
    (32, 8),
    (16, 64),
    (64, 16)
]

# Special "OR" config options. If an option is listed as a key here, it's an
# alias that should be written out as the the OR of its values.
__CONFIG_ORS = {
    'CHROMA_2X2_OR_SUB8X8': {'CHROMA_2X2', 'CHROMA_SUB8X8'}
}

# A list of which block sizes are conditional on which flags. A key of X means
# we need to hide that block size behind #if CONFIG_X
__BLOCK_SIZE_SWITCHES = {
    'EXT_PARTITION': {(64, 128), (128, 64), (128, 128)},
    'CHROMA_2X2_OR_SUB8X8': {(2, 2), (2, 4), (4, 2)},
}

# The transform sizes from enums.h
__TX_SIZES = [
    (2, 2),
    (4, 4),
    (8, 8),
    (16, 16),
    (32, 32),
    (64, 64),
    (4, 8),
    (8, 4),
    (8, 16),
    (16, 8),
    (16, 32),
    (32, 16),
    (4, 16),
    (16, 4),
    (8, 32),
    (32, 8)
]

# A list of which transform sizes are conditional on which flags (the same
# format as __BLOCK_SIZE_SWITCHES).
__TX_SIZE_SWITCHES = {
    'CHROMA_2X2': {(2, 2)},
    'TX64X64': {(64, 64)}
}

def reverse_map_to_set(dictionary):
    '''Reverse a dictionary that was mapping keys to sets of values.

    The new map will take values to sets of keys'''
    ret = {}
    for (k, vals) in dictionary.items():
        for val in vals:
            if val not in ret:
                ret[val] = {k}
            else:
                ret[val].add(k)
    return ret

# The reverse maps for __BLOCK_SIZE_SWITCHES and __TX_SIZE_SWITCHES. These are
# what we need for looking flags up, but it is easier to write them down in the
# other direction. __BLOCK_SIZE_REVERSE is keyed by block size and (if there is
# a value) the value is a set of switches behind which that block size
# hides. __TX_SIZE_REVERSE is similar, but keyed by transform size.
__BLOCK_SIZE_REVERSE = reverse_map_to_set(__BLOCK_SIZE_SWITCHES)
__TX_SIZE_REVERSE = reverse_map_to_set(__TX_SIZE_SWITCHES)

def size_to_name(size_type, size):
    '''Turn a size (pair of numbers) into a name (a string) with size_type'''
    assert len(size) == 2
    return '{}_{}X{}'.format(size_type, size[0], size[1])

def if_single_condition(flag, is_true):
    '''Format the condition to pass to #if for flag being true/false.'''
    if flag in __CONFIG_ORS:
        parts = sorted(__CONFIG_ORS[flag])
        joined = ' || '.join(['CONFIG_' + p for p in parts])
        if is_true:
            return joined
        else:
            return '!({})'.format(joined)
    else:
        if is_true:
            return 'CONFIG_' + flag
        else:
            return '!CONFIG_' + flag


def if_condition(yes_flags, no_flags=None):
    '''Format a condition suitable for #if from yes_flags/no_flags'''
    assert yes_flags or no_flags
    cfg_options = [if_single_condition(f, True) for f in yes_flags]
    if no_flags:
        for no_flag in no_flags:
            cfg_options.append(if_single_condition(no_flag, False))

    return ' && '.join(cfg_options)


def print_if_line(stream, yes_flags, no_flags=None, extra_comment=None):
    '''Write the #if or #endif line for the given flags to stream'''
    assert yes_flags or no_flags
    stream.write('#if {}'.format(if_condition(yes_flags, no_flags)))
    if extra_comment:
        stream.write('  // {}'.format(extra_comment))
    stream.write('\n')


def format_tx_line(tx_sizes):
    '''Generate a line of TX_xxx sizes suitable for the table'''
    assert len(tx_sizes) == 4
    return ('      {{ {{ {}, {} }}, {{ {}, {} }} }},'
            .format(size_to_name('TX', tx_sizes[0]),
                    size_to_name('TX', tx_sizes[1]),
                    size_to_name('TX', tx_sizes[2]),
                    size_to_name('TX', tx_sizes[3])))

def fits_in(size_a, size_b):
    '''True if size_a fits inside size_b'''
    return size_a[0] <= size_b[0] and size_a[1] <= size_b[1]

def greatest_txforms_within(orig_bsize, bsize, txform, base_flags):
    '''Find the largest transform less that fits in both txform and bsize.

    Also don't make the result less square than the original transform unless
    the resulting transform is the same shape as the original block size.

    If bsize is too small for anything to fit (1X2 or 2X1), use 2X2 if
    CONFIG_CHROMA_2X2 is true, else 4X4.

    Even if CONFIG_CHROMA_2X2 is true, don't produce an output transform
    smaller than 4X4 unless either the original block size was strictly smaller
    than 4x4 or the original block size fits in 8x8 and the transform is <=4
    wherever the original block size is.

    Since the answer might depend on config options, we return a list of pairs,
    (flags, size), where earlier entries are more preferred. flags is the set
    of config flags that must be true for this to be the right answer. size is
    the txform size. The last entry in the list will always be (set(), size)
    for some size.

    '''
    best = None
    best_area = 0
    def thinness(x):
        return max(x) / min(x)

    thinness_allowed = thinness(txform)

    for tx_size in __TX_SIZES:
        if not fits_in(tx_size, txform):
            continue
        if not fits_in(tx_size, bsize):
            continue

        tx_thinness = thinness(tx_size)
        if tx_thinness > thinness_allowed:
            if thinness_allowed == 2:
                continue
            if tx_size[0]*orig_bsize[1] != orig_bsize[0]*tx_size[1]:
                continue

        area = tx_size[0] * tx_size[1]
        if best and area <= best_area:
            continue

        best = tx_size
        best_area = area

    # If we didn't find anything that fits, choose 2x2.
    if not best:
        best = (2, 2)

    # Only allow 2x2 if either:
    #
    #   1. The original block size is strictly smaller than 4x4
    #
    #   2. The original block size is <= 8x8 and the transform is <= 4 wherever
    #      the original block size <= 4
    if min(best) < 4:
        case1 = max(orig_bsize) <= 4 and min(orig_bsize) < 4
        case2 = ((max(orig_bsize) <= 8) and
                 (orig_bsize[0] > 4 or txform[0] <= 4) and
                 (orig_bsize[1] > 4 or txform[1] <= 4))

        if not (case1 or case2):
            best = (4, 4)
        else:
            best = (2, 2)

    # If min(best) < 4, this is only valid when CONFIG_CHROMA_2X2, so we have
    # to return two values, unless CHROMA_2X2 is in base_flags, in which case
    # we can return best unconditionally
    if min(best) < 4 and 'CHROMA_2X2' not in base_flags:
        return [({'CHROMA_2X2'}, best), (set(), (4, 4))]

    return [(set(), best)]


def txsize_options(block_size, tx_size, base_flags):
    '''Merge the results of greatest_txforms_within to find the possible values
    for the array entry.

    The result is a pair (FS, TXS) where TXS has length 4 and TX[i] is the
    result of calling greatest_txforms_within for the i'th divided transform
    size and FS is the set of all flags that appear as a key in at least one of
    TX[0]..TX[3] but that are not in base_flags.

    '''
    bsizes = [(block_size[0], block_size[1]),
              (block_size[0], block_size[1] / 2),
              (block_size[0] / 2, block_size[1]),
              (block_size[0] / 2, block_size[1] / 2)]

    quad = [greatest_txforms_within(block_size, bs, tx_size, base_flags)
            for bs in bsizes]

    flags = set()
    for lst in quad:
        for (flagset, _) in lst:
            flags |= flagset
    flags -= base_flags

    return (flags, quad)


def flags_for_index(flagset, idx, base_flags=None):
    '''Return base_flags plus the set of flags from flagset selected by idx in
    binary

    In order to show #if P before #if ! P, we reverse the standard meaning, so
    a flag is true if the corresponding bit is not set.

    '''
    ret = base_flags.copy() if base_flags else set()
    num_bits = len(flagset)
    for j in range(num_bits):
        if (idx >> j) & 1 == 0:
            ret.add(flagset[j])
    return ret


def tx_for_flags(txforms, flags):
    '''Return the first match for flags in txforms (which should be of the
    format returned from greatest_txforms_within).

    '''
    for (cflags, txform) in txforms:
        if cflags <= flags:
            return txform
    # We should never get here because lists of choices from
    # greatest_txforms_within should all end with a choice with no flags.
    assert False


def merge_txsize_options(opts, base_flags):
    '''Merge a list of results from txsize_options (corresponding to successive
    entries in the output array).

    The result is a pair (FS, ARRAYS) where FS is an ordered list of flags and
    ARRAYS has length 2**len(FS). Entry ARRAYS[i] is the set of rows to use
    with the set of flags where F[j] holds if and only if (i >> j) & 1 != 0.

    '''
    all_flags = set()
    for (flags, _) in opts:
        all_flags |= flags

    ordered_flags = list(all_flags)

    arrays = []
    for i in range(2**len(ordered_flags)):
        flags_for_i = flags_for_index(ordered_flags, i, base_flags)
        rows = []
        for (_, quad) in opts:
            rows.append([tx_for_flags(txforms, flags_for_i) for txforms in quad])
        arrays.append(rows)
    return (ordered_flags, arrays)


def print_block_size_lines(stream, block_size, bs_flags):
    '''Print the lines of TX choices for a given block_size'''

    # Firstly, get hold of the options for each input transform. base_flags is
    # the flags that have to hold for the transform to be written out (so we
    # can pass them as base flags to txsize_options)
    options = []
    for tx_size in __TX_SIZES:
        tx_flags = __TX_SIZE_REVERSE.get(tx_size, set())
        base_flags = tx_flags | bs_flags
        options.append(txsize_options(block_size, tx_size, base_flags))

    # Now merge these together. This time the set of base flags is just
    # bs_flags. Note that if a result for a given tx_size depends on a flag
    # that was set in tx_flags on that row, that flag won't appear in the flags
    # (options[i][0] for some i) so we won't get a spurious #if.
    opt_flags, arrays = merge_txsize_options(options, bs_flags)

    # Finally, write out the different possibilities.
    for i, rows in enumerate(arrays):
        case_flags = set()
        case_no_flags = set()
        if opt_flags:
            case_flags = flags_for_index(opt_flags, i)
            case_no_flags = set(opt_flags) - case_flags
            print_if_line(stream, case_flags, case_no_flags,
                          extra_comment=('case {}/{}'
                                         .format(i + 1, len(arrays))))

        for tx_size, tx_row in zip(__TX_SIZES, rows):
            # extra_tx_flags is the set of flags that must hold in order for
            # this input transformation size to exist. Since we're at the
            # current block size and in the current case, we can subtract any
            # flags that apply from that.
            tx_flags = __TX_SIZE_REVERSE.get(tx_size, set())
            extra_tx_flags = tx_flags - bs_flags - case_flags

            # If extra_tx_flags has intersection with case_no_flags, we know
            # this case can't happen.
            if extra_tx_flags & case_no_flags:
                continue

            if extra_tx_flags:
                print_if_line(stream, extra_tx_flags)

            stream.write(format_tx_line(tx_row))
            stream.write('\n')

            if extra_tx_flags:
                stream.write('#endif\n')

        if opt_flags:
            stream.write('#endif\n')


def main():
    '''Main entry point'''
    print('// An array of transform sizes to use when partitioning blocks.\n'
          '//\n'
          '// Entries are keyed by block and initial transform size and give\n'
          '// the transform sizes to use with: no split, a horizontal\n'
          '// split, a vertical split and a split both ways.\n'
          '//\n'
          '// This table was generated by tools/gen_uv_txsize_lookup.py\n'
          'static const TX_SIZE '
          'uv_txsize_lookup[BLOCK_SIZES_ALL][TX_SIZES_ALL][2][2] = {')
    last_bs_flags = None
    for bsize in __BLOCK_SIZES:
        bs_flags = __BLOCK_SIZE_REVERSE.get(bsize, set())
        if bs_flags != last_bs_flags:
            if last_bs_flags:
                print('#endif')
            if bs_flags:
                print_if_line(sys.stdout, bs_flags)

        last_bs_flags = bs_flags
        print('  {{\n'
              '// {}'
              .format(size_to_name('BLOCK', bsize)))

        print_block_size_lines(sys.stdout, bsize, bs_flags)

        print('  },')

    if last_bs_flags:
        print('#endif')

    print('};')


if __name__ == "__main__":
    main()
