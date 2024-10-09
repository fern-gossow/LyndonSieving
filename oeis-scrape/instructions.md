# Instructions for OEIS-scrape

Here we provide instructions for how to use `oeis-scrape.py` to find examples of sequences satisfying Gauss congruence. Ensure that you have [python3](https://www.python.org/) installed. 

The program reads the compressed forms of the OEIS database, which can be [downloaded here](https://oeis.org/wiki/Compressed_Versions_of_the_Database). Ensure you have these files (unzipped) in the same folder as the program, which will be called `names` and `stripped`. Running the program will then create a new file consisting of the sequence numbers and shifts, the name of the sequence, and the first few values.

The code has a number of options:
- `MAX_ITER` gives the number of lines in the OEIS to iterate through. Set to -1 if you want to read the entire database.
- `CHECK` is the number of terms that Gauss congruence is checked for. Setting this number too low will give additional terms due to random change, but setting this value too high will increase runtime and remove sequences that do not have enough values present.
- `MAX_VALUE` gives the maximum value in the sequence of checked values to allow for.
- `REQUIRE_POS_COLOURS` determines whether the colour values need to be nonnegative.
- `REQUIRE_POS_PARAMS` determines whether the Lyndon parameters need to be nonnegative.
- `IGNORE_REPETITION` determines whether to ignore sequences that contain many repeated values, such as the sequence of zeros.
- `SHIFT` gives the starting indexes that Gauss congruence will be checked at.
- `NAME` gives the filename to write the sequences to.

A sample file is provided with the first 100000 entries of the OEIS, requiring positive colours and ignoring repetition.