# Scrape the OEIS and extract potential Gauss sequences

# Make sure you have the files 'names' and 'stripped' unzipped and in this folder
# Download these from https://oeis.org/wiki/Compressed_Versions_of_the_Database

# How many entries of the OEIS to read, set -1 for infinite
MAX_ITER = 10000

# How many values to check the congruence for
CHECK = 20
# Largest value allowed in the checked terms
MAX_VALUE = 10**8
# Whether we require the Lyndon colours to be nonnegative
REQUIRE_POS_COLOURS = True
# Whether we require the Lyndon parmaters to be nonnegative
REQUIRE_POS_PARAMS = False
# Whether we ignore sequences that are mostly the same value
IGNORE_REPETITION = True
# Shift values that are checked (must be nonnegative)
SHIFT = [0,1,2]
# Write file name
NAME = "gauss-seqs.txt"


# Mobius-mu values and divisors for n = 0,1,...,12
DIV = [[0]] * (CHECK + 1)
MU = [0] * (CHECK + 1)
MU[1] = -1
for n in range(1, CHECK + 1):
    div_n = [d for d in range(1,n+1) if n % d == 0]
    DIV[n] = div_n
    MU[n] = -sum([MU[d] for d in div_n])


# Read and write files
read_nums = open('stripped', 'r', encoding = "utf8")
read_names = open('names', 'r', encoding = "utf8")
w = open(NAME, 'w', encoding = "utf8")


# Check if the sequence [a_1,...,a_n] satisfies the Gauss congruence
def is_gauss_seq(vals):
    if len(vals) > 1+min(len(MU), len(DIV)):
        return None
    for n in range(1,len(vals)+1):
        if not (sum([MU[round(n/d)]*vals[d-1] for d in DIV[n]]) % n == 0):
            return False
    return True

# Find the colour sequence for a given Gauss-congruent sequence [a_1,...,a_n]
def lyndon_colours(vals):
    cols = [0] * len(vals)
    for n in range(1,len(vals)+1):
        cols[n-1] = (vals[n-1] - sum([vals[n-k-1]*cols[k-1] for k in range(1,n)]))/n
    return cols

# Find the Lyndon parameters for a given Gauss-congruent sequence [a_1,...,a_n]
def lyndon_params(vals):
    return [round(sum([MU[round(n/d)]*vals[d-1] for d in DIV[n]])/n) for n in range(1,len(vals)+1)]


# Check if the entries [a_1,...,a_n] satisfy the required properties
def is_valid_entry(vals):
    if len(vals) < CHECK:
        return False
    if max(vals) > MAX_VALUE:
        return False
    if not is_gauss_seq(vals):
        return False
    cols = lyndon_colours(vals)
    if REQUIRE_POS_COLOURS and sum([c >= 0 for c in cols]) < len(vals):
        return False
    params = lyndon_params(vals)
    if REQUIRE_POS_PARAMS and sum([p >= 0 for p in params]) < len(vals):
        return False
    if IGNORE_REPETITION and sum([v == vals[-1] for v in vals]) > len(vals)/2-1:
        return False
    # If all checks are passed, return true
    return True


iter = 0
while iter < MAX_ITER + 4 or MAX_ITER < 0:
    nums = read_nums.readline()
    if len(nums) == 0:
        break
    name = read_names.readline().strip().split(' ')
    if nums[0] == 'A':
        entries = nums.split(',')
        # Check there are enough entries, and they aren't too big
        for s in SHIFT:
            if len(entries) > s+2:
                vals = [int(x) for x in entries[s+1:-1]]
                if is_valid_entry(vals[:CHECK]):
                    str_name = name[0] + "[{}]".format(s)
                    # Pad the string description, then 
                    str_desc_long = " ".join(name[1:]) + " "*100
                    str_desc = str_desc_long[:100]
                    str_vals = ",".join([str(v) for v in vals[1:CHECK]])
                    print(str_name)
                    w.write(" ".join([str_name, str_desc, str_vals]) + "\n")
    iter += 1

w.close()
read_names.close()
read_nums.close()