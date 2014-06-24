from numpy import *

# Prices based on gen9 "summer sale", valid until Jul 7 2014:
# http://www.gen9bio.com/why-gen9/price/summer-sale/
# = 20,000 bp minimum order =
# 500 bp to 1000 bp: $149
# 1001 bp to 2000 bp: $249
# 2001 bp to 3000 bp: $349
# = 200,000 bp minimum order =
# 500 bp to 1000 bp: $99
# 1001 bp to 2000 bp: $199
# 2001 bp to 3000 bp: $299

def construct_cost(construct_length, total_nbps):
    if construct_length > 3000 or construct_length < 500:
        raise Exception, 'Construct length (%d bps) outside of bounds supported by gen9 (500-3000 bps).' % construct_length

    if total_nbps < 20000:
        raise Exception, 'Total bps (%d) less than gen9 minimum order size of 20,000 bp.' % total_nbps

    if total_nbps >= 20000 and total_nbps < 200000:
        if construct_length <= 1000:
            cost = 149
        elif construct_length <= 2000:
            cost = 249
        elif construct_length <= 3000:
            cost = 349
    elif total_nbps > 200000:
        if construct_length <= 1000:
            cost = 99
        elif construct_length <= 2000:
            cost = 199
        elif construct_length <= 3000:
            cost = 299
    return cost

# copied/pasted data from "Gen9_Seq_Submission_Short Form-choderalab-41_sgc_kinases.xlsx"
sgc_constructs_bp_lens = array([843,
1002,
909,
879,
879,
924,
1077,
654,
870,
897,
900,
906,
903,
1008,
894,
1086,
1125,
1230,
813,
1011,
864,
987,
927,
879,
870,
870,
966,
987,
864,
1074,
885,
906,
915,
978,
933,
984,
969,
984,
864,
1185,
1272])

total_bps = sum(sgc_constructs_bp_lens)

costs = array(map(construct_cost, sgc_constructs_bp_lens, [total_bps] * len(sgc_constructs_bp_lens)))
total_price = sum(costs)

print 'Total number of constructs: %d' % len(sgc_constructs_bp_lens)
print 'Total number of bps: %d' % sum(sgc_constructs_bp_lens)

for unique_cost in sorted(list(set(costs))):
    nconstructs = sum(costs == unique_cost)
    print 'Total constructs at cost $%d: %d' % (unique_cost, nconstructs)

print 'Total price: $%d' % total_price

