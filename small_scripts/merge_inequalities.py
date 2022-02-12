""" Merging two inequality files """

import argparse

# parse the inputs
parser = argparse.ArgumentParser()
parser.add_argument(dest='input_file_1', type=str, help='First input file.')
parser.add_argument(dest='input_file_2', type=str, help='Second input file.')
parser.add_argument(dest='output_file', type=str, help='Output file.')
args = parser.parse_args()

# load files
with open(args.input_file_1) as f:
    f1 = f.readlines()
with open(args.input_file_2) as f:
    f2 = f.readlines()

# Find indices
start1 = 0
for i in range(len(f1)):
    if 'Inequalities' in f1[i]:
        start1 = i
        break
start2 = 0
for i in range(len(f2)):
    if 'Inequalities' in f2[i]:
        start2 = i
        break


# generate new file
res  = f1[start1:] + f2[start2 + 1:]

# save new file
outfile = open(args.output_file, 'w+')
outfile.writelines(res)
outfile.close()
print('Wrote result to file: {}'.format(args.output_file))

