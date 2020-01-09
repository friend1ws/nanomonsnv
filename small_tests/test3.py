#! /usr/bin/env python

def my_sum(output_file, *args):

    hout = open(output_file, 'w')
    print(sum(args), file = hout)
    hout.close()


my_sum("sum.test.txt", 1, 3, 4)

# import argparse


# parser = argparse.ArgumentParser()

# parser.add_argument("output_file", type = str, help = "the path to output file")
# parser.add_argument("values", type = int, nargs = "*", help = "the value to be summed")


