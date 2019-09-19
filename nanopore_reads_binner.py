#!/usr/local/env python3.6

import os
import gzip
import sys
from argparse import ArgumentParser
from time import time
from collections import defaultdict
from dateutil.parser import parse
import pathlib
import operator
import numpy as np
from math import ceil


class FastqObjects(object):
    def __init__(self, header, seq, qual, time_string, bin):
        # Create seq object with its attributes
        self.header = header
        self.seq = seq
        self.qual = qual
        self.time_string = time_string
        self.bin = bin


class NanoReadsBinner(object):
    def __init__(self, args):
        self.args = args
        self.input_folder = args.fastq
        self.output_folder = args.output
        self.interval = args.interval
        self.prefix = args.prefix
        self.cpu = args.threads

        # Shared data structure(s)
        self.sample_dict = defaultdict()

        # Create a list of fastq files in input folder
        self.input_fastq_list = list()

        # Time tracking
        self.total_time = list()

        # Unit used
        self.bin_size = None
        self.units = None

        # run the script
        self.run()

    def run(self):
        """
        Run everything
        :return:
        """
        # Initial time stamp
        self.total_time.append(time())

        # Check if correct argument combination is being used
        self.check_args()

        # Parse fastq files into dictionary
        if self.input_folder:
            print("Parsing fastq files...", end="", flush=True)
            start_time = time()
            read_counter = NanoReadsBinner.parse_fastq_list_parallel(self.input_fastq_list, self.sample_dict, self.cpu)
            end_time = time()
            interval = end_time - start_time
            print(" took %s for %d reads" % (NanoReadsBinner.elapsed_time(interval), read_counter))

            # Check if there is data
            if not self.sample_dict:
                raise Exception('No data!')
            else:
                self.convert_datetimes(self.sample_dict)  # convert datetimes to elapsed time
                self.assign_bins(self.sample_dict)

        self.total_time.append(time())
        print("\n Total run time: {}".format(NanoReadsBinner.elapsed_time(self.total_time[1] - self.total_time[0])))

    def check_args(self):
        """
        Check if correct argument combination is being used
        Check if output folder exists
        :return:
        """

        # Check output folder
        if not self.output_folder:
            print('Please specify an output folder ("-o")')
            parser.print_help(sys.stderr)
            sys.exit(1)
        else:
            pathlib.Path(self.output_folder).mkdir(parents=True, exist_ok=True)  # Create if if it does not exist

        # Check input folder
        if self.input_folder:
            for root, directories, filenames in os.walk(self.input_folder):
                for filename in filenames:
                    absolute_path = os.path.join(root, filename)
                    if os.path.isfile(absolute_path) and filename.endswith(('.fastq', '.fastq.gz', 'fq', 'fq.gz')):
                        self.input_fastq_list.append(absolute_path)

            # check if input_fastq_list is not empty
            if not self.input_fastq_list:
                raise Exception("No fastq files found in {}!".format(self.input_folder))

        # Check interval
        good_units_list = ('h', 'm', 's')
        self.bin_size = float(self.interval[:-1])
        self.units = self.interval[-1]
        if self.units not in good_units_list:
            raise Exception("Please use one of the following unit for the interval: {}".format(
                ', '.join(good_units_list)))
        # TODO -> check if interval size is an integer

    @staticmethod
    def elapsed_time(seconds):
        """
        Transform a time value into a string
        :param seconds: Elapsed time in seconds
        :return: Formatted time string
        """
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
        time_string = ''.join('{}{}'.format(int(np.round(value)), name) for name, value in periods if value)

        return time_string

    @staticmethod
    def parse_fastq_list_parallel(l, d, num_threads):
        """
        Parse fastq files in parallel
        :param l: list of files
        :param d: dictionary
        :param num_threads: number of files to process in parallel
        :return: Number or reads
        """
        # with futures.ProcessPoolExecutor(max_workers=num_threads) as pool:
        #     results = pool.map(NanoReadsBinner.parse_file, l, chunksize=1)

        results = list()
        for f in l:
            result = NanoReadsBinner.parse_file(f)
            results.append(result)

        # Update dictionary with results from every chunk
        read_counter = 0
        for dictionary in results:
            read_counter += len(dictionary.keys())
            d.update(dictionary)  # Do the merge

        return read_counter

    @staticmethod
    def parse_fastq_to_dict(l, my_dict):
        """
        Get info and stats from header, sequence and quality lines . Update master dictionary with info.
        :param l: a list of 4 items, for each line of the fastq entry
        :param my_dict: an empty dictionary to store results
        :return:
        """
        # Put list items (fastq lines) into variables
        header, seq, extra, qual = l  # get each component of list in a variable

        # Split header at space character into a list
        items = header.split()

        # Sequence ID
        seq_id = items[0]

        # Read Time stamp
        time_string = next((t for t in items if b'start_time=' in t), None)
        time_string = time_string.split(b'=')[1]
        time_string = parse(time_string)

        # Put relevant info into object
        seq_object = FastqObjects(header, seq, qual, time_string, None)

        # Add to dictionary
        my_dict[seq_id] = seq_object

    @staticmethod
    def parse_file(f):
        """
        Open file and read 4 lines (one fastq entry) and call parse_fastq_to_dict to compute metrics
        :param f: Fastq file to parse
        :return: dictionary
        """

        # Check file size, skip empty fastq
        size = os.path.getsize(f)
        if size == 0:
            return  # Exit function and don't process that file

        # Parse
        my_dict = {}
        with gzip.open(f, 'rb', 1024 * 1024) if f.endswith('gz') else open(f, 'rb', 1024 * 1024) as file_handle:
            lines = []  # a list to store the 4 lines of a fastq entry in order
            for line in file_handle:
                if not line:  # end of file?
                    break
                line = line.rstrip()
                if len(lines) == 4:  # full entry. Assume that the file is OK
                    NanoReadsBinner.parse_fastq_to_dict(lines, my_dict)
                    lines = []
                lines.append(line)
            # Parse the last entry of the file
            if len(lines) == 4:
                NanoReadsBinner.parse_fastq_to_dict(lines, my_dict)

        return my_dict

    def convert_datetimes(self, d):
        """

        :param d:
        :return:
        """
        # Find the smallest datetime
        t_min = min(y.time_string for x, y in d.items())

        for read, info in d.items():
            elapsed = info.time_string - t_min
            # Convert to unit requested (h, m or s)
            if self.units == 'h':
                elapsed = elapsed.days * 24 + elapsed.seconds / 3600  # Convert to hours (float)
            elif self.units == 'm':
                elapsed = elapsed.days * 1440 + elapsed.seconds / 60  # Convert to minutes (float)
            else:  # 's'
                elapsed = elapsed.days * 86400 + elapsed.seconds  # Convert to hours (float)
            elapsed = round(elapsed, 2)
            info.time_string = elapsed

    def assign_bins(self, d):
        """
        Reorganize the reads based on sequencing time intervals
        :param d: dictionary with all the reads of a sequencing run
        :return:
        """
        # Find the smallest and biggest datetime
        t_max = max(y.time_string for x, y in d.items())

        # Bin data
        # https://stackoverflow.com/questions/6163334/binning-data-in-python-with-scipy-numpy
        num_bins = ceil(t_max / self.bin_size) + 1
        t_max = num_bins * self.bin_size
        bins = np.linspace(0, t_max, num=num_bins, endpoint=False)

        # Sort dictionary the object attribute value
        for read in sorted(d.values(), key=operator.attrgetter('time_string')):
            # Assign bin read
            read.bin = np.digitize(read.time_string, bins)

        # TODO -> parallel output all the bins instead of having to go through the
        #         whole dictionary for each bin
        # TODO -> check for any empty file and remove them

        # TODO -> no need to do the last bin, just copy and rename the original file

        counter = 0
        for i in bins[1:]:  # skip first one, which is zero
            print('\tBinning reads from 0 to {}{}... '.format(int(i), self.units), end="", flush=True)
            with gzip.open(self.output_folder + '/' + self.prefix
                           + '_0-' + str(int(i)) + self.units + '.fastq.gz', 'wt') as out_fh:
                for ident, info in d.items():
                    if info.time_string < i:
                        counter = counter + 1
                        out_fh.write("{}\n{}\n{}\n{}\n".format(info.header, info.seq, '+', info.qual))
            print('{} reads'.format(counter))
            counter = 0


if __name__ == '__main__':
    from multiprocessing import cpu_count

    cpu = cpu_count()

    parser = ArgumentParser(description='Bin MinION reads by (cumulative) intervals')
    parser.add_argument('-f', '--fastq', metavar='/basecalled/folder/',
                        required=True,
                        help='Input folder with fastq file(s), gzipped or not')
    parser.add_argument('-o', '--output', metavar='/qc/',
                        required=True,
                        help='Output folder')
    parser.add_argument('-t', '--threads', metavar='{}'.format(cpu),
                        required=False, default=cpu,
                        help='Number of CPU'
                             'Default: {}'.format(cpu))
    parser.add_argument('-i', '--interval', metavar='1h',
                        required=True,
                        default='1h',
                        help='Time interval to create the bins.'
                             'Use integer values (no decimal)'
                             'Can use \'h\', \'m\' or \'s\''
                             'Bins are cumulative. If interval is set to at 1h, second bin will'
                             'also contains the reads from the first bin.'
                             'Default: 1h')
    parser.add_argument('-p', '--prefix', metavar='my_output_name',
                        required=False,
                        default='interval',
                        help='Output file prefix'
                             'If using "my_sample" as prefix, files will be named "my_sample_1h.fastq.gz", etc.'
                             'Default: interval')

    # Get the arguments into an object
    arguments = parser.parse_args()

    NanoReadsBinner(arguments)
