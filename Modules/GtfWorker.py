#!/usr/bin/env python3
from collections import defaultdict
import gzip
import pandas as pd
import re
import pdb

class Gtf_file():
    def __init__(self, gtf_file):
        self.gtf_file = gtf_file
        self.GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score','strand', 'frame']
        self.R_SEMICOLON = re.compile(r'\s*;\s*')
        self.R_COMMA     = re.compile(r'\s*,\s*')
        self.R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')

    def dataframe(self):
        filename = self.gtf_file
        """Open an optionally gzipped GTF file and return a pandas.DataFrame.
        """
        # Each column is a list stored as a value in this dict.
        result = defaultdict(list)

        for i, line in enumerate(lines(filename)):
            for key in line.keys():
                # This key has not been seen yet, so set it to None for all
                # previous lines.
                if key not in result:
                    result[key] = [None] * i

            # Ensure this row has some value for each column.
            for key in result.keys():
                result[key].append(line.get(key, None))

        return pd.DataFrame(result)


    def lines(self):
        filename = self.gtf_file
        """Open an optionally gzipped GTF file and generate a dict for each line.
        """
        fn_open = gzip.open if filename.endswith('.gz') else open
        with fn_open(filename,'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                else:
                    yield self.parse(line)


    def parse(self,line):
        """Parse a single GTF line and return a dict.
        """
        result = {}

        fields = line.rstrip().split('\t')

        for i, col in enumerate(self.GTF_HEADER):
            result[col] = self._get_value(fields[i])

        # INFO field consists of "key1=value;key2=value;...".
        infos = [x for x in re.split(self.R_SEMICOLON, fields[8]) if x.strip()]

        for i, info in enumerate(infos, 1):
            # It should be key="value".
            try:
                key, _, value = re.split(self.R_KEYVALUE, info, 1)
            # But sometimes it is just "value".
            except ValueError:
                key = 'INFO{}'.format(i)
                value = info
            # Ignore the field if there is no value.
            if value:
                result[key] = self._get_value(value)

        return result


    def _get_value(self,value):
        if not value:
            return None

        # Strip double and single quotes.
        value = value.strip('"\'')

        # Return a list if the value has a comma.
        if ',' in value:
            value = re.split(self.R_COMMA, value)
        # These values are equivalent to None.
        elif value in ['', '.', 'NA']:
            return None

        return value
