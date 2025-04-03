"""File handling functions for BED-like files"""

# pylint: disable=too-many-nested-blocks,too-many-branches,too-many-statements
import re

import pandas as pd

from conifer.io.files import is_gzip, sniff_infile

CATEGORY_COLUMNS = [
    "CHROM",
    "chrom",
    "NAME",
    "name",
    "STRAND",
    "strand",
    "TYPE",
    "type",
    "SOURCE",
    "source",
]

COLUMN_MAP = {"CHROM": "chrom", "POS": "end", "PI": "score", "N_CHR": "score"}

BED3COLUMNS = ["chrom", "begin", "end"]
BED4COLUMNS = ["chrom", "begin", "end", "score"]
BED6COLUMS = ["chrom", "begin", "end", "name", "score", "strand"]


def read_fai_genomefile(file) -> pd.DataFrame:
    """Read genome file in fai format"""
    delimiter, _ = sniff_infile(file.name)
    compression = "gzip" if is_gzip(file.name) else None
    df = pd.read_csv(
        file.name, compression=compression, header=None, sep=delimiter
    )
    if len(df.columns) == 5 and re.search(r"(.fai|.fai.gz)", file.name):
        df.columns = ["chrom", "length", "offset", "linebases", "linewidth"]
        df["begin"] = 0
        df["end"] = df["length"]
        return df
    return None


class BufferedBedReader:  # pylint: disable=too-many-instance-attributes
    """Read BED-like files"""

    def __init__(  # pylint: disable=too-many-arguments
        self,
        file,
        chunksize=1000000,
        groupby=None,
        delimiter=None,
        columns=None,
        **kwargs,
    ):
        self.file = file
        self.chunksize = chunksize
        self.kwargs = kwargs
        if delimiter is None or columns is None:
            self.delimiter, self.columns = sniff_infile(file)
        else:
            self.delimiter = delimiter
            self.columns = columns
        self.names = [COLUMN_MAP.get(x, x) for x in self.columns]
        self.compression = "gzip" if is_gzip(file) else None
        self._buffer = []
        self._setup_groupby(groupby)
        self._setup_iterator()
        self._result = pd.DataFrame()

    def _setup_groupby(self, groupby):
        if groupby is not None:
            if groupby not in self.columns:
                raise ValueError(
                    f"Groupby column {groupby} not in columns {self.columns}"
                )
            self.groupby = groupby
        else:
            if COLUMN_MAP.get(self.columns[0]) in CATEGORY_COLUMNS:
                self.groupby = COLUMN_MAP.get(self.columns[0])
            else:
                self.groupby = 0

    @property
    def has_header(self):
        """Check whether file has header"""
        return any(not isinstance(x, int) for x in self.columns)

    def _setup_iterator(self):
        kwargs = self.kwargs.copy()
        kwargs["compression"] = self.compression
        kwargs["sep"] = self.delimiter
        kwargs["header"] = 0 if self.has_header else None
        kwargs["names"] = self.names
        kwargs["dtype"] = {
            COLUMN_MAP.get(x, x): "category"
            for x in self.columns
            if x in CATEGORY_COLUMNS
        }
        self.iterator = pd.read_csv(self.file, iterator=True, **kwargs)

    def __iter__(self):
        return self

    def __next__(self):
        try:
            result = self._buffer.pop(0)
        except IndexError:
            result = [None, pd.DataFrame()]
        while True:
            if len(self._buffer) > 0:
                break
            try:
                chunk = self.iterator.get_chunk(self.chunksize)
                chunk["begin"] = chunk["end"] - 1
            except StopIteration:
                break
            if chunk.empty:
                continue
            grouped = chunk.groupby(self.groupby, observed=True)
            for name, group in grouped:
                group[self.groupby] = group[
                    self.groupby
                ].cat.remove_unused_categories()
                if len(result[1]) == 0:
                    result[0] = name
                    result[1] = pd.concat([result[1], group])
                else:
                    if name == result[0]:
                        result[1] = pd.concat([result[1], group])
                    else:
                        if len(self._buffer) > 0:
                            if name == self._buffer[-1][0]:
                                self._buffer[-1][1] = pd.concat(
                                    [self._buffer[-1][1], group]
                                )
                            else:
                                self._buffer.append([name, group])
                        else:
                            self._buffer.append([name, group])

            if len(self._buffer) > 0:
                break
        if len(result[1]) == 0:
            raise StopIteration
        result[1] = result[1].reindex(BED4COLUMNS, axis=1)
        return result[0], result[1]
