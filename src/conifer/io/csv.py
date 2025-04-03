"""File handling functions for CSV-like files"""

import pandas as pd

from conifer.io.files import is_gzip, sniff_infile


class BufferedCsvReader:
    """Read CSV-like files."""

    def __init__(self, file, chunksize=1000000, **kwargs):
        delimiter = kwargs.get("sep")
        columns = kwargs.get("names")
        names = columns
        if delimiter is None or columns is None:
            delimiter, names = sniff_infile(file)
            columns = names if columns is None else columns
        self.kwargs = {
            "compression": "gzip" if is_gzip(file) else None,
            "sep": delimiter,
            **kwargs,
        }
        if names is not None and "header" not in self.kwargs:
            self.kwargs["header"] = 0
        self.chunksize = chunksize
        self.file = file
        self.columns = columns
        self._buffer = pd.DataFrame()
        self._setup_iterator()

    def _setup_iterator(self):
        self.iterator = pd.read_csv(self.file, iterator=True, **self.kwargs)

    def __iter__(self):
        return self

    def __next__(self):
        return self.iterator.get_chunk(self.chunksize)

    def close(self):
        """Flush and close this stream."""
        self.iterator.close()


def _process_group(df, groupby, sortby=None):
    if len(df) == 0:
        return None, None, 0
    if not isinstance(groupby, list):
        groupby = [groupby]
    sortindex = 0 if sortby is None else groupby.index(sortby)
    sorting_key = df[groupby[sortindex]].unique().tolist()

    def _sort_group(x):
        if isinstance(x, tuple):
            return sorting_key.index(x[sortindex])
        return sorting_key.index(x)

    grouped = df.groupby(groupby, observed=True, sort=False)
    group_name = list(sorted(grouped.groups.keys(), key=_sort_group))[0]
    return group_name, grouped.get_group(group_name), len(grouped)


class BufferedCsvGroupReader(BufferedCsvReader):
    """Read CSV-like files. Return chunks grouped by column of choice."""

    def __init__(
        self,
        file,
        chunksize=1000000,
        groupby=None,
        **kwargs,
    ):
        super().__init__(file, chunksize, **kwargs)
        self._setup_groupby(groupby)

    def _setup_groupby(self, groupby):
        if groupby is not None:
            if groupby not in self.columns:
                raise ValueError(
                    f"Groupby column {groupby} not in columns {self.columns}"
                )
            self.groupby = groupby
        else:
            self.groupby = 0 if self.columns is None else self.columns[0]

    def _process_group(self, df, groupby, sortby=None):
        return _process_group(df, groupby, sortby)

    def _append_to_buffer(self):
        success = False
        try:
            data = super().__next__()
            if data is not None:
                if isinstance(data, tuple):
                    _, data = data
                self._buffer = pd.concat([self._buffer, data])
                success = True
        except StopIteration:
            success = False
        return success

    def _can_return(self):
        _, _, n_groups = self._process_group(self._buffer, self.groupby)
        if n_groups > 1:
            return True
        if self.chunksize is not None:
            if len(self._buffer) >= self.chunksize:
                return True
        return False

    def _return_value(self):
        group_name, group_result, _ = self._process_group(
            self._buffer, self.groupby
        )
        if self.chunksize is not None:
            imax = min(self.chunksize, len(group_result))
        else:
            imax = len(group_result)
        result = self._buffer.iloc[:imax]
        self._buffer = self._buffer.iloc[imax:]
        return group_name, result

    def __next__(self):
        while not self._can_return():
            success = self._append_to_buffer()
            if not success:
                if len(self._buffer) > 0:
                    return self._return_value()
                raise StopIteration
            if len(self._buffer) == 0:
                raise StopIteration
        return self._return_value()
