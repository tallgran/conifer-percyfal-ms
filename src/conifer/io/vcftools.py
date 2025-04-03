"""File handling functions for vcftools-like files"""

from conifer.io.csv import BufferedCsvGroupReader


class BufferedVcftoolsReader(BufferedCsvGroupReader):  # pylint: disable=too-few-public-methods
    """Read vcftools output files.

    Read vcftools output files. Return chunks grouped by chromosome.
    Returned chunks can also be limited by chromosome position.


    """

    def __init__(
        self,
        file,
        chunksize=1000000,
        position_chunksize=None,
        **kwargs,
    ):
        super().__init__(file, chunksize, **kwargs)
        self._position_chunksize = position_chunksize

    def _update_buffer(self, group_name):
        if len(self._buffer) == 0:
            return
        self._buffer = self._buffer.drop(
            self._buffer[
                (self._buffer[self.groupby] == group_name[0])
                & (self._buffer["PCHUNK"] == group_name[1])
            ].index
        )

    def _can_return(self):
        self._buffer["PCHUNK"] = 0
        if len(self._buffer) > 0 and self._position_chunksize is not None:
            self._buffer["PCHUNK"] = (
                self._buffer["POS"] // self._position_chunksize
            )
        _, _, n_groups = self._process_group(
            self._buffer, [self.groupby, "PCHUNK"]
        )
        if n_groups > 1:
            return True
        return False

    def _return_value(self):
        group_name, group_result, _ = self._process_group(
            self._buffer, [self.groupby, "PCHUNK"]
        )
        self._update_buffer(group_name)
        return group_name[0], group_result.drop("PCHUNK", axis=1)
