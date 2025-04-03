"""Common file handling functions"""

import csv
import gzip
import logging
import os
import sys
from io import TextIOWrapper

logger = logging.getLogger(__name__)


common_delimeters = {repr(" "), repr("\t"), repr(",")}


def is_gzip(file):
    """Check whether file is gzipped"""
    ext = os.path.splitext(file)[1]
    return ext == ".gz"


def sniff_infile(file) -> str:
    """Sniff infile header and determine delimiter and column names."""
    has_header = False
    header = None
    fopen = gzip.open if is_gzip(file) else open
    with fopen(file, "rb") as f:
        dialect = csv.Sniffer().sniff(f.read(2048).decode())
        delimiter = dialect.delimiter
        if repr(delimiter) not in common_delimeters:
            logger.warning(
                "%s has uncommon delimiter '%s'; assuming tab", file, delimiter
            )
            delimiter = "\t"
        f.seek(0)
        if csv.Sniffer().has_header(f.read(2048).decode()):
            logger.debug("%s has a header", file)
            has_header = True
    if has_header:
        with fopen(file, "r") as fh:
            if is_gzip(file):
                reader = csv.reader(TextIOWrapper(fh), dialect=dialect)
            else:
                reader = csv.reader(fh, dialect=dialect)
            header = next(reader)
    return delimiter, header


def _sizeof_fmt(num, suffix="B"):
    for unit in ["", "K", "M", "G", "T", "P", "E", "Z"]:
        if abs(num) < 1024.0:
            return f"{num:.2f} {unit}{suffix}"
        num /= 1024.0
    return f"{num:.2f} Y{suffix}"


def _size(obj):
    return sys.getsizeof(obj)


def size(obj, humanize=True):
    """Return the size of an object"""
    if humanize:
        return _sizeof_fmt(_size(obj))
    return _size(obj)
