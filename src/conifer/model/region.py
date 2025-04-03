"""Model of genomic regions"""

import logging
import math
from dataclasses import InitVar, asdict, dataclass
from typing import List, Type, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

HARMONIC_NUMBER_MAX = 2500

# Harmonic number. For large n, the naive function is slow and an
# approximation should be used (cf
# https://stackoverflow.com/questions/404346/python-program-to-calculate-harmonic-series).
# The largest n we see is <2100 so we can calculate a lookup
# dictionary to speed up calculations.
HARMONIC_NUMBER = np.zeros(HARMONIC_NUMBER_MAX + 1)
HARMONIC_NUMBER_SQUARED = np.zeros(HARMONIC_NUMBER_MAX + 1)
HARMONIC_NUMBER[0] = np.nan
HARMONIC_NUMBER_SQUARED[0] = np.nan
for ii in np.arange(1, HARMONIC_NUMBER_MAX):
    HARMONIC_NUMBER[ii] = np.sum(1 / np.arange(1, ii))
    HARMONIC_NUMBER_SQUARED[ii] = np.sum(1 / (np.arange(1, ii) ** 2))


@dataclass
class Region:
    """Region class. Represents a genomic region. Intervals are
    0-based, half-open."""

    chrom: str
    begin: int
    end: int

    @property
    def width(self):
        """Region width"""
        return self.end - self.begin

    @property
    def pos(self):
        """Region positions in 0-based coordinates."""
        return np.arange(self.begin, self.end)

    def overlaps(self, other):
        """Check whether region overlaps with other region"""
        if self.chrom != other.chrom:
            return False
        return (self.begin <= other.begin < self.end) or (
            other.begin <= self.begin < other.end
        )

    def is_disjoint(self, other):
        """Check whether region is disjoin with other region"""
        return not self.overlaps(other)

    def contains(self, other):
        """Check whether region contains other region"""
        return (self.begin <= other.begin) and (self.end >= other.end)

    def trim(self, other):
        """Trim region to other region"""
        return Region(
            self.chrom, max(self.begin, other.begin), min(self.end, other.end)
        )

    def as_region(self):
        """Return self as Region"""
        return Region(self.chrom, self.begin, self.end)

    def __str__(self):
        return f"{self.chrom}:{self.begin}-{self.end}"

    def __len__(self):
        return self.width

    def __lt__(self, other):
        return self.begin < other.begin

    def __gt__(self, other):
        return self.begin > other.begin

    def __le__(self, other):
        return self.begin <= other.begin

    def __ge__(self, other):
        return self.begin >= other.begin

    def __add__(self, other):
        if not self.overlaps(other):
            raise ValueError(f"Regions {self} and {other} do not overlap")
        return Region(
            self.chrom, min(self.begin, other.begin), max(self.end, other.end)
        )

    def __sub__(self, other):
        if not self.overlaps(other):
            raise ValueError(f"Regions {self} and {other} do not overlap")
        if self == other:
            return None
        if self.begin < other.begin:
            return Region(self.chrom, self.begin, other.begin)
        return Region(self.chrom, other.end, self.end)

    def __or__(self, other):
        if not self.overlaps(other):
            raise ValueError(f"Regions {self} and {other} do not overlap")
        return Region(
            self.chrom, min(self.begin, other.begin), max(self.end, other.end)
        )

    def __and__(self, other):
        if not self.overlaps(other):
            raise ValueError(f"Regions {self} and {other} do not overlap")
        return Region(
            self.chrom, max(self.begin, other.begin), min(self.end, other.end)
        )

    def __getitem__(self, item):
        if not isinstance(item, slice):
            item = slice(item, item + 1)
        if not self.contains(
            Region(self.chrom, self.begin + item.start, self.begin + item.stop)
        ):
            raise ValueError(
                f"Region {self} ({type(self)}) does not contain {item}"
            )
        kw = asdict(self)
        kw["begin"] = self.begin + item.start
        kw["end"] = self.begin + item.stop
        return type(self)(**kw)

    def __eq__(self, other):
        return (
            (self.chrom == other.chrom)
            and (self.begin == other.begin)
            and (self.end == other.end)
        )


@dataclass
class Bed4Region(Region):
    """BED4 region class. Represents a region with a name"""

    name: str

    def __str__(self):
        return f"{self.chrom}:{self.begin}-{self.end} ({self.name})"

    def trim(self, other):
        """Trim region to other region"""
        r = super().trim(other)
        return type(self)(r.chrom, r.begin, r.end, self.name)


@dataclass
class Bed6Region(Bed4Region):
    """BED6 region class."""

    score: float
    strand: str

    def __post_init__(self):
        if self.strand is None:
            self.strand = "."
        if self.strand not in ["+", "-", "."]:
            raise ValueError(
                f"Strand must be '+', '-' or '.', not {self.strand}"
            )

    def __str__(self):
        return (
            f"{self.chrom}\t{self.begin}\t{self.end}\t"
            f"{self.name}\t{self.score}\t{self.strand}"
        )


@dataclass
class Score(Bed6Region):
    """Score class. Extended bed format with an extra column that
    holds the number of observed chromosomes for a site.

    """

    n_chr: int = 0

    def __str__(self):
        return (
            f"{self.chrom}\t{self.begin}\t{self.end}\t"
            f"{self.name}\t{self.score}\t{self.strand}\t{self.n_chr}"
        )


@dataclass
class Mask(Bed4Region):
    """Mask class"""

    mask: InitVar[np.ndarray] = None
    invert: InitVar[bool] = False

    def __post_init__(self, mask, invert):
        if mask is None:
            fun = np.zeros if invert else np.ones
            self.mask = fun(self.width, dtype=np.uint8)
        else:
            if isinstance(mask, (np.uint8, np.int8)):
                mask = np.array([mask], dtype=np.uint8)
            assert len(mask) == self.width, (
                f"passed mask length ({len(mask)}) must be equal in length to"
                f"region ({self.width})"
            )
            self.mask = mask

    @property
    def invert_mask(self):
        """Invert mask"""
        return np.logical_not(self.mask).astype(np.uint8)

    @property
    def n_unmasked(self):
        """Number of unmasked (0) sites."""
        return np.sum(self.invert_mask)

    @property
    def n_masked(self):
        """Number of masked (1) sites"""
        return np.sum(self.mask)

    def unmask(self, region: Region, invert: bool = False):
        """Invert the mask based on a list of regions"""
        assert issubclass(
            type(region), Region
        ), f"region must subclass Region; {type(region)}"
        maskval = 1 if invert else 0
        begin = max(region.begin - self.begin, 0)
        end = min(region.end - self.begin, self.end - self.begin)
        self.mask[begin:end] = maskval
        return self

    def unmask_r(self, regions: List[Region], invert: bool = False):
        """Unmask a set of regions"""
        logger.debug("unmasking mask %s from %i regions", self, len(regions))
        for r in regions:
            self.unmask(r, invert)
        return self

    def trim(self, other):
        """Trim region to other region"""
        r = Region.trim(self, other)
        return type(self)(
            r.chrom,
            r.begin,
            r.end,
            self.name,
            self.mask[(r.begin - self.begin) : (r.end - self.begin)],
        )

    def __repr__(self):
        return repr(self.mask)

    def __or__(self, other):
        if not Region.__eq__(self, other):
            raise ValueError(f"Region of mask {self} must be equal to {other}")
        m = super().__or__(other)
        name = f"{self.name}|{other.name}"
        return Mask(
            m.chrom,
            m.begin,
            m.end,
            name,
            np.logical_or(self.mask, other.mask).astype(np.uint8),
        )

    def __and__(self, other):
        if not Region.__eq__(self, other):
            raise ValueError(f"Region of mask {self} must be equal to {other}")
        m = super().__and__(other)
        name = f"{self.name}&{other.name}"
        return Mask(
            m.chrom,
            m.begin,
            m.end,
            name,
            np.logical_and(self.mask, other.mask).astype(np.uint8),
        )

    def __eq__(self, other):
        return Region.__eq__(self, other) & np.array_equal(
            self.mask, other.mask
        )

    def __getitem__(self, item):
        m = super().__getitem__(item)
        return type(self)(m.chrom, m.begin, m.end, self.name, self.mask[item])


@dataclass
class Feature(Mask):
    """Feature mask class"""

    def __repr__(self):
        return repr(self.mask)


RegionLike = Type[
    Union[
        Region,
        Bed4Region,
        Bed6Region,
        Score,
        Mask,
        Feature,
    ]
]


@dataclass
class ScoreVector(Bed4Region):
    """Score vector class that flattens a list of Score objects into
    numpy arrays.

    Holds a score vector and corresponding positions for segregating
    sites. In addition, may hold information about number of
    chromosomes, feature and accessibility mask in region.

    """

    __slots__ = [
        "_score",
        "_score_pos",
        "_n_chr",
        "_unmaskary",
    ]

    score_vector: InitVar[list[Score]] = None

    def __post_init__(
        self,
        score_vector: list[Score],
    ):
        if score_vector is None:
            score_vector = []
        num_overlap = int(np.sum([self.overlaps(s) for s in score_vector]))
        self._score = np.zeros(num_overlap, dtype=np.float32)
        self._score_pos = np.zeros(num_overlap, dtype=np.uint32)
        self._n_chr = np.zeros(num_overlap, dtype=np.uint16)
        i = 0
        for s in score_vector:
            if not self.contains(s):
                continue
            self.score[i] = s.score
            self.score_pos[i] = s.begin
            self.n_chr[i] = s.n_chr
            i += 1
        self._unmaskary = np.ones(len(self.score), dtype=np.uint8).astype(bool)

    @property
    def score(self):
        """Score"""
        return self._score

    @score.setter
    def score(self, value):
        """Set score"""
        self._score = value

    @property
    def score_pos(self):
        """Score position in 0-based coordinates."""
        return self._score_pos

    @score_pos.setter
    def score_pos(self, value):
        """Set score position"""
        self._score_pos = value

    @property
    def n_chr(self):
        """Number of chromosomes"""
        return self._n_chr

    @n_chr.setter
    def n_chr(self, value):
        """Set number of chromosomes"""
        self._n_chr = value

    @property
    def score_pos_mask(self):
        """Get score positions as a mask, where 0 indicates the score
        position for the entire region"""
        mask = Mask(
            chrom=self.chrom, begin=self.begin, end=self.end, name=self.name
        )
        mask.mask[self.score_pos - self.begin] = 0
        return mask

    @property
    def as_scores(self) -> List[Score]:
        """Return scores as list of Scores without specific name and
        strand annotations."""
        result = []
        for pos, score, n_chr in zip(self.score_pos, self.score, self.n_chr):
            result.append(
                Score(self.chrom, pos, pos + 1, self.name, score, ".", n_chr)
            )
        return result

    def __getitem__(self, item):
        sv = super().__getitem__(item)
        return ScoreVector(
            sv.chrom, sv.begin, sv.end, self.name, self.as_scores
        )


@dataclass
class RegionFeatureScore(Bed4Region):
    """Region score class

    Holds scores in a region for a given feature. The scores are
    stored in a numpy array along with a vector that holds the
    positions in genomic coordinates.

    Mandatory feature masks and optional accessibilty masks are taken
    into account to calculate


    In addition to the scores, the class also holds information on
    regional masks."""

    accessibility_mask: Mask
    feature_mask: Mask
    region_mask: Mask
    score_vector: ScoreVector


@dataclass
class AnalysisChunk(Bed4Region):
    """Analysis chunk class.

    Holds necessary information to generate FeatureScore that combine
    scores, regions, accessibility masks and feature masks.

    """

    __slots__ = ["_score", "_mask", "_feature", "_region"]

    def __post_init__(self):
        self._score = pd.DataFrame({})
        self._mask = {}
        self._feature = {}
        self._region = {}

    @property
    def feature(self):
        """Feature dictionary consisting of Mask objects."""
        return self._feature

    @property
    def mask(self):
        """Mask dictionary consisting of Mask objects."""
        return self._mask

    def _cmask(self, logical_and=True):
        """Return combination of masks"""
        m = None
        for v in self.mask.values():
            if m is None:
                m = v
            else:
                if logical_and:
                    m = v & m
                else:
                    m = v | m
        assert (
            len(m) == self.width
        ), f"Mask length ({len(m)} not equal to region width ({self.width})"
        return m

    @property
    def accessibility_mask(self):
        """Return combination of all masks. Kept for backwards
        compatibility."""
        return self._cmask(logical_and=False)

    @property
    def amask(self):
        """Return AND combination of all masks"""
        return self._cmask()

    @property
    def omask(self):
        """Return OR combination of all masks"""
        return self._cmask(False)

    @property
    def region(self) -> List[RegionLike]:
        """Region dictionary consisting of List[RegionLike] values."""
        return self._region

    def region_mask(self, key) -> Mask:
        """Return mask for region"""
        r = self._region[key]
        return Mask(self.chrom, self.begin, self.end, key).unmask_r(r)

    @property
    def score(self) -> pd.DataFrame:
        """Score data frame"""
        return self._score

    @score.setter
    def score(self, value):
        if not isinstance(value, pd.DataFrame):
            raise TypeError(f"Score must be a DataFrame, not {type(value)}")
        self._score = value

    def add_mask(self, label: str, df: pd.DataFrame, invert: bool = False):
        """Add mask"""
        mask = dataframe_to_region_list(df, cls=Bed4Region, name=label)
        logger.info("adding mask %s from %i regions", label, len(mask))
        self._mask[label] = Mask(
            self.chrom, self.begin, self.end, self.name, invert=invert
        ).unmask_r(mask, invert=invert)

    def add_feature(self, label: str, df: pd.DataFrame, invert: bool = False):
        """Add feature"""
        feature = dataframe_to_region_list(df, cls=Bed4Region, name=label)
        logger.info("adding feature %s from %i regions", label, len(feature))
        ft = Feature(
            self.chrom, self.begin, self.end, self.name, invert=invert
        )
        if len(feature) > 0:
            self._feature[label] = ft.unmask_r(feature, invert=invert)
        else:
            self._feature[label] = ft

    def add_region(self, label: str, df: pd.DataFrame):
        """Add region"""
        region = dataframe_to_region_list(df, cls=Bed4Region)
        logger.info("adding region %s from %i regions", label, len(region))
        self._region[label] = region

    def describe(self):
        """Describe the chunk"""
        s = "\n\n" + str(self) + "\n\n"
        data = []
        for k, v in self._mask.items():
            data.append(["Mask", k, len(v)])
        for k, v in self._feature.items():
            data.append(["Feature", k, len(v)])
        for k, v in self._region.items():
            data.append(["Region", k, len(v)])
        s += repr(pd.DataFrame(data, columns=["Type", "Name", "Count"]))
        return s

    def __getitem__(self, item):
        region = super().__getitem__(item)
        result = type(self)(region.chrom, region.begin, region.end, self.name)
        result.score = subset_dataframe_by_region(self.score, result)

        for label, mask in self._mask.items():
            result._mask[label] = mask[item.start : item.stop]
        for label, ft in self._feature.items():
            result._feature[label] = ft[item.start : item.stop]
        return result


def calculate_summary_stats(  # pylint: disable=too-many-locals,too-many-arguments
    region: Region,
    mask: Mask,
    feature: Feature,
    score_vector: np.ndarray,
    n_chr: np.ndarray,
    score_pos: np.ndarray,
) -> np.ndarray:
    """Calculate summary statistics for region

    Calculate and report the following statistics:

    score - sum of scores for accessible sites within region and feature
    total_score - sum of scores for all sites within region and feature
    site_score - score divided by accessible sites within region and feature
    n_sites - total number of sites in region and feature, segregating
        and monomorphic
    n_accessible - total number of accessible sites in region and feature,
        segregating and monomorphic
    n_segregating_sites - total number of segregating sites in region
        and feature
    n_segregating_sites_accessible - total number of accessible
        segregating sites in region and feature
    theta_W - Watterson's theta
    tajimas_D - tajimas D

    Importantly, the masks must be inverted, i.e., 0 denotes sites to
    be removed, 1 denotes sites to keep!

    """
    allsites = feature
    i_score_all = allsites.mask[score_pos - region.begin].astype(bool)

    unmasked = mask & allsites
    i_score_accessible = unmasked.mask[score_pos - region.begin].astype(bool)

    n_sites = np.sum(allsites.mask)
    n_accessible = np.sum(unmasked.mask)
    n_segregating_sites = len(score_vector[i_score_all])
    n_segregating_sites_accessible = len(score_vector[i_score_accessible])

    score = np.nansum(score_vector[i_score_accessible])
    total_score = np.nansum(score_vector[i_score_all])
    site_score = score / n_accessible
    try:
        theta = wattersons_theta(
            n_chr[i_score_accessible],
            n_accessible,
            n_segregating_sites_accessible,
        )
    except TypeError:
        theta = np.nan
    try:
        theta_abs = wattersons_theta(
            n_chr[i_score_accessible],
            n_accessible,
            n_segregating_sites_accessible,
            True,
        )
    except TypeError:
        theta_abs = np.nan
    except ValueError:
        taj_d = np.nan
    try:
        taj_d = tajimas_d(
            n_chr[i_score_accessible],
            total_score,
            theta_abs,
            n_segregating_sites_accessible,
        )
    except TypeError:
        taj_d = np.nan
    except ValueError:
        taj_d = np.nan

    return (  # pylint: disable=duplicate-code
        score,
        total_score,
        site_score,
        n_sites,
        n_accessible,
        n_segregating_sites,
        n_segregating_sites_accessible,
        theta,
        taj_d,
    )


def wattersons_theta(
    n_chr: np.array, n: int, s: int, absolute: bool = False
) -> float:
    """Calculate Watterson's theta"""
    if len(n_chr) == 0:
        return np.nan
    if np.all(n_chr == 0):
        return np.nan
    denominator = np.nanmean(np.array([HARMONIC_NUMBER[x] for x in n_chr]))
    if not absolute:
        denominator = denominator * n
    return np.divide(s, denominator)


def tajimas_d(
    n_chr: np.array, score: np.array, theta: np.array, s: int
) -> float:
    """Calculate Tajima's D over a region"""
    if np.any(n_chr == 0) or len(n_chr) == 0:
        return np.nan
    a1 = np.nanmean(np.array([HARMONIC_NUMBER[x] for x in n_chr]))
    a2 = np.nanmean(np.array([HARMONIC_NUMBER_SQUARED[x] for x in n_chr]))
    n = np.nanmean(n_chr)
    b1 = np.divide(n + 1, 3 * (n - 1))
    b2 = np.divide(2 * (n**2 + n + 3), 9 * n * (n - 1))
    c1 = b1 - np.divide(1, a1)
    c2 = b2 - np.divide(n + 2, a1 * n) + np.divide(a2, a1**2)
    e1 = np.divide(c1, a1)
    e2 = np.divide(c2, a1**2 + a2)
    d = np.divide(
        score - theta,
        math.sqrt(e1 * s + e2 * s * (s - 1)),
    )
    return d


def dataframe_to_region_list(
    df: pd.DataFrame, *, cls: RegionLike = Region, **kw
) -> List[RegionLike]:
    """Convert dataframe to list of Region-like objects"""
    result = []
    try:
        result = [cls(**{**x, **kw}) for _, x in df.iterrows()]
    except TypeError as e:
        logger.warning(e)
        logger.info("Reverting to Region class")
        result = [Region(**{**x}) for _, x in df.iterrows()]
    return result


def subset_dataframe_by_region(
    df: pd.DataFrame,
    region: RegionLike,
    overlap: bool = True,
    trim: bool = True,
) -> pd.DataFrame:
    """Subset dataframe by region."""
    if df is None or len(df) == 0:
        return df
    if overlap:
        i = (
            ((df.begin <= region.begin) & (region.begin < df.end))
            | ((region.begin <= df.begin) & (df.begin < region.end))
        ) & (df.chrom == region.chrom)
    else:
        i = ((df.begin >= region.begin) & (df.end < region.end)) & (
            df.chrom == region.chrom
        )
    df = df[i].copy()
    if trim:
        x = df.begin.apply(lambda x: max(x, region.begin))
        df.begin = x
        x = df.end.apply(lambda x: min(x, region.end))
        df.end = x

    return df


__all__ = [
    "Region",
    "Bed4Region",
    "Bed6Region",
    "Score",
    "Mask",
    "Feature",
    "RegionFeatureScore",
    "AnalysisChunk",
    "dataframe_to_region_list",
    "subset_dataframe_by_region",
]
