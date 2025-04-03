"""Presence-absence coverage analysis helper classes and functions."""

import itertools
import logging
from dataclasses import dataclass, field
from pathlib import Path

import holoviews as hv
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn import model_selection
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA as PrinComp  # noqa: N811
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import RepeatedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def load_coverage_data(
    samples,
    feature="gene",
    template="data/presabs/{feature}/{sample}.sum.tsv",
    home=Path("."),
    chrom_regex=r"PA_chr\d+_P?G",
    gene_filter=None,
):
    """Load coverage data from multiple samples and return a data frame."""

    def load_sample(sample, chrom_index=None):
        infile = str.format(
            str(home / template), feature=feature, sample=sample
        )
        columns = ["feature", "coverage", "width"]
        df = pd.read_table(infile, header=None, names=columns)
        df["gene"] = df["feature"].map(lambda x: x.split("|")[0])
        if gene_filter is not None:
            df = df[~df["gene"].isin(gene_filter)]
        df["feature"] = df["feature"].map(lambda x: x.split("|")[1])
        df.set_index(["gene", "feature", "width"], inplace=True)
        if chrom_index is None:
            return df
        return df[chrom_index]

    chrom_index = (
        load_sample(samples[0])
        .index.get_level_values("feature")
        .str.match(chrom_regex)
    )
    data = pd.concat(
        [load_sample(s, chrom_index=chrom_index) for s in samples], axis=1
    )
    data.columns = samples
    return data


@dataclass
class Annotation:
    data: pd.DataFrame = field(default_factory=pd.DataFrame)
    columns: list[str] = field(default_factory=list)

    def __post_init__(self):
        if len(self.columns) == 0:
            self.columns = self.data.columns.values

    @property
    def annotation(self):
        return self.data[self.columns]

    @property
    def empty(self):
        return self.data.empty


@dataclass
class Coverage:
    """Container class for coverage data."""

    data: pd.DataFrame
    feature: str = "gene"
    type: str = "raw"
    randomized: bool = False
    samples: list[str] = field(default_factory=list)
    samplesets: dict = None
    annotation: Annotation = field(default_factory=Annotation)
    columns: list[str] = field(default_factory=list)
    n_features: int = None

    def __post_init__(self):
        assert isinstance(self.data, pd.DataFrame), (
            "data argument must be of type %s; saw %s"
            % (pd.DataFrame, type(self.data))
        )
        if self.type == "raw":
            self.n_features = self.data.shape[0]
            by = "gene"
            width = (
                self.data.index.to_frame(index=False)
                .groupby(by)
                .agg("sum")["width"]
            )
            _data = self.data.groupby(by).agg("sum")
            _data["width"] = width
            self.data = _data.reset_index().set_index(["gene", "width"])

        if (
            len(self.data.columns.names) == 1
            and self.data.columns.names[0] is None
        ):
            self.data.columns.names = ["sample"]
        self.samples = self.data.columns.get_level_values("sample")
        self._genes = self.data.index.get_level_values("gene")
        if self.samplesets is None:
            self.samplesets = {"all": self.samples}
        revmap = {
            item: k
            for k, v in self.samplesets.items()
            for item in self.samplesets[k]
        }
        ss = [revmap[s] for s in self.samples]
        if self.data.columns.names != ["sample", "sampleset"]:
            self.data.columns = pd.MultiIndex.from_tuples(
                zip(self.samples, ss), names=["sample", "sampleset"]
            )
        self.columns = (
            None if self.annotation.empty else self.annotation.columns
        )

    def normalize(self):
        """Divide values by length of feature."""
        if self.type == "normalized":
            return self
        if self.type != "raw":
            raise Exception("Can only normalize raw coverage data")
        _data = self.data.div(
            self.data.index.get_level_values("width"), axis=0
        )
        return Coverage(
            _data,
            feature=self.feature,
            type="normalized",
            n_features=self.n_features,
            annotation=self.annotation,
        )

    def discretize(self, threshold=0.2):
        """Discretize by gene"""
        return Coverage(
            (self.normalize().data > threshold).astype(np.int8),
            feature=self.feature,
            type="discretized",
            n_features=self.n_features,
            annotation=self.annotation,
        )

    def randomize(self, seed=42):
        """Randomize values of every row independently."""
        np.random.seed(seed)
        data = pd.DataFrame(
            np.apply_along_axis(np.random.permutation, axis=1, arr=self.data)
        )
        data.columns = self.data.columns
        data.index = self.data.index
        return Coverage(
            data,
            type=self.type,
            feature=self.feature,
            randomized=True,
            n_features=self.n_features,
            annotation=self.annotation,
        )

    @property
    def genes(self):
        return self._genes

    @property
    def n_genes(self):
        return len(self._genes)

    def subset(self, index):
        """Subset data and return new Coverage object. Indexes must align."""
        assert (
            index.names == self.data.index.names
        ), "subsetting with malaligned index"
        data = self.data.loc[index]
        return Coverage(
            data,
            type=self.type,
            feature=self.feature,
            randomized=self.randomized,
            n_features=self.n_features,
            annotation=self.annotation,
        )

    @property
    def adata(self):
        """Show data with annotations in index"""
        self.annotation.columns = self.columns
        annotation = (
            pd.MultiIndex.from_frame(
                self.data.index.to_frame().join(self.annotation.annotation)
            )
            .to_frame()
            .set_index(["gene"])
        )
        return annotation

    def mean(self, axis=1, by="sampleset"):
        data = self.data.groupby(by, axis=axis).agg(np.mean)
        return data

    def sum(self, axis=1, by="sampleset"):
        data = self.data.groupby(by, axis=axis).agg(np.sum)
        return data

    def summary(self, axis=1, by="sampleset"):
        if self.type == "discretized":
            return self.sum(axis=axis, by=by)
        return self.mean(axis=axis, by=by)


@dataclass
class PCBase:
    """Base class for PCA and PLS"""

    coverage: Coverage = field(default_factory=Coverage)
    n_components: int = 10

    _label = "PC"

    def __post_init__(self):
        assert isinstance(self.coverage, Coverage), (
            "coverage argument must be of type %s; saw %s"
            % (Coverage, type(self.coverage))
        )
        self._model = None
        self._model_m = pd.DataFrame()
        self._model_loadings = pd.DataFrame()
        self.model()

    @property
    def data(self):
        return self.coverage

    @property
    def pcs(self):
        return self._model_m

    @property
    def loadings(self):
        """Retrieve the loadings only"""
        columns = [
            x
            for x in self._model_loadings.columns
            if x.startswith(self._label)
        ]
        return self._model_loadings[columns]

    @property
    def samples(self):
        return self.coverage.samples

    @property
    def genes(self):
        return self.coverage.genes

    @property
    def feature(self):
        return self.coverage.feature

    @property
    def type(self):
        return self.coverage.type

    @property
    def randomized(self):
        return self.coverage.randomized

    @property
    def n_features(self):
        return self.coverage.n_features

    @property
    def n_genes(self):
        return self.coverage.n_genes

    def model(self):
        raise NotImplementedError


@dataclass
class PCA(PCBase):
    """PCA class. Perform PCA on coverage object"""

    def model(self):
        columns = [f"{self._label}{i+1}" for i in np.arange(self.n_components)]
        self._model = PrinComp(n_components=self.n_components)
        m = self.coverage.data
        x = m.T
        x = StandardScaler().fit_transform(x)
        pcs = self._model.fit_transform(x)
        self._model_m = pd.DataFrame(data=pcs, columns=columns)
        self._model_m["sample"] = self.samples
        self._model_m.set_index(["sample"], inplace=True)
        self._set_loadings()

    def _set_loadings(self):
        data = self._model_m
        columns = [x for x in data if x.startswith(self._label)]
        self._model_loadings = pd.DataFrame(
            self._model.components_.T, columns=columns, index=self.genes
        )

    def explained_variance_ratio(self, comp):
        return self._model.explained_variance_ratio_[comp]


@dataclass
class PLSDA(PCBase):
    """PLS-DA class. Perform PLS-DA on coverage object"""

    _label = "PLS"

    def setup_data(self):
        m = self.coverage.data
        x = m.T
        y = np.unique(
            self.coverage.data.columns.get_level_values("sampleset"),
            return_inverse=True,
        )[1]
        return x, y

    def model(self):
        columns = [f"{self._label}{i+1}" for i in np.arange(self.n_components)]
        self._model = PLSRegression(
            n_components=self.n_components, scale=False
        )
        x, y = self.setup_data()
        pcs = self._model.fit_transform(
            X=StandardScaler().fit_transform(x), y=y
        )
        self._model_m = pd.DataFrame(data=pcs[0], columns=columns)
        self._model_m["sample"] = self.samples
        self._model_m.set_index(["sample"], inplace=True)
        self._set_loadings()

    def _set_loadings(self):
        data = self._model_m
        columns = [x for x in data if x.startswith(self._label)]
        self._model_loadings = pd.DataFrame(
            self._model.x_loadings_, columns=columns, index=self.genes
        )

    def cross_validation(
        self, n_repeats=10, max_comp=20, scoring="neg_mean_squared_error", **kw
    ):
        x, y = self.setup_data()
        n_splits = kw.pop("n_splits", len(x))
        random_state = kw.pop("random_state", 42)
        cv = RepeatedKFold(
            n_splits=n_splits, n_repeats=n_repeats, random_state=random_state
        )
        loss = []
        predict = []
        logger.info(
            (
                "Running cross-validation with parameters n_splits=%d, "
                "n_repeats=%d, max_comp=%d"
            ),
            n_splits,
            n_repeats,
            max_comp,
        )
        for i in tqdm(np.arange(1, max_comp)):
            pls = PLSRegression(n_components=i, scale=False)
            score = (
                -1
                * model_selection.cross_val_score(
                    pls,
                    StandardScaler().fit_transform(x),
                    y,
                    cv=cv,
                    scoring=scoring,
                ).mean()
            )
            pred = model_selection.cross_val_predict(
                pls, StandardScaler().fit_transform(x), y
            )
            predict.append(pred)
            loss.append(score)

        return loss, predict

    def rmse(self, n_comp, test_size=0.2):
        x, y = self.setup_data()
        x_train, x_test, y_train, y_test = train_test_split(
            x, y, test_size=test_size
        )
        pls = PLSRegression(n_components=n_comp, scale=False)
        pls.fit(StandardScaler().fit_transform(x_train), y_train)
        y_pred = pls.predict(StandardScaler().fit_transform(x_test))
        return mean_squared_error(y_test, y_pred, squared=False)


@dataclass
class Loadings:
    model: PCBase = field(default_factory=PCBase)

    def __post_init__(self):
        self.data = self.model.loadings
        self._l_columns = [
            x for x in self.data.columns if x.startswith(self.model._label)
        ]
        self.reset()

    def plot_settings(self, plot, element):
        p = plot.state
        p.axis.axis_label_text_font_style = "bold"
        p.axis.minor_tick_line_color = None
        p.outline_line_width = 1
        p.outline_line_color = "black"

    @property
    def loadings(self):
        return self.data[self._l_columns]

    def hist(self, pcs=None, **kwargs):
        """Plot histogram of loadings"""
        if pcs is None:
            pcs = self.loadings.columns
        overlay = self.loadings[pcs].hvplot.hist(subplots=True, **kwargs)
        for sp in overlay:
            sp.opts(hooks=[self.plot_settings])
        return overlay

    def _top_loading_indices(self, pc, signif):
        data = np.sqrt(np.square(self.loadings[pc]).agg(np.sum, axis=1))
        upper = data.quantile(1 - signif / 2)
        lower = data.quantile(signif / 2)
        self._magnitude = data

        self._i = (data > upper) | (data < lower)

    @property
    def i(self):
        return self._i

    @property
    def index(self):
        return self._index

    def top_loadings(self, pc, signif=0.01):
        assert isinstance(pc, list), "pc argument must be list of components"
        assert all([x.startswith(self.model._label) for x in pc]), (
            "pc must have label %s; saw %s" % (self.model._label, pc)
        )
        self._top_loading_indices(pc=pc, signif=signif)
        retval = self.model.data.adata[self.i.values].join(
            self.data[self.i][pc]
        )
        retval["magnitude"] = self._magnitude[self.i]
        summary = self.model.coverage.summary()
        retval = retval.join(summary)
        if len(summary) > 1:
            pairs = itertools.combinations(summary, 2)
            for p in pairs:
                retval["-".join(p)] = summary[p[0]] - summary[p[1]]
        self._index = retval.index
        j = np.flip(np.argsort(np.square(retval[pc]).agg(np.sum, axis=1)))
        return retval.iloc[j]

    @property
    def cov(self):
        return self.model.coverage.subset(self.index)

    def reset(self):
        self._index = list(range(self.model.coverage.data.shape[0]))
        self._i = None


@dataclass
class PCALoadings(Loadings):
    model: PCA = field(default_factory=PCA)

    def __post_init__(self):
        assert isinstance(self.model, PCA), (
            "data must be of type %s; saw %s"
            % (
                PCA,
                type(self.model),
            )
        )
        super().__post_init__()


@dataclass
class PLSDALoadings(Loadings):
    model: PLSDA = field(default_factory=PLSDA)

    def __post_init__(self):
        assert isinstance(self.model, PLSDA), (
            "data must be of type %s; saw %s"
            % (
                PLSDA,
                type(self.model),
            )
        )
        super().__post_init__()


@dataclass
class Plotter:
    data: pd.DataFrame = field(default_factory=pd.DataFrame)
    width: int = 600
    height: int = 600
    samplesets: dict = None
    bgcolor: str = "#E5E5E5"
    line_color: str = "black"
    alpha: float = 1.0
    colormap: dict = None
    annotation: Annotation = field(default_factory=Annotation)
    backend_opts: dict = field(default_factory=dict)
    kwargs: dict = field(default_factory=dict)

    _data: pd.DataFrame = field(init=False, repr=False)

    def __post_init__(self):
        assert isinstance(self.data, pd.DataFrame), (
            "data must be of type %s; saw %s" % (pd.DataFrame, type(self.data))
        )
        if self.samplesets is None:
            self.samplesets = {"all": self.samples}
        revmap = {
            item: k
            for k, v in self.samplesets.items()
            for item in self.samplesets[k]
        }
        self.data["sampleset"] = [revmap[s] for s in self.samples]
        if not self.annotation.empty:
            self.data = self.data.join(self.annotation.data)
        self.backend_opts.update(
            {
                "plot.outline_line_width": 1,
                "plot.outline_line_color": "black",
                "plot.xaxis.minor_tick_line_color": None,
                "plot.yaxis.minor_tick_line_color": None,
                "plot.axis.axis_label_text_font_style": "bold",
                "plot.title.text_font_style": "bold",
                "plot.title.text_font_size": "28pt",
            }
        )

    def plot_settings(self, plot, element):
        p = plot.state
        p.axis.axis_label_text_font_style = "bold"
        p.axis.minor_tick_line_color = None
        p.outline_line_width = 1
        p.outline_line_color = "black"

    @property
    def data(self):  # noqa: F811
        return self._data

    @data.setter
    def data(self, data: pd.DataFrame):
        self._data = data

    @property
    def samples(self):
        if "sample" in self.data.index.names:
            return self.data.samples
        return self.data.columns.get_level_values("sample")


@dataclass
class CovPlotter(Plotter):
    cov: Coverage = field(default_factory=Coverage)
    xrotation: int = 45

    def __post_init__(self):
        assert isinstance(self.cov, Coverage), (
            "data must be of type %s; saw %s" % (Coverage, type(self.cov))
        )
        self.data = self.cov.data

    def boxplot(
        self,
        *,
        x: str = ["sampleset", "sample"],
        y: str = "value",
        color: str = "sampleset",
        group: str = None,
        label: str = None,
        **kwargs,
    ):
        """Make boxplot over samples"""
        self.kwargs.update(**kwargs)
        width = self.kwargs.pop("width", self.width)
        height = self.kwargs.pop("height", self.height)
        xrotation = self.kwargs.pop("xrotation", self.xrotation)
        group = self.kwargs.pop("group", self.cov.type)
        label = self.kwargs.pop("label", self.cov.feature)

        return hv.BoxWhisker(
            pd.melt(self.data), x, y, group=group, label=label
        ).opts(
            width=width,
            height=height,
            xrotation=xrotation,
            bgcolor=self.bgcolor,
            box_fill_color=color,
            cmap=self.colormap,
            **self.kwargs,
        )

    def hist(self, *, by: str = "sampleset", **kwargs):
        self.kwargs.update(**kwargs)
        group_label = self.kwargs.pop(
            "group_label", f"{self.cov.type} / {self.cov.feature}"
        )
        label = self.kwargs.pop("label", self.cov.feature)
        width = self.kwargs.pop("width", self.width)
        height = self.kwargs.pop("height", self.height)
        if self.colormap is not None:
            try:
                fill_color = [
                    self.colormap[x]
                    for x in set(self.data.columns.get_level_values(by))
                ]
                self.kwargs["fill_color"] = self.kwargs.get(
                    "fill_color", fill_color
                )
            except KeyError:
                logger.warning("Mismatch between colormap and grouping level")
        return pd.melt(self.data).hvplot.hist(
            by=by,
            width=width,
            height=height,
            group_label=group_label,
            label=label,
            **self.kwargs,
        )

    def dhist(self, *, pair: tuple, **kwargs):
        """Plot histogram of difference between two columns"""
        assert len(pair) == 2, "must supply two sampleset names for difference"
        group_label = self.kwargs.pop(
            "group_label", f"{self.cov.type} / {self.cov.feature}"
        )
        width = self.kwargs.pop("width", self.width)
        height = self.kwargs.pop("height", self.height)
        data = (
            self.data.groupby("sampleset", axis=1)
            .agg(np.mean)[pair]
            .diff(axis=1)
            .drop(pair[0], axis=1)
            .assign(diff=lambda x: -x.get(pair[1]))
            .drop(pair[1], axis=1)
        )
        return data.hvplot.hist(
            group_label=group_label, width=width, height=height, **kwargs
        ).opts(
            xlabel=f"diff ({pair[0]} - {pair[1]})", hooks=[self.plot_settings]
        )

    def _cluster(self, *, method, col_cluster=False):
        row_linkage = linkage(self.data, method=method)
        if col_cluster:
            col_linkage = linkage(self.data.T, method=method)
        else:
            col_linkage = None
        return row_linkage, col_linkage

    def heatmap(self, *, method="average", col_cluster=False, **kwargs):
        self.kwargs.update(**kwargs)
        self.kwargs["width"] = self.kwargs.pop("width", self.width)
        self.kwargs["height"] = self.kwargs.pop("height", self.height)
        row_linkage, col_linkage = self._cluster(
            method=method, col_cluster=col_cluster
        )
        irow = dendrogram(row_linkage, no_plot=True)["leaves"]
        irow.reverse()
        if col_cluster:
            icol = dendrogram(col_linkage, no_plot=True)["leaves"]
        else:
            icol = list(range(self.data.shape[1]))
        columns = self.data.columns.values[icol]
        data = pd.melt(
            self.data.iloc[irow][columns]
            .droplevel("sampleset", axis=1)
            .reset_index(level="gene"),
            id_vars=["gene"],
        )
        p = hv.HeatMap(data, ["sample", "gene"]).opts(
            hv.opts.HeatMap(
                labelled=[],
                tools=["hover"],
                cmap="viridis",
                colorbar=True,
                yaxis=None,
                xaxis=None,
                title=f"nrow x ncol={self.data.shape}",
                **self.kwargs,
            ),
        )
        return p

    def clustermap(
        self,
        annotate=False,
        *,
        method="average",
        col_cluster=False,
        cbar_pos=[1.0, 0.25, 0.03, 0.5],
        **kwargs,
    ):
        self.kwargs.update(**kwargs)
        if not annotate:
            self.kwargs.update(
                **{
                    "xticklabels": False,
                    "yticklabels": False,
                }
            )
            cbar_pos = None
        row_linkage, col_linkage = self._cluster(
            method=method, col_cluster=col_cluster
        )
        try:
            col_colors = [
                self.colormap[x]
                for x in self.data.columns.get_level_values("sampleset")
            ]
        except KeyError:
            logger.warning("Failed to set colormap")
            col_colors = None
        p = sns.clustermap(
            self.data,
            row_linkage=row_linkage,
            col_linkage=col_linkage,
            col_cluster=col_cluster,
            col_colors=col_colors,
            cbar_pos=cbar_pos,
            **self.kwargs,
        )
        if annotate:
            p.fig.suptitle(f"nrow x ncol={self.data.shape}")
        else:
            p.ax_heatmap.set_xlabel(None)
            p.ax_heatmap.set_ylabel(None)
            p.ax_row_dendrogram.set_visible(False)
        return p


@dataclass
class PCAPlotter(Plotter):
    model: PCA = field(default_factory=PCA)

    def __post_init__(self):
        assert isinstance(self.model, PCA), (
            "data must be of type %s; saw %s"
            % (
                PCA,
                type(self.model),
            )
        )
        self.data = self.model.pcs
        super().__post_init__()

    def scatter(
        self,
        title: str = None,
        *,
        data: pd.DataFrame = None,
        x: str = "PC1",
        y: str = "PC2",
        color: str = "sampleset",
        **kwargs,
    ):
        if data is None:
            data = self.data
        hover_cols = [x for x in data.columns if not x.startswith("PC")]
        xi = data.columns.values.tolist().index(x)
        yi = data.columns.values.tolist().index(y)
        xvar = self.model.explained_variance_ratio(xi)
        yvar = self.model.explained_variance_ratio(yi)
        xlab = f"{x} ({xvar * 100:.2f}%)"
        ylab = f"{y} ({yvar * 100:.2f}%)"
        if title is None:
            title = (
                f"{self.model.type}/{self.model.feature} "
                f"({self.model.n_features} features, "
                f"{self.model.n_genes} genes)"
            )
            if self.model.randomized:
                title = f"RND {title}"
        self.kwargs.update(**kwargs)
        self.kwargs["hover_cols"] = self.kwargs.get("hover_cols", hover_cols)
        return data.hvplot.scatter(
            x=x,
            y=y,
            title=title,
            color=color,
            xlabel=xlab,
            ylabel=ylab,
            cmap=self.colormap,
            bgcolor=self.bgcolor,
            alpha=self.alpha,
            line_color=self.line_color,
            **self.kwargs,
        ).opts(
            backend_opts=self.backend_opts,
        )

    @property
    def samples(self):
        return self.model.samples


@dataclass
class PLSDAPlotter(Plotter):
    model: PLSDA = field(default_factory=PLSDA)

    def __post_init__(self):
        assert isinstance(self.model, PLSDA), (
            "data must be of type %s; saw %s"
            % (
                PLSDA,
                type(self.model),
            )
        )
        self.data = self.model.pcs
        super().__post_init__()

    def scatter(
        self,
        title: str = None,
        *,
        data: pd.DataFrame = None,
        x: str = "PLS1",
        y: str = "PLS2",
        color: str = "sampleset",
        **kwargs,
    ):
        if data is None:
            data = self.data
        hover_cols = [x for x in data.columns if not x.startswith("PC")]
        xlab = f"{x}"
        ylab = f"{y}"
        if title is None:
            title = (
                f"{self.model.type}/{self.model.feature} "
                f"({self.model.n_features} features, "
                f"{self.model.n_genes} genes)"
            )
            if self.model.randomized:
                title = f"RND {title}"
        self.kwargs.update(**kwargs)
        self.kwargs["hover_cols"] = self.kwargs.get("hover_cols", hover_cols)
        return data.hvplot.scatter(
            x=x,
            y=y,
            title=title,
            color=color,
            xlabel=xlab,
            ylabel=ylab,
            cmap=self.colormap,
            bgcolor=self.bgcolor,
            alpha=self.alpha,
            line_color=self.line_color,
            **self.kwargs,
        ).opts(
            backend_opts=self.backend_opts,
        )

    @property
    def samples(self):
        return self.model.samples
