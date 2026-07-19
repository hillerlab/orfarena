#!/usr/bin/env python3
"""Format TransDecoder and GeneMarkS-T predictions for ORF Arena."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Final, TypeAlias, cast

import numpy as np
import numpy.typing as npt

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/hillerlab/orfarena"
__version__ = "0.0.1"

# This formatter only uses pandas' core tabular operations. Prevent pandas from
# probing optional compiled accelerators during import: an incompatible optional
# pyarrow/bottleneck install should not make this TSV-only executable unusable.
_disabled_optional_modules: list[str] = []
for _module_name in ("pyarrow", "bottleneck"):
    if _module_name not in sys.modules:
        sys.modules[_module_name] = cast(ModuleType, None)
        _disabled_optional_modules.append(_module_name)
try:
    import pandas as pd
finally:
    for _module_name in _disabled_optional_modules:
        sys.modules.pop(_module_name, None)


DataFrame: TypeAlias = pd.DataFrame
Series: TypeAlias = pd.Series
FloatArray: TypeAlias = npt.NDArray[np.float64]
Column: TypeAlias = int | str

COL_CHROM: Final[int] = 0
COL_START: Final[int] = 1
COL_NAME: Final[int] = 3
COL_STRAND: Final[int] = 5
COL_THICK_START: Final[int] = 6
COL_THICK_END: Final[int] = 7
COL_BLOCK_SIZES: Final[int] = 10
COL_BLOCK_STARTS: Final[int] = 11
COL_QUERY_START: Final[int] = 12
COL_QUERY_END: Final[int] = 13
COL_SCORE: Final[str] = "score"
COL_KEY: Final[str] = "key"
COL_PROBABILITY: Final[str] = "coding_prob"
BOUNDARY: float | None = None


@dataclass(frozen=True, slots=True)
class CommonOptions:
    """Arguments shared by all formatters."""

    outdir: Path
    prefix: str
    ground: Path


@dataclass(frozen=True, slots=True)
class TransDecoderOptions(CommonOptions):
    """Arguments for the TransDecoder formatter."""

    td_bed: Path
    mapped_cds: Path


@dataclass(frozen=True, slots=True)
class GeneMarkOptions(CommonOptions):
    """Arguments for the GeneMarkS-T formatter."""

    gff: Path
    mapped_cds: Path


@dataclass(frozen=True, slots=True)
class Td2Options(CommonOptions):
    """Arguments for the TD2 formatter."""

    mapped_cds: Path
    td2: Path | None


@dataclass(frozen=True, slots=True)
class PsauronOptions(CommonOptions):
    """Arguments for the PSAURON formatter."""

    psauron: Path
    noncoding: bool


def _integers(value: object) -> list[int]:
    """Parse a comma-separated BED integer list."""
    return [int(item) for item in str(value).rstrip(",").split(",") if item]


def cds_length(row: Series) -> int:
    """Return the coding length after intersecting BED blocks with the CDS."""
    chrom_start = int(row[COL_START])
    thick_start = int(row[COL_THICK_START])
    thick_end = int(row[COL_THICK_END])
    length = 0

    block_sizes = _integers(row[COL_BLOCK_SIZES])
    block_starts = _integers(row[COL_BLOCK_STARTS])
    if len(block_sizes) != len(block_starts):
        raise ValueError("BED blockSizes and blockStarts have different lengths")

    for block_size, block_offset in zip(block_sizes, block_starts, strict=True):
        block_start = chrom_start + block_offset
        block_end = block_start + block_size
        overlap_start = max(block_start, thick_start)
        overlap_end = min(block_end, thick_end)
        if overlap_end > overlap_start:
            length += overlap_end - overlap_start

    return length


def transcript_length(row: Series) -> int:
    """Return the full spliced transcript length in nucleotides."""
    return sum(_integers(row[COL_BLOCK_SIZES]))


def add_features(frame: DataFrame) -> DataFrame:
    """Add coding-length and length-normalized score columns."""
    result = frame.copy()
    result["cds_len"] = result.apply(cds_length, axis=1)
    result["tx_len"] = result.apply(transcript_length, axis=1)
    safe_cds_length = result["cds_len"].clip(lower=1)
    result["score_per_nt"] = result[COL_SCORE] / safe_cds_length
    result["score_per_codon"] = result[COL_SCORE] / (safe_cds_length / 3.0)
    return result


def find_boundary(values: Series | FloatArray, fallback: float = 0.0) -> float:
    """Find the crossover between two Gaussian components in a score series."""
    scores = np.asarray(values, dtype=np.float64)
    scores = scores[np.isfinite(scores)]
    if scores.size < 2 or np.unique(scores).size < 2:
        return fallback

    # Fit a small, deterministic 1-D Gaussian mixture here instead of taking a
    # compiled scikit-learn dependency for this one calculation.
    means = np.quantile(scores, [0.25, 0.75]).astype(np.float64)
    overall_variance = max(float(np.var(scores)), 1e-12)
    minimum_variance = max(overall_variance * 1e-6, 1e-12)
    variances = np.full(2, overall_variance, dtype=np.float64)
    weights = np.full(2, 0.5, dtype=np.float64)

    for _ in range(200):
        old_means = means.copy()
        log_probabilities = np.empty((scores.size, 2), dtype=np.float64)
        for component in range(2):
            log_probabilities[:, component] = (
                np.log(weights[component])
                - 0.5 * np.log(2.0 * np.pi * variances[component])
                - (scores - means[component]) ** 2 / (2.0 * variances[component])
            )
        row_maxima = log_probabilities.max(axis=1, keepdims=True)
        probabilities = np.exp(log_probabilities - row_maxima)
        responsibilities = probabilities / probabilities.sum(axis=1, keepdims=True)
        component_totals = responsibilities.sum(axis=0)
        if np.any(component_totals <= 1e-12):
            return fallback

        weights = component_totals / scores.size
        means = (responsibilities * scores[:, np.newaxis]).sum(
            axis=0
        ) / component_totals
        for component in range(2):
            residuals = scores - means[component]
            variance = float(
                np.dot(responsibilities[:, component], residuals**2)
                / component_totals[component]
            )
            variances[component] = max(variance, minimum_variance)
        if np.allclose(means, old_means, rtol=1e-8, atol=1e-10):
            break

    component_order = np.argsort(means)
    low_component, high_component = map(int, component_order)
    grid = np.linspace(means[low_component], means[high_component], 2_001)
    weighted_densities = np.empty((grid.size, 2), dtype=np.float64)
    for component in range(2):
        weighted_densities[:, component] = (
            weights[component]
            * np.exp(-((grid - means[component]) ** 2) / (2.0 * variances[component]))
            / np.sqrt(2.0 * np.pi * variances[component])
        )
    difference = (
        weighted_densities[:, low_component] - weighted_densities[:, high_component]
    )
    crossings = np.flatnonzero(np.diff(np.sign(difference)))
    return fallback if crossings.size == 0 else float(grid[crossings[0]])


def to_unit(
    normalized_scores: Series | FloatArray,
    boundary: float,
    steepness: float | None = None,
) -> FloatArray:
    """Map normalized scores to [0, 1], with 0.5 at ``boundary``."""
    scores = np.asarray(normalized_scores, dtype=np.float64)
    if steepness is None:
        standard_deviation = float(np.nanstd(scores))
        steepness = 1.0 / standard_deviation if standard_deviation > 0 else 1.0
    exponent = np.clip(-steepness * (scores - boundary), -709.0, 709.0)
    return cast(FloatArray, 1.0 / (1.0 + np.exp(exponent)))


def _score_from_attribute(value: object) -> float:
    """Extract the final ``name=value`` score from a predictor attribute."""
    final_attribute = str(value).rsplit(",", maxsplit=1)[-1]
    if "=" not in final_attribute:
        raise ValueError(f"could not extract a score from {value!r}")
    return float(final_attribute.rsplit("=", maxsplit=1)[-1])


def _prediction_key(
    frame: DataFrame,
    chrom: Column,
    thick_start: Column,
    thick_end: Column,
    name: Column,
) -> Series:
    """Build the stable key used to join predictions to ground records."""
    return (
        frame[chrom].astype(str)
        + ":"
        + frame[thick_start].astype(str)
        + "-"
        + frame[thick_end].astype(str)
        + "@"
        + frame[name].astype(str)
    )


def transdecoder_df_builder(
    mapped_cds: DataFrame,
    transdecoder_bed: DataFrame,
) -> DataFrame:
    """Attach TransDecoder scores to mapped BED12 CDS records."""
    mapped = mapped_cds.copy()
    scores = transdecoder_bed.copy()
    scores[COL_SCORE] = scores[COL_NAME].map(_score_from_attribute)
    scores[COL_KEY] = (
        scores[COL_CHROM].astype(str)
        + ":"
        + scores[COL_THICK_START].astype(str)
        + "-"
        + scores[COL_THICK_END].astype(str)
    )
    mapped[COL_KEY] = (
        mapped[COL_NAME].astype(str)
        + ":"
        + mapped[COL_QUERY_START].astype(str)
        + "-"
        + mapped[COL_QUERY_END].astype(str)
    )
    result = mapped.merge(
        scores[[COL_KEY, COL_SCORE]],
        on=COL_KEY,
        how="left",
        validate="many_to_one",
    )
    result[COL_SCORE] = result[COL_SCORE].fillna(0.0)
    return result.drop(columns=[COL_KEY, COL_QUERY_START, COL_QUERY_END])


def genemark_df_builder(
    genemark_gff: DataFrame,
    mapped_cds: DataFrame,
) -> DataFrame:
    """Attach GeneMarkS-T scores to mapped BED12 CDS records."""
    scores = genemark_gff.copy()
    mapped = mapped_cds.copy()
    scores[COL_SCORE] = scores[8].map(_score_from_attribute)
    scores = (
        scores[[COL_CHROM, COL_SCORE]]
        .sort_values(COL_SCORE, ascending=False)
        .drop_duplicates(COL_CHROM)
    )
    result = mapped.merge(
        scores,
        left_on=COL_NAME,
        right_on=COL_CHROM,
        how="left",
        validate="many_to_one",
    )
    result[COL_SCORE] = result[COL_SCORE].fillna(0.0)
    result = result.drop(columns=[f"{COL_CHROM}_y"])
    result = result.rename(columns={f"{COL_CHROM}_x": COL_CHROM})

    strands = set(result[COL_STRAND].dropna().astype(str))
    invalid_strands = strands.difference({"+", "-"})
    if invalid_strands:
        invalid = ", ".join(sorted(invalid_strands))
        raise ValueError(f"mapped CDS contains invalid strand value(s): {invalid}")

    plus = result[result[COL_STRAND] == "+"].copy()
    plus[COL_THICK_START] = plus[COL_THICK_START] - 1
    minus = result[result[COL_STRAND] == "-"].copy()
    minus[COL_THICK_END] = minus[COL_THICK_END] + 1
    return (
        pd.concat([plus, minus], ignore_index=True)
        .sort_values(COL_SCORE, ascending=False)
        .drop_duplicates(COL_NAME)
    )


def _td2_orf_coordinates(value: object) -> tuple[int, int]:
    """Extract transcript-relative ORF coordinates from a TD2 BED name."""
    name = str(value)
    try:
        coordinates = name.rsplit(":", maxsplit=1)[-1].split("(", maxsplit=1)[0]
        fields = coordinates.split("-")
        return int(fields[0]), int(fields[-1])
    except (IndexError, ValueError) as error:
        raise ValueError(
            f"could not extract TD2 ORF coordinates from {name!r}"
        ) from error


def td2_df_builder(
    mapped_cds: DataFrame,
    td2_scores: DataFrame | None,
) -> DataFrame:
    """Attach scaled TD2 scores to mapped BED12 CDS records."""
    mapped = mapped_cds.copy()

    if td2_scores is None:
        if mapped.shape[1] < 15:
            raise ValueError(
                "TD2 formatting requires --td2 because a standard 14-column "
                "mapped CDS BED does not contain TD2 scores"
            )
        raw_scores = pd.to_numeric(mapped.iloc[:, 14], errors="raise")
        mapped = mapped.iloc[:, :14].copy()
        mapped[COL_SCORE] = raw_scores.to_numpy()
    else:
        scores = td2_scores.copy()
        coordinates = scores[COL_NAME].map(_td2_orf_coordinates)
        scores["orf_start"] = coordinates.map(lambda value: value[0])
        scores["orf_end"] = coordinates.map(lambda value: value[1])
        scores["inner_key"] = (
            scores[COL_CHROM].astype(str)
            + ":"
            + scores["orf_start"].astype(str)
            + "-"
            + scores["orf_end"].astype(str)
        )
        scores[COL_SCORE] = pd.to_numeric(scores[4], errors="raise")

        mapped["inner_key"] = (
            mapped[COL_NAME].astype(str)
            + ":"
            + mapped[COL_QUERY_START].astype(str)
            + "-"
            + mapped[COL_QUERY_END].astype(str)
        )
        mapped = mapped.merge(
            scores[["inner_key", COL_SCORE]],
            on="inner_key",
            how="left",
        )

    mapped[COL_SCORE] = mapped[COL_SCORE] / 1_000.0

    strands = set(mapped[COL_STRAND].dropna().astype(str))
    invalid_strands = strands.difference({"+", "-"})
    if invalid_strands:
        invalid = ", ".join(sorted(invalid_strands))
        raise ValueError(f"mapped CDS contains invalid strand value(s): {invalid}")

    plus = mapped[mapped[COL_STRAND] == "+"].copy()
    plus[COL_THICK_START] = plus[COL_THICK_START] - 1
    minus = mapped[mapped[COL_STRAND] == "-"].copy()
    minus[COL_THICK_END] = minus[COL_THICK_END] + 1
    result = pd.concat([plus, minus], ignore_index=True)
    result[COL_KEY] = _prediction_key(
        result,
        COL_CHROM,
        COL_THICK_START,
        COL_THICK_END,
        COL_NAME,
    )
    return result


def _read_table(path: Path, minimum_columns: int, description: str) -> DataFrame:
    """Read a headerless tab-separated predictor table and validate its width."""
    with path.open(encoding="utf-8") as stream:
        metadata_rows = [
            line_number
            for line_number, line in enumerate(stream)
            if line.startswith(("track ", "browser "))
        ]
    try:
        frame = pd.read_csv(
            path,
            sep="\t",
            header=None,
            comment="#",
            skiprows=metadata_rows,
        )
    except pd.errors.EmptyDataError as error:
        raise ValueError(f"{description} is empty: {path}") from error
    if frame.shape[1] < minimum_columns:
        raise ValueError(
            f"{description} must contain at least {minimum_columns} tab-separated "
            f"columns; found {frame.shape[1]} in {path}"
        )
    return frame


def _read_ground(path: Path) -> DataFrame:
    """Read ground records and ensure that they contain a merge key."""
    with path.open(encoding="utf-8") as stream:
        first_data_line = next(
            (
                line.rstrip("\n")
                for line in stream
                if line.strip() and not line.startswith("#")
            ),
            None,
        )
    if first_data_line is None:
        raise ValueError(f"ground file is empty: {path}")

    fields = first_data_line.split("\t")
    has_header = COL_KEY in fields
    ground = pd.read_csv(
        path,
        sep="\t",
        header=0 if has_header else None,
        comment="#",
    )
    if COL_KEY in ground.columns:
        ground[COL_KEY] = ground[COL_KEY].astype(str)
        return ground
    if ground.shape[1] == 1:
        ground = ground.rename(columns={ground.columns[0]: COL_KEY})
        ground[COL_KEY] = ground[COL_KEY].astype(str)
        return ground
    if ground.shape[1] < 12:
        raise ValueError(
            "ground must be BED12, a one-column key file, or a table with a "
            f"'key' header; found {ground.shape[1]} columns in {path}"
        )
    ground[COL_KEY] = _prediction_key(
        ground,
        COL_CHROM,
        COL_THICK_START,
        COL_THICK_END,
        COL_NAME,
    )
    return ground


def _score_predictions(predictions: DataFrame) -> tuple[DataFrame, float]:
    """Length-normalize scores and calculate coding probabilities."""
    features = add_features(predictions)
    boundary = find_boundary(features["score_per_nt"], fallback=0.0)

    if BOUNDARY is not None:
        boundary = BOUNDARY

    features[COL_PROBABILITY] = to_unit(
        features["score_per_nt"],
        boundary=boundary,
    )
    return features.sort_values(COL_SCORE), boundary


def _save_distribution(probabilities: Series, output_path: Path) -> None:
    """Save the coding-probability distribution as a PNG histogram."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(figsize=(8.0, 5.0), constrained_layout=True)
    axes.hist(probabilities.to_numpy(dtype=float), bins=100)
    axes.axvline(0.5, color="black", linestyle="--", linewidth=1.0)
    axes.set(xlabel="Coding probability", ylabel="Count", title="Score distribution")
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def _save_results(
    ground: DataFrame,
    predictions: DataFrame,
    output_path: Path,
) -> DataFrame:
    """Merge probabilities onto ground records and save a headerless BED file."""
    probabilities = (
        predictions[[COL_KEY, COL_PROBABILITY]]
        .sort_values(COL_PROBABILITY, ascending=False)
        .drop_duplicates(COL_KEY)
    )
    results = ground.merge(
        probabilities,
        on=COL_KEY,
        how="left",
        validate="many_to_one",
    )
    results[COL_PROBABILITY] = results[COL_PROBABILITY].fillna(0.0)
    results.to_csv(output_path, sep="\t", index=False, header=False)
    return results


def _write_outputs(
    predictions: DataFrame,
    ground: DataFrame,
    options: CommonOptions,
) -> DataFrame:
    """Score predictions and write both required output files."""
    scored, boundary = _score_predictions(predictions)
    scored[COL_KEY] = _prediction_key(
        scored,
        COL_CHROM,
        COL_THICK_START,
        COL_THICK_END,
        COL_NAME,
    )

    distribution_path = options.outdir / f"{options.prefix}.distribution.png"
    results_path = options.outdir / f"{options.prefix}.results.bed"
    # _save_distribution(scored[COL_PROBABILITY], distribution_path)
    results = _save_results(ground, scored, results_path)

    positives = int((scored[COL_PROBABILITY] >= 0.5).sum())
    negatives = int((scored[COL_PROBABILITY] < 0.5).sum())
    print(f"boundary: {boundary}")
    print(f"positives: {positives}")
    print(f"negatives: {negatives}")
    print(f"distribution: {distribution_path}")
    print(f"results: {results_path}")
    return results


def run_transdecoder(options: TransDecoderOptions) -> DataFrame:
    """Run the TransDecoder formatting workflow."""
    transdecoder_bed = _read_table(options.td_bed, 8, "TransDecoder BED")
    mapped_cds = _read_table(options.mapped_cds, 14, "mapped CDS BED")
    ground = _read_ground(options.ground)
    predictions = transdecoder_df_builder(mapped_cds, transdecoder_bed)
    return _write_outputs(predictions, ground, options)


def run_genemark(options: GeneMarkOptions) -> DataFrame:
    """Run the GeneMarkS-T formatting workflow."""
    gff = _read_table(options.gff, 9, "GeneMark GFF")
    mapped_cds = _read_table(options.mapped_cds, 14, "mapped CDS BED")
    ground = _read_ground(options.ground)
    predictions = genemark_df_builder(gff, mapped_cds)
    return _write_outputs(predictions, ground, options)


def run_td2(options: Td2Options) -> DataFrame:
    """Run the TD2 formatting workflow."""
    mapped_cds = _read_table(options.mapped_cds, 14, "mapped CDS BED")
    td2_scores = (
        _read_table(options.td2, 5, "TD2 BED") if options.td2 is not None else None
    )
    ground = _read_ground(options.ground)
    predictions = td2_df_builder(mapped_cds, td2_scores)
    results = (
        ground.merge(
            predictions[[COL_KEY, COL_SCORE]],
            on=COL_KEY,
            how="left",
        )
        .fillna(0.0)
        .sort_values(COL_SCORE)
    )
    results_path = options.outdir / f"{options.prefix}.results.bed"
    results.to_csv(results_path, sep="\t", index=False, header=False)
    print(f"results: {results_path}")
    return results


def run_psauron(options: PsauronOptions) -> DataFrame:
    """Run the PSAURON formatting workflow."""
    ground = _read_table(options.ground, 4, "ground BED")
    try:
        scores = pd.read_csv(options.psauron, comment="#", skiprows=3)
    except pd.errors.EmptyDataError as error:
        raise ValueError(f"PSAURON CSV is empty: {options.psauron}") from error

    ground["key"] = (
        ground[COL_CHROM].astype(str)
        + ":"
        + ground[COL_THICK_START].astype(str)
        + "-"
        + ground[COL_THICK_END].astype(str)
        + "@"
        + ground[COL_NAME]
    )

    required_columns = {"description", "in_frame_score"}
    missing_columns = required_columns.difference(scores.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"PSAURON CSV is missing required column(s): {missing}")

    scores["in_frame_score"] = pd.to_numeric(
        scores["in_frame_score"],
        errors="raise",
    )

    results = (
        ground.merge(
            scores[["description", "in_frame_score"]],
            left_on=COL_NAME,
            right_on="description",
            how="left",
        )
        .drop(columns=["description"])
        .fillna(0.0)
    )

    if options.noncoding:
        results = (
            results.sort_values("in_frame_score", ascending=False)
            .drop_duplicates(COL_NAME)
            .reset_index(drop=True)
        )



    results_path = options.outdir / f"{options.prefix}.results.bed"
    results.to_csv(results_path, sep="\t", index=False, header=False)
    print(f"results: {results_path}")
    return results


def _add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """Add common options while allowing them on either side of a subcommand."""
    parser.add_argument(
        "--outdir",
        type=Path,
        default=argparse.SUPPRESS,
        metavar="PATH",
        help="directory for generated files",
    )
    parser.add_argument(
        "--prefix",
        default=argparse.SUPPRESS,
        metavar="STR",
        help="prefix for generated file names",
    )
    parser.add_argument(
        "-g",
        "--ground",
        type=Path,
        default=argparse.SUPPRESS,
        metavar="PATH",
        help="ground-truth BED or key table",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="show program's version number and exit",
    )


def build_parser() -> argparse.ArgumentParser:
    """Construct the command-line parser."""
    parser = argparse.ArgumentParser(
        description="Format coding-potential predictions for ORF Arena.",
    )
    _add_common_arguments(parser)
    subparsers = parser.add_subparsers(dest="command", required=True)

    transdecoder = subparsers.add_parser(
        "transdecoder",
        help="format TransDecoder predictions",
    )
    _add_common_arguments(transdecoder)
    transdecoder.add_argument(
        "--td-bed",
        required=True,
        type=Path,
        metavar="PATH",
        help="TransDecoder BED file",
    )
    transdecoder.add_argument(
        "--mapped-cds",
        required=True,
        type=Path,
        metavar="PATH",
        help="mapped CDS BED file",
    )

    genemark = subparsers.add_parser(
        "genemark",
        help="format GeneMarkS-T predictions",
    )
    _add_common_arguments(genemark)
    genemark.add_argument(
        "--gff",
        required=True,
        type=Path,
        metavar="PATH",
        help="GeneMark GFF file",
    )
    genemark.add_argument(
        "--mapped-cds",
        required=True,
        type=Path,
        metavar="PATH",
        help="mapped CDS BED file",
    )

    td2 = subparsers.add_parser(
        "td2",
        help="format TD2 predictions",
    )
    _add_common_arguments(td2)
    td2.add_argument(
        "--mapped-cds",
        "--mapped_cds",
        dest="mapped_cds",
        required=True,
        type=Path,
        metavar="PATH",
        help="mapped CDS BED file",
    )
    td2.add_argument(
        "--td2",
        type=Path,
        metavar="PATH",
        help="raw TD2 BED containing transcript-relative coordinates and scores",
    )

    psauron = subparsers.add_parser(
        "psauron",
        help="format PSAURON predictions",
    )
    _add_common_arguments(psauron)
    psauron.add_argument(
        "--psauron",
        required=True,
        type=Path,
        metavar="PATH",
        help="PSAURON score CSV",
    )
    psauron.add_argument(
        "--noncoding",
        action="store_true",
        help="sort scores descending and keep the best row per transcript",
    )
    return parser


def _common_options(namespace: argparse.Namespace) -> CommonOptions:
    """Validate and extract common CLI options."""
    missing = [
        option
        for option in ("outdir", "prefix", "ground")
        if not hasattr(namespace, option)
    ]
    if missing:
        rendered = ", ".join(f"--{option}" for option in missing)
        raise ValueError(f"the following arguments are required: {rendered}")

    outdir = cast(Path, namespace.outdir)
    prefix = cast(str, namespace.prefix)
    ground = cast(Path, namespace.ground)
    if not prefix:
        raise ValueError("--prefix cannot be empty")
    for path, label in ((ground, "--ground"),):
        if not path.is_file():
            raise ValueError(f"{label} is not a readable file: {path}")
    outdir.mkdir(parents=True, exist_ok=True)
    return CommonOptions(outdir=outdir, prefix=prefix, ground=ground)


def _require_file(path: Path, option: str) -> Path:
    """Ensure that a CLI input path is a file."""
    if not path.is_file():
        raise ValueError(f"{option} is not a readable file: {path}")
    return path


def main(argv: Sequence[str] | None = None) -> int:
    """Run the selected formatter."""
    parser = build_parser()
    namespace = parser.parse_args(argv)
    try:
        common = _common_options(namespace)
        if namespace.command == "transdecoder":
            options = TransDecoderOptions(
                outdir=common.outdir,
                prefix=common.prefix,
                ground=common.ground,
                td_bed=_require_file(namespace.td_bed, "--td-bed"),
                mapped_cds=_require_file(namespace.mapped_cds, "--mapped-cds"),
            )
            run_transdecoder(options)
        elif namespace.command == "genemark":
            options = GeneMarkOptions(
                outdir=common.outdir,
                prefix=common.prefix,
                ground=common.ground,
                gff=_require_file(namespace.gff, "--gff"),
                mapped_cds=_require_file(namespace.mapped_cds, "--mapped-cds"),
            )
            run_genemark(options)
        elif namespace.command == "td2":
            td2_path = (
                _require_file(namespace.td2, "--td2")
                if namespace.td2 is not None
                else None
            )
            options = Td2Options(
                outdir=common.outdir,
                prefix=common.prefix,
                ground=common.ground,
                mapped_cds=_require_file(namespace.mapped_cds, "--mapped-cds"),
                td2=td2_path,
            )
            run_td2(options)
        elif namespace.command == "psauron":
            options = PsauronOptions(
                outdir=common.outdir,
                prefix=common.prefix,
                ground=common.ground,
                psauron=_require_file(namespace.psauron, "--psauron"),
                noncoding=namespace.noncoding,
            )
            run_psauron(options)
        else:  # pragma: no cover - argparse enforces the command choices.
            parser.error(f"unknown command: {namespace.command}")
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(status=2, message=f"format.py: error: {error}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
