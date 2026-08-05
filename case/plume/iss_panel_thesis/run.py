"""Run the ISS-representative solar-panel plume-impingement studies.

Each study is a YAML file in ``study/`` layered on this case's
``config.ini``; running one is the ordinary package-level call

    TradeStudy.from_config(<yaml>).run()

so this script only picks the file, optionally redirects the output, and
prints where the artifacts landed.

Usage
-----
    python run.py                        # list the available studies
    python run.py baseline-simplified    # centered source, Simplified model
    python run.py baseline-full-cai      # centered source, full Cai model
    python run.py sweep                  # 3 distances x 5 u offsets
    python run.py all                    # every study, in that order

    python run.py sweep --output-dir /tmp/scratch   # run elsewhere
    python run.py sweep --no-plots                  # skip the figures

Nothing here handles DSMC: these studies produce ANALYTICAL PyRPOD datasets
only. Comparing them with externally generated DSMC results is a separate
workflow.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
# Allow `python run.py` from anywhere without installing the package.
sys.path.insert(0, str(CASE_DIR.parents[2]))

from pyrpod.mdao.TradeStudy import TradeStudy  # noqa: E402

#: Study name -> configuration file, in the order `all` runs them.
STUDIES: dict[str, str] = {
    "baseline-simplified": "iss_panel_baseline_simplified.yaml",
    "baseline-full-cai": "iss_panel_baseline_full_cai.yaml",
    "sweep": "iss_panel_offset_distance_sweep.yaml",
}


def run_study(name: str, output_dir: str | None = None,
              plots: bool = True) -> None:
    """Run one named study and report its artifacts."""
    config_path = CASE_DIR / "study" / STUDIES[name]
    print(f"\n=== {name}: {config_path.name} ===")

    study = TradeStudy.from_config(config_path, output_dir=output_dir)
    results = study.run()

    print(f"cases        : {len(results)} records "
          f"({study.study_config.n_cases} swept cases)")
    print(f"summary CSV  : {results.summary_csv_path}")
    print(f"metadata JSON: {results.metadata_path}")

    first = results.cases[0]
    print(f"model        : {first.plume_model}")
    print(f"normal force : {first.normal_force:.4g} N "
          "(positive = pressed into the panel)")
    print(f"peak pressure: {first.max_pressure:.4g} Pa")
    if first.knudsen_number is not None:
        print(f"Kn (metadata): {first.knudsen_number:.4g} "
              f"[{first.knudsen_definition}]")
    if first.surface_distribution_path:
        print(f"distribution : {first.surface_distribution_path}")
    if first.vtk_path:
        print(f"VTK          : {first.vtk_path}")

    if plots:
        written = study.plot()
        print(f"plots        : {len(written)} figures")
        for path in written[:6]:
            print(f"               {path}")
        if len(written) > 6:
            print(f"               ... and {len(written) - 6} more")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("study", nargs="?", default=None,
                        choices=[*STUDIES, "all"],
                        help="which study to run; omit to list them")
    parser.add_argument("--output-dir", default=None,
                        help="override the configured output directory")
    parser.add_argument("--no-plots", action="store_true",
                        help="skip figure generation")
    parser.add_argument("--verbose", action="store_true",
                        help="log study progress at INFO level")
    args = parser.parse_args(argv)

    if args.verbose:
        logging.basicConfig(level=logging.INFO,
                            format="%(levelname)s %(name)s: %(message)s")

    if args.study is None:
        print(__doc__)
        print("Available studies:")
        for name, filename in STUDIES.items():
            print(f"  {name:20s} study/{filename}")
        return 0

    names = list(STUDIES) if args.study == "all" else [args.study]
    for name in names:
        # With several studies and one --output-dir, keep them apart.
        output_dir = (f"{args.output_dir}/{name}"
                      if args.output_dir and len(names) > 1
                      else args.output_dir)
        run_study(name, output_dir=output_dir, plots=not args.no_plots)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
