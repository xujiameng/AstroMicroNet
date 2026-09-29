# Validation

Run commands from the repository root. NumPy and SciPy suffice for input regeneration and AstroMicroNet replay.

| Command | What it checks |
|---|---|
| `python validation/check_synthetic_generators.py` | Regenerates 80 weighted networks, 500 null networks and 40 waveform inputs using their archived seeds |
| `python validation/reproduce.py` | Recomputes V01_baseline; compares exact final labels and support values |
| `python validation/reproduce.py --bank null` | Recomputes NULL0001 |
| `python validation/reproduce.py --bank waveform` | Recomputes the first waveform input under all three connectivity definitions |

Add `--all` to replay every input in a bank. These scripts use the archived settings and never retune parameters. The synthesis code was extracted from `review_fair_tuning.py`, `review_expand_null.py` and `review_waveform_benchmark.py` in the September 2026 analysis archive. Paths were adapted for this repository.

## BE comparison

The complete numerical BE checks, original/checked kernels, inputs and expected outputs are in `supplementary/supporting_data/BE_checks/`. To replay the displayed fifth-percentile binary setting:

```text
python -m pip install -r supplementary/supporting_data/BE_checks/requirements.txt
python supplementary/supporting_data/BE_checks/run_crosscheck.py --input supplementary/supporting_data/BE_checks/inputs/V01_baseline.npz --output outputs/be_example
```

Omit `--input` to run all 80 matrices. Optional `--implementation independent-weighted` runs the separate weighted diagnostic. The displayed binary BE choice is exploratory, is not the best-performing BE setting, and does not establish a universal method ranking. Full alternative diagnostic results are retained.

`method_expected/` preserves the earlier comparison outputs: its `BE_labels` use the earlier adaptive threshold. **They are not the updated S8 BE labels.** The updated labels are in `supplementary/supporting_data/BE_checks/expected_q05_seed0/labels/`; updated plotted scores are in `supplementary/data/method_comparison.csv`.

Rombach labels/scores, waveform comparator results and development-grid values are archived for inspection. The original cpnet 0.0.21 comparator sources (BE, Rombach and helpers) and their bundled license/metadata are in `third_party/cpnet/`. This is a source archive, not a complete pip-installable cpnet package. The default replay script reruns AstroMicroNet; it does not rerun those other comparators or repeat model selection.

The BE reference environment is Python 3.12, NumPy 2.3.5, SciPy 1.15.3 and Numba 0.63.1. Use a separate environment for its pinned requirements. The core software and plotting checks were run in Python 3.13 with NumPy 2.2.5 and SciPy 1.15.3. Package compatibility and exact random solutions can differ across environments.

The paper's raw-imaging workflow and omitted real-cell GUI regression fixtures require separate authorized inputs. This selected release does not claim full raw-data reproduction.
