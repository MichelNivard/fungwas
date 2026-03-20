from pathlib import Path

import numpy as np

from fungwas_stage1.cli import (
    build_rif_residuals,
    load_loco_predictions,
    resolve_loco_files,
)
from fungwas_stage1 import core


def test_resolve_loco_files_single_file(tmp_path: Path):
    loco_file = tmp_path / "bmi.loco"
    loco_file.write_text("FID_IID sample1\n22 0.0\n")

    mapping = resolve_loco_files(str(loco_file), ["bmi"])

    assert mapping == {"bmi": str(loco_file)}


def test_load_loco_predictions_aligns_fid_iid_aliases(tmp_path: Path):
    loco_file = tmp_path / "bmi.loco"
    loco_file.write_text(
        "FID_IID IEU2_IEU2 IEU1_IEU1\n"
        "21 0.2 0.1\n"
        "22 0.4 0.3\n"
    )

    preds, valid_mask = load_loco_predictions(
        str(loco_file),
        ["bmi"],
        [["IEU1", "IEU1_IEU1"], ["IEU2", "IEU2_IEU2"]],
    )

    np.testing.assert_array_equal(valid_mask, np.array([True, True]))
    np.testing.assert_allclose(preds["bmi"]["21"], np.array([0.1, 0.2]))
    np.testing.assert_allclose(preds["bmi"]["22"], np.array([0.3, 0.4]))


def test_load_loco_predictions_drops_samples_missing_from_loco(tmp_path: Path):
    loco_file = tmp_path / "bmi.loco"
    loco_file.write_text(
        "FID_IID IEU2_IEU2 IEU1_IEU1\n"
        "22 0.4 0.3\n"
    )

    preds, valid_mask = load_loco_predictions(
        str(loco_file),
        ["bmi"],
        [["IEU1", "IEU1_IEU1"], ["IEU2", "IEU2_IEU2"], ["IEU3", "IEU3_IEU3"]],
    )

    np.testing.assert_array_equal(valid_mask, np.array([True, True, False]))
    np.testing.assert_allclose(preds["bmi"]["22"], np.array([0.3, 0.4]))


def test_loco_adjustment_changes_residualized_rif():
    taus = np.array([0.25, 0.5, 0.75])
    x = np.column_stack([np.ones(5), np.linspace(-1.0, 1.0, 5)])
    y_base = {"bmi": np.array([0.0, 1.0, 2.0, 3.0, 4.0])}
    y_loco = {"bmi": y_base["bmi"] - np.array([0.0, 0.2, 0.2, 0.0, 0.0])}

    rif_base, _, q_base = build_rif_residuals(y_base, ["bmi"], x, taus)
    rif_loco, _, q_loco = build_rif_residuals(y_loco, ["bmi"], x, taus)

    assert not np.allclose(rif_base[0], rif_loco[0])
    assert not np.allclose(q_base[0], q_loco[0])


def test_cpp_requirement_fails_closed_when_extension_missing():
    original_have_cpp = core.HAVE_CPP
    original_allow = core.ALLOW_NUMPY_FALLBACK
    try:
        core.HAVE_CPP = False
        core.set_numpy_fallback_allowed(False)
        try:
            core.require_cpp_extension("test")
        except RuntimeError as exc:
            assert "--allow-numpy-fallback" in str(exc)
        else:
            raise AssertionError("Expected RuntimeError when C++ extension is unavailable")
    finally:
        core.HAVE_CPP = original_have_cpp
        core.ALLOW_NUMPY_FALLBACK = original_allow
