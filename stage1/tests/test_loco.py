from pathlib import Path

import numpy as np

from fungwas_stage1.cli import (
    build_rif_residuals,
    compute_chr_rif_residuals,
    load_loco_predictions,
    normalize_chr_label,
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


def test_normalize_chr_label_handles_numeric_and_special_labels():
    assert normalize_chr_label("04") == "4"
    assert normalize_chr_label("4") == "4"
    assert normalize_chr_label("chr04") == "4"
    assert normalize_chr_label("CHR22") == "22"
    assert normalize_chr_label("x") == "X"
    assert normalize_chr_label("chrxy") == "XY"
    assert normalize_chr_label("mt") == "MT"


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


def test_compute_chr_rif_residuals_accepts_zero_padded_chr_labels():
    taus = np.array([0.25, 0.5, 0.75])
    x = np.column_stack([np.ones(4), np.linspace(-1.0, 1.0, 4)])
    y_values = {"bmi": np.array([1.0, 2.0, 3.0, 4.0])}
    loco_predictions = {"bmi": {"4": np.array([0.1, 0.2, 0.3, 0.4])}}
    valid_mask = np.array([True, True, True, True])

    rif_resid_list, q_chr, q_tau_list = compute_chr_rif_residuals(
        ["bmi"], y_values, x, taus, loco_predictions, "04", valid_mask
    )

    assert len(rif_resid_list) == 1
    assert rif_resid_list[0].shape == (4, 3)
    assert q_chr.shape == x.shape
    assert len(q_tau_list) == 1
    assert q_tau_list[0].shape == (3,)


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


def test_numpy_fallback_warning_is_explicit(caplog):
    original_have_cpp = core.HAVE_CPP
    original_allow = core.ALLOW_NUMPY_FALLBACK
    try:
        core.HAVE_CPP = False
        with caplog.at_level("WARNING"):
            core.set_numpy_fallback_allowed(True)
        warning_text = "\n".join(record.message for record in caplog.records)
        assert "NUMPY FALLBACK ENABLED FOR STAGE 1" in warning_text
        assert "SUPER SUPER SLOW" in warning_text
        assert "DO NOT USE THIS FOR PRODUCTION GWAS" in warning_text
        assert "--allow-numpy-fallback" in warning_text
    finally:
        core.HAVE_CPP = original_have_cpp
        core.ALLOW_NUMPY_FALLBACK = original_allow


def test_multi_phenotype_block_scores_match_single_path_under_numpy_fallback():
    original_have_cpp = core.HAVE_CPP
    original_allow = core.ALLOW_NUMPY_FALLBACK
    try:
        core.HAVE_CPP = False
        core.set_numpy_fallback_allowed(True)

        rng = np.random.default_rng(123)
        n_samples = 120
        n_snps = 6
        n_taus = 3
        n_blocks = 8

        g = rng.binomial(2, 0.3, size=(n_samples, n_snps)).astype(np.float64)
        x = np.column_stack([np.ones(n_samples), rng.normal(size=(n_samples, 2))])
        q, _ = np.linalg.qr(x)
        block_ids = rng.integers(0, n_blocks, size=n_samples, dtype=np.int32)

        rif_a = rng.normal(size=(n_samples, n_taus))
        rif_b = rng.normal(size=(n_samples, n_taus))

        multi_stats = core.compute_block_scores_multi(
            g, [rif_a, rif_b], q, block_ids, n_blocks
        )
        single_stats_a = core.compute_block_scores(
            g, rif_a, q, block_ids, n_blocks
        )
        single_stats_b = core.compute_block_scores(
            g, rif_b, q, block_ids, n_blocks
        )

        assert len(multi_stats) == 2
        np.testing.assert_allclose(multi_stats[0], single_stats_a)
        np.testing.assert_allclose(multi_stats[1], single_stats_b)
    finally:
        core.HAVE_CPP = original_have_cpp
        core.ALLOW_NUMPY_FALLBACK = original_allow
