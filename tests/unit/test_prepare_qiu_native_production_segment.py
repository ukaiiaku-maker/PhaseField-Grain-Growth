from pathlib import Path

import numpy as np
import pytest

from scripts.prepare_qiu_native_production_segment import accepted_step, wrapper_text


def test_segment_wrapper_encodes_first_segment_without_replay():
    wrapper=wrapper_text(0,5000,"","")
    assert "QIU_SEGMENT_TARGET_STEP=5000" in wrapper
    assert "unset QIU_PRODUCTION_RESTART" in wrapper
    assert "native-restart-step005000.npz" in wrapper
    assert "baseline_eligible'] is True" in wrapper
    assert "NUMBA_NUM_THREADS=1 QIU_NUMBA_THREADS=1" in wrapper


def test_segment_wrapper_binds_exact_parent_hash():
    wrapper=wrapper_text(5000,10000,"parent/native-restart-step005000.npz","a"*64)
    assert "a"*64 in wrapper
    assert 'QIU_PRODUCTION_RESTART="$PWD/parent/native-restart-step005000.npz"' in wrapper
    assert "parent_checkpoint_sha256':'" + "a"*64 in wrapper


def test_parent_accepted_step_is_read_without_pickle(tmp_path):
    checkpoint=tmp_path/"restart.npz";np.savez(checkpoint,accepted_step=np.asarray(7500))
    assert accepted_step(checkpoint)==7500


def test_segment_bounds_are_preregistered_in_wrapper():
    wrapper=wrapper_text(199750,200000,"parent/restart.npz","b"*64)
    assert "target=200000" in wrapper
    assert "range(((start//250)+1)*250,target+1,250)" in wrapper
