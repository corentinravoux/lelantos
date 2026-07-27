"""Unit tests for :mod:`lelantos.task_manager` (no real solver execution)."""

import os

import pytest

from lelantos import task_manager


# --- Machine log parsing ---------------------------------------------------


def test_bash_machine_is_finished_and_error():
    machine = task_manager.Bash(ending_str="Total time", error_str=["srun:"])
    assert machine.is_finished(["setup", "Total time", "done"]) is True
    assert machine.is_finished(["setup", "running"]) is False
    assert machine.gives_error(["srun: error happened"]) is True
    assert machine.gives_error(["all good"]) is False


def test_bash_create_launcher_is_noop():
    machine = task_manager.Bash()
    assert machine.create_launcher("pwd", "dir", "cmd", "name") == ()


def test_load_cluster_optional_arguments_defaults():
    machine = task_manager.Nersc(N=2, n=4, c=8)
    assert machine.load_cluster_optional_arguments() == (2, 4, 8)
    assert task_manager.Nersc().load_cluster_optional_arguments() == (1, 1, 1)


# --- Dachshund input writing ----------------------------------------------


def test_dachshund_create_input_writes_cfg(workdir):
    solver = task_manager.Dachshund()
    assert solver.name == "dachshund"
    assert solver.ending_str == "Total time"
    params = [
        {
            "lx": 10.0, "ly": 12.0, "lz": 14.0, "npix": 100,
            "nx": 5, "ny": 6, "nz": 7, "sigmaf": 0.1, "lperp": 13, "lpar": 13,
            "namepixel": "pixel.bin", "namemap": "map.bin",
        }
    ]
    solver.create_input([str(workdir)], params, ["dachshund.cfg"])
    cfg = workdir / "dachshund.cfg"
    assert cfg.is_file()
    text = cfg.read_text()
    assert "lx = 10.0" in text
    assert "map_nx = 5" in text
    assert "corr_var_s = 0.1" in text
    assert "pixel.bin" in text and "map.bin" in text


# --- Manager directory helpers --------------------------------------------


def test_manager_dir_helpers(workdir):
    task_manager.TomographyManager.create_python_dir(str(workdir))
    assert (workdir / "Tmp").is_dir()
    task_manager.TomographyManager.create_dir(str(workdir), ["000000", "000001"])
    assert (workdir / "Tmp" / "000000").is_dir()
    assert (workdir / "Tmp" / "000001").is_dir()
    fname = str(workdir / "log.txt")
    task_manager.TomographyManager.create_file(fname)
    assert os.path.isfile(fname)


def test_manager_init_bash(workdir):
    manager = task_manager.TomographyManager(
        str(workdir), "dachshund", "bash", "pixel.bin", "launch.pickle", n=1
    )
    assert isinstance(manager.software, task_manager.Dachshund)
    assert isinstance(manager.machine, task_manager.Bash)
    # the bash machine inherits the solver completion marker
    assert manager.machine.ending_str == "Total time"


def test_manager_unknown_backends_return_keyerror(workdir):
    manager = task_manager.TomographyManager(
        str(workdir), "dachshund", "bash", "pixel.bin", "launch.pickle"
    )
    assert isinstance(manager.init_sofware("bogus"), KeyError)
    assert isinstance(manager.init_machine("bogus"), KeyError)
