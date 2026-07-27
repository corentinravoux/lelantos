"""Unit tests for the config-parsing helpers in :mod:`lelantos.interface`."""

import pytest

from lelantos import interface


def test_parse_int_tuple():
    assert interface.parse_int_tuple("2, 2, 300") == (2, 2, 300)
    assert interface.parse_int_tuple("None") is None


def test_parse_float_tuple():
    assert interface.parse_float_tuple("7, 50") == (7.0, 50.0)
    assert interface.parse_float_tuple("None") is None


def test_parse_str_tuple():
    assert interface.parse_str_tuple("snr , delta") == ("snr", "delta")
    assert interface.parse_str_tuple("None") is None


def test_parse_dict_single_quotes():
    assert interface.parse_dict("{'n': 8, 'method': 'SPHERICAL'}") == {
        "n": 8,
        "method": "SPHERICAL",
    }
    assert interface.parse_dict("None") is None


def test_parse_float():
    assert interface.parse_float("0.3147") == pytest.approx(0.3147)
    assert interface.parse_float("None") is None


def test_parse_int():
    assert interface.parse_int("42") == 42
    assert interface.parse_int("None") is None


def test_parse_string():
    assert interface.parse_string("dachshund") == "dachshund"
    assert interface.parse_string("None") is None
