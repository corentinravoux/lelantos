#!/usr/bin/env python3
"""Script: run the full lelantos pipeline from a config file. Usage:
``lelantos_interface.py interface.ini`` (see interface_explanatory.ini)."""

from lelantos import interface
import sys

if __name__ == "__main__":
    interface.main(sys.argv[1])
