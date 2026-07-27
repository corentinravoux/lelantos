#!/usr/bin/env python3
"""Script: print the approximate map shape/size implied by a config file
(lelantos.interface.print_approximate_shape_size_from_interface_file)."""

from lelantos import interface
import sys


if __name__ == "__main__":
    interface.print_approximate_shape_size_from_interface_file(sys.argv[1])
