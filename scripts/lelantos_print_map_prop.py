#!/usr/bin/env python3
from lelantos import interface
import sys


if __name__ == "__main__":
    interface.print_approximate_shape_size_from_interface_file(sys.argv[1])
