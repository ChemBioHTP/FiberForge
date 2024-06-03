"""
Unit and regression test for the fiberForge package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import fiberForge


def test_fiberForge_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "fiberForge" in sys.modules
