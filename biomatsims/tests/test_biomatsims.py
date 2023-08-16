"""
Unit and regression test for the biomatsims package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import biomatsims


def test_biomatsims_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "biomatsims" in sys.modules
