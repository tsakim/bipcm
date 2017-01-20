.. _testing:

Testing
--------------------------------------------------------------------------------

The methods in the ``bipcm`` module  have been implemented using `doctests
<https://docs.python.org/2/library/doctest.html>`_. To run the tests,
execute::

    >>> python -m doctest bipcm_tests.txt

from the folder `src` in the command line. If you want to run the tests in
verbose mode, use::

    >>> python -m doctest -v bipcm_tests.txt

Note that `bipcm.py` and `bipcm_tests.txt` have to be in the same directory to
run the test.

