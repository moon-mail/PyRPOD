# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 12-06-23

# ========================
# PyRPOD: test_all.py
# ========================
# The script will read and run all test case scripts in numerical order.
# At minimum, this script should be run at the begining and end of each
# day to ensure new development dont break any old code.

# Good practice calls for re-running this any time code developments need to be tested or debugged.
# Can be multiple times a day or even multiple times per hour depending on the specific task.
# (Re-word to be more clear^?)

import test_header
import unittest, os, sys
def run_tests():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover(os.getcwd(), pattern='test_case_*.py')
    test_runner = unittest.TextTestRunner(verbosity = 2)
    result = test_runner.run(test_suite)
    return result

if __name__ == '__main__':
    result = run_tests()
    print(result)
    if not result.wasSuccessful():
        sys.exit(1)