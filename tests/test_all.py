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
def run_tests(cat, group):
    # Run test suite of suplied group.
    test_loader = unittest.TestLoader()
    cwd = os.getcwd() + '/' + group
    pattern = group + '_'+ cat + '_test_*.py'
    print('cwd', cwd)
    print('pattern', pattern)
 #   input()
    test_suite = test_loader.discover(cwd, pattern= pattern)
    test_runner = unittest.TextTestRunner(verbosity = 2)
    result = test_runner.run(test_suite)

    # Report results, exit if there is an error.
    if not result.wasSuccessful():
        sys.exit(1)

    # Return number of tests ran and errors encountered.
    result_str = str(result)
    tests_run = int(result_str.split()[1][-1])
    errors = int(result_str.split()[2][-1])
    # errors = int(result_str)
    return tests_run, errors

if __name__ == '__main__':

    test_cats = ['unit', 'integration', 'verification']
    test_groups = ['plume', 'rpod', 'mission', 'mdao']

    # Run all tests in their groups.
    total_tests = 0
    total_errors = 0
    for group in test_groups:
        for cat in test_cats:
            tests_run, errors = run_tests(cat, group)
            total_tests += tests_run
            total_errors += errors

    # Report cummulative results.
    print(total_tests, "tests ran", total_errors, "total errors")
