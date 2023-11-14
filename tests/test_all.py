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