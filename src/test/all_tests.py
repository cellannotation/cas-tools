import os
import sys
import unittest

loader = unittest.TestLoader()
suite = loader.discover(start_dir=os.getcwd(), pattern="*_test.py")

runner = unittest.TextTestRunner()
result = runner.run(suite).wasSuccessful()

sys.exit(not result)
