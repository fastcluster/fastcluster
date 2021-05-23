import unittest

class fastcluster_test(unittest.TestCase):
    def test(self):
        from tests.test import test
        test(10)

    def test_nan(self):
        from tests.nantest import test
        test()

    def test_vector(self):
        from tests.vectortest import test
        test(10)
