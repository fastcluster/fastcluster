import unittest

class fastcluster_test(unittest.TestCase):
    def test(self):
        from tests.test import test
        self.assertTrue(test(10))

    def test_nan(self):
        from tests.nantest import test
        self.assertTrue(test())

    def test_vector(self):
        from tests.vectortest import test
        self.assertTrue(test(10))
