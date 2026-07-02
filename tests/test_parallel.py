import unittest

from peptacular.sequence.parallel import parallel_apply_internal


def _simple_double(x: int) -> int:
    """Simple function to double a number."""
    return x * 2


def _add_value(x: int, value: int = 0) -> int:
    """Function that adds a value to x."""
    return x + value


def _square(x: int) -> int:
    """Simple function to square a number."""
    return x**2


class TestParallelApply(unittest.TestCase):
    def test_sequential_method(self):
        """Test sequential execution method."""
        items = [1, 2, 3, 4, 5]
        result = parallel_apply_internal(_simple_double, items, method="sequential")
        self.assertEqual(result, [2, 4, 6, 8, 10])

    def test_process_method(self):
        """Test multiprocessing execution method."""
        items = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        result = parallel_apply_internal(_simple_double, items, method="process", n_workers=2, chunksize=2)
        self.assertEqual(result, [2, 4, 6, 8, 10, 12, 14, 16, 18, 20])

    def test_thread_method(self):
        """Test threading execution method."""
        items = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        result = parallel_apply_internal(_simple_double, items, method="thread", n_workers=2, chunksize=2)
        self.assertEqual(result, [2, 4, 6, 8, 10, 12, 14, 16, 18, 20])

    def test_auto_method(self):
        """method=None auto-selects a backend and still returns correct, ordered results.

        Regression: auto-detection (thread on free-threaded Python, process otherwise)
        was documented but never wired up; method=None now routes through it.
        """
        from peptacular.sequence.parallel import _get_optimal_method

        self.assertIn(_get_optimal_method(), ("thread", "process"))
        items = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        result = parallel_apply_internal(_simple_double, items, method=None, n_workers=2)
        self.assertEqual(result, [2, 4, 6, 8, 10, 12, 14, 16, 18, 20])

    def test_with_kwargs(self):
        """Test parallel execution with keyword arguments."""
        items = [1, 2, 3, 4, 5]
        result = parallel_apply_internal(_add_value, items, method="sequential", value=10)
        self.assertEqual(result, [11, 12, 13, 14, 15])

    def test_order_preservation(self):
        """Test that results maintain input order."""
        items = list(range(20))
        result = parallel_apply_internal(_square, items, method="process", n_workers=4)
        expected = [x**2 for x in items]
        self.assertEqual(result, expected)

    def test_empty_input(self):
        """Test with empty input list."""
        items = []
        result = parallel_apply_internal(_simple_double, items, method="sequential")
        self.assertEqual(result, [])

    def test_single_item(self):
        """Test with single item."""
        items = [42]
        result = parallel_apply_internal(_simple_double, items, method="process")
        self.assertEqual(result, [84])

    def test_default_method(self):
        """Test with default method (should use process)."""
        items = [1, 2, 3, 4, 5]
        result = parallel_apply_internal(_simple_double, items)
        self.assertEqual(result, [2, 4, 6, 8, 10])

    def test_auto_chunksize(self):
        """Test automatic chunksize calculation."""
        items = list(range(100))
        result = parallel_apply_internal(_square, items, method="process", n_workers=4)
        expected = [x**2 for x in items]
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
