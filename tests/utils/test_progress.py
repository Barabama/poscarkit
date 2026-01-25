import time
import unittest
import io
import sys
from unittest.mock import patch, MagicMock

from src.utils.progress import SimpleProgress, progress, progress_context


class TestSimpleProgress(unittest.TestCase):
    """Test SimpleProgress class."""

    def test_fixed_length_iteration(self):
        """Test iteration over fixed-length iterable."""
        data = list(range(10))
        result = []

        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            for item in progress(data, desc="Test"):
                result.append(item)
                time.sleep(0.01)

        self.assertEqual(result, data)
        output = fake_out.getvalue()
        self.assertIn("Test:", output)

    def test_unknown_length_iteration(self):
        """Test iteration over unknown-length iterable."""
        def generator():
            for i in range(5):
                yield i
                time.sleep(0.01)

        result = []

        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            for item in progress(generator(), desc="Generator"):
                result.append(item)

        self.assertEqual(result, list(range(5)))
        output = fake_out.getvalue()
        self.assertIn("Generator:", output)

    def test_manual_update(self):
        """Test manual update with context manager."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            with progress(total=10, desc="Manual") as p:
                for i in range(10):
                    p.update(1)
                    time.sleep(0.01)

        output = fake_out.getvalue()
        self.assertIn("Manual:", output)

    def test_manual_update_unknown_total(self):
        """Test manual update with unknown total."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            with progress(total=None, desc="Unknown") as p:
                for i in range(5):
                    p.update(1)
                    time.sleep(0.01)

        output = fake_out.getvalue()
        self.assertIn("Unknown:", output)

    def test_wrap_method_fixed_length(self):
        """Test wrap method with fixed-length iterable."""
        data = list(range(10))
        result = []

        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            for item in SimpleProgress(desc="Wrap").wrap(data):
                result.append(item)
                time.sleep(0.01)

        self.assertEqual(result, data)
        output = fake_out.getvalue()
        self.assertIn("Wrap:", output)

    def test_wrap_method_unknown_length(self):
        """Test wrap method with unknown-length iterable."""
        def generator():
            for i in range(5):
                yield i
                time.sleep(0.01)

        result = []

        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            for item in SimpleProgress(desc="WrapGen").wrap(generator()):
                result.append(item)

        self.assertEqual(result, list(range(5)))
        output = fake_out.getvalue()
        self.assertIn("WrapGen:", output)

    def test_wrap_method_with_iterable_param(self):
        """Test wrap method with iterable parameter."""
        data = list(range(8))
        result = []

        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            p = SimpleProgress(desc="WrapParam")
            for item in p.wrap(data):
                result.append(item)
                time.sleep(0.01)

        self.assertEqual(result, data)
        output = fake_out.getvalue()
        self.assertIn("WrapParam:", output)

    def test_wrap_method_no_iterable(self):
        """Test wrap method with no iterable provided."""
        with patch('sys.stdout', new=io.StringIO()):
            p = SimpleProgress(desc="NoIter")
            with self.assertRaises(ValueError) as cm:
                list(p.wrap(None))
            self.assertIn("No iterable provided", str(cm.exception))

    def test_close(self):
        """Test close method."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            p = SimpleProgress(total=10, desc="Close")
            p.update(5)
            p.close()

        output = fake_out.getvalue()
        self.assertIn("Close:", output)

    def test_context_manager(self):
        """Test context manager usage."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            with SimpleProgress(total=10, desc="Context") as p:
                p.update(5)

        output = fake_out.getvalue()
        self.assertIn("Context:", output)

    def test_progress_function(self):
        """Test progress() function."""
        data = list(range(5))
        result = []

        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            for item in progress(data, desc="Func"):
                result.append(item)

        self.assertEqual(result, data)
        output = fake_out.getvalue()
        self.assertIn("Func:", output)

    def test_progress_context_function(self):
        """Test progress_context() function."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            with progress_context(total=10, desc="ContextFunc") as p:
                for i in range(10):
                    p.update(1)

        output = fake_out.getvalue()
        self.assertIn("ContextFunc:", output)

    def test_empty_iterable(self):
        """Test empty iterable."""
        data = []
        result = []

        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            for item in progress(data, desc="Empty"):
                result.append(item)

        self.assertEqual(result, [])
        output = fake_out.getvalue()
        self.assertIn("Empty:", output)

    def test_custom_unit(self):
        """Test custom unit."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            p = SimpleProgress(total=100, unit="bytes")
            p.update(50)

        output = fake_out.getvalue()
        self.assertIn("bytes", output)

    def test_no_description(self):
        """Test without description."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            for item in progress(range(5)):
                pass

        output = fake_out.getvalue()
        self.assertIn("it/s", output)

    def test_unit_scale(self):
        """Test unit scaling."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            data = range(1000000)
            for item in progress(data, desc="Scale", unit_scale=True):
                if item % 100000 == 0:
                    time.sleep(0.01)

        output = fake_out.getvalue()
        self.assertIn("Scale:", output)

    def test_ncols_parameter(self):
        """Test ncols parameter."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            p = SimpleProgress(total=100, desc="Width", ncols=80)
            p.update(50)

        output = fake_out.getvalue()
        self.assertIn("Width:", output)

    def test_mininterval_parameter(self):
        """Test mininterval parameter."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            p = SimpleProgress(total=100, mininterval=0.5)
            p.update(50)
            time.sleep(0.6)
            p.update(50)

        output = fake_out.getvalue()
        self.assertIn("100%", output)

    def test_dynamic_ncols_parameter(self):
        """Test dynamic_ncols parameter."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            p = SimpleProgress(total=100, desc="Dynamic", dynamic_ncols=True)
            p.update(50)

        output = fake_out.getvalue()
        self.assertIn("Dynamic:", output)

    def test_leave_parameter(self):
        """Test leave parameter."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            p = SimpleProgress(total=100, desc="Leave", leave=False)
            p.update(100)

        output = fake_out.getvalue()
        self.assertIn("Leave:", output)

    def test_update_multiple(self):
        """Test updating by multiple steps."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            p = SimpleProgress(total=100, desc="Multi")
            p.update(10)
            p.update(20)
            p.update(30)

        output = fake_out.getvalue()
        self.assertIn("Multi:", output)

    def test_context_manager_exception(self):
        """Test context manager with exception."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            try:
                with SimpleProgress(total=10, desc="Exception") as p:
                    p.update(5)
                    raise ValueError("Test exception")
            except ValueError:
                pass

        output = fake_out.getvalue()
        self.assertIn("Exception:", output)

    def test_tqdm_import_error(self):
        """Test ImportError when tqdm is not available."""
        with patch('src.utils.progress.TqdmProgress', None):
            with self.assertRaises(ImportError) as cm:
                SimpleProgress(total=10)
            self.assertIn("tqdm is required", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
