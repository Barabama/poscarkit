import unittest

from src.utils.elements import canonicalize_element


class TestCanonicalizeElement(unittest.TestCase):
    """Test element symbol canonicalization."""

    def test_all_caps_two_letter(self):
        """All-caps export symbols are title-cased."""
        self.assertEqual(canonicalize_element("CO"), "Co")
        self.assertEqual(canonicalize_element("CU"), "Cu")
        self.assertEqual(canonicalize_element("NI"), "Ni")
        self.assertEqual(canonicalize_element("CR"), "Cr")
        self.assertEqual(canonicalize_element("FE"), "Fe")

    def test_single_letter_unchanged(self):
        """Single-letter symbols are unaffected."""
        self.assertEqual(canonicalize_element("V"), "V")
        self.assertEqual(canonicalize_element("C"), "C")
        self.assertEqual(canonicalize_element("N"), "N")

    def test_proper_case_unchanged(self):
        """Canonical symbols pass through unchanged."""
        self.assertEqual(canonicalize_element("Co"), "Co")
        self.assertEqual(canonicalize_element("Cu"), "Cu")
        self.assertEqual(canonicalize_element("Ni"), "Ni")
        self.assertEqual(canonicalize_element("Au"), "Au")

    def test_lowercase_input(self):
        """Lowercase input is corrected."""
        self.assertEqual(canonicalize_element("co"), "Co")
        self.assertEqual(canonicalize_element("ni"), "Ni")

    def test_mixed_case(self):
        """Mixed-case input is corrected to canonical."""
        self.assertEqual(canonicalize_element("cO"), "Co")

    def test_empty_string(self):
        """Empty string passes through unchanged."""
        self.assertEqual(canonicalize_element(""), "")


if __name__ == "__main__":
    unittest.main()
