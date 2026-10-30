"""Tests for the terminal production-card form and its saved TOML."""
from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import tomllib
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts/production"))
import run_pilot  # noqa: E402


class ProductionFormTests(unittest.TestCase):
    def test_dump_round_trip(self):
        card = {"schema": "test", "count": 2, "enabled": True,
                "selected": ["bb0nu", "single_electrons"],
                "gas": {"pressure_bar": 15.0}}
        self.assertEqual(tomllib.loads(run_pilot.toml_dump(card)), card)

    def test_defaults_generate_five_cards(self):
        values = run_pilot.collect_values(run_pilot.templates(), interactive=False)
        self.assertEqual(values["j_variance_ratio"], 0.01)
        self.assertEqual(values["j_status"], "oliveira_fig6_user_reading")
        with tempfile.TemporaryDirectory() as directory:
            path = run_pilot.write_cards(values, Path(directory))
            self.assertEqual(len(list(path.parent.glob("*.toml"))), 5)
            run = tomllib.loads(path.read_text())
            self.assertEqual(run["selected_samples"], ["bb0nu", "single_electrons"])
            self.assertEqual(run["pilot_only"], True)
            realisation = tomllib.loads((path.parent / "realisation.toml").read_text())
            self.assertEqual(realisation["chamber_length_cm"], 150.0)
            electrons = tomllib.loads((path.parent / "electrons.toml").read_text())
            self.assertEqual(electrons["j_variance_ratio"], 0.01)

    def test_figure_six_anchors_and_intermediate_value(self):
        self.assertEqual(run_pilot.figure6_j_reading(1.5), 0.1)
        self.assertEqual(run_pilot.figure6_j_reading(2.0), 0.05)
        self.assertEqual(run_pilot.figure6_j_reading(3.0), 0.01)
        self.assertIsNone(run_pilot.figure6_j_reading(2.5))
        values = run_pilot.collect_values(run_pilot.templates(), interactive=False)
        values["reduced_field_kV_cm_bar"] = 2.5
        with self.assertRaisesRegex(ValueError, "matching E/P anchor"):
            run_pilot.validate_values(values)
        values["j_status"] = "user_specified"
        values["j_variance_ratio"] = 0.03
        run_pilot.validate_values(values)

    def test_intermediate_field_requires_explicit_j_in_form(self):
        def answer(prompt):
            return "2.5" if prompt.startswith("EL E/P") else ""
        with patch("builtins.input", side_effect=answer):
            with self.assertRaisesRegex(ValueError, "Enter J explicitly"):
                run_pilot.collect_values(run_pilot.templates(), interactive=True)
        def specified_answer(prompt):
            if prompt.startswith("EL E/P"):
                return "2.5"
            return "0.03" if prompt.startswith("EL variance ratio J") else ""
        with patch("builtins.input", side_effect=specified_answer):
            values = run_pilot.collect_values(run_pilot.templates(), interactive=True)
        self.assertEqual(values["j_status"], "user_specified")
        self.assertEqual(values["j_variance_ratio"], 0.03)

    def test_sample_and_domain_validation(self):
        values = run_pilot.collect_values(run_pilot.templates(), interactive=False)
        values["samples"] = ["single_electrons"]
        run_pilot.validate_values(values)
        values["reduced_field_kV_cm_bar"] = 3.1
        with self.assertRaisesRegex(ValueError, "E/P"):
            run_pilot.validate_values(values)
        values["reduced_field_kV_cm_bar"] = 3.0
        values["min_ion_drift_cm"] = 150.0
        with self.assertRaisesRegex(ValueError, "drift bounds"):
            run_pilot.validate_values(values)


if __name__ == "__main__":
    unittest.main()
