"""Regenerate a saved matched report without a working-tree configuration."""
import importlib.util
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/elm_matched_production"


@unittest.skipUnless((SOURCE / "matched_analysis.h5").exists(), "requires local matched reduction")
class SavedReportTest(unittest.TestCase):
    def test_saved_branch_cards_and_missing_timing(self):
        import h5py
        sys.path.insert(0, str(ROOT / "scripts/production"))
        spec = importlib.util.spec_from_file_location(
            "matched_presentation", ROOT / "scripts/production/present_matched_analysis.py")
        report = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(report)
        with tempfile.TemporaryDirectory(dir=ROOT / "work") as directory:
            base = Path(directory)
            for name in ("matched_analysis.h5", "electron_analysis.h5", "manifest.toml"):
                shutil.copy2(SOURCE / name, base / name)
            with h5py.File(base / "matched_analysis.h5", "r+") as handle:
                del handle["metadata/run_timings_toml"]
                handle["metadata/run_timings_toml"] = ""
                old = handle["metadata/electron_branch"][()].decode()
                new = "tau_saved_snapshot/" + old.split("/")[1]
                del handle["metadata/electron_branch"]
                handle["metadata/electron_branch"] = new
            with h5py.File(base / "electron_analysis.h5", "r+") as handle:
                handle.move(old.split("/")[0], "tau_saved_snapshot")
                card = handle[new + "/realisation_card"][()].decode()
            report.ROOT = base  # No config/ directory exists here.
            report.OUT = base
            report.LATEX = base / "report"
            with patch.object(report, "build_report"):
                report.main()
            html = (base / "elm_matched_production.html").read_text()
            tex = (base / "report/elm_matched_production.tex").read_text()
            self.assertIn("Workflow timing was not recorded", html)
            self.assertIn(card, tex)
            self.assertTrue((base / "matched_pulses.png").exists())


if __name__ == "__main__":
    unittest.main()
