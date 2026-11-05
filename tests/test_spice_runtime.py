"""Exercise simulator selection without launching or installing a simulator."""
from pathlib import Path
import hashlib
import plistlib
import tempfile
import unittest

from ion_signal.spice_runtime import resolve_simulator, simulator_identity, source_path
from ion_signal.spice_storage import source_manifest_matches


class SimulatorSelection(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.executables = []
        for name in ("explicit", "environment", "local"):
            path = self.root / name
            path.write_text(name)
            path.chmod(0o755)
            self.executables.append(path)
        self.card = self.root / "config/local/tools.toml"
        self.card.parent.mkdir(parents=True)
        self.card.write_text('[ltspice]\nexecutable = "local"\n')

    def test_precedence_and_relative_local_path(self):
        environment = {"LTSPICE": str(self.executables[1])}
        self.assertEqual(resolve_simulator(self.executables[0], root=self.root,
                                          environment=environment), self.executables[0])
        self.assertEqual(resolve_simulator(root=self.root, environment=environment), self.executables[1])
        self.assertEqual(resolve_simulator(root=self.root, environment={}), self.executables[2])

    def test_pin_detects_replaced_binary(self):
        pin = hashlib.sha256(self.executables[2].read_bytes()).hexdigest()
        self.card.write_text(f'[ltspice]\nexecutable = "local"\nexpected_sha256 = "{pin}"\n')
        self.assertEqual(resolve_simulator(root=self.root, environment={}), self.executables[2])
        self.executables[2].write_text("different version")
        with self.assertRaisesRegex(ValueError, "pinned"):
            resolve_simulator(root=self.root, environment={})

    def test_missing_configuration_and_bad_explicit_path(self):
        self.card.unlink()
        self.assertIsNone(resolve_simulator(root=self.root, environment={}, required=False))
        with self.assertRaisesRegex(FileNotFoundError, "Configure"):
            resolve_simulator(root=self.root, environment={})
        with self.assertRaisesRegex(FileNotFoundError, "missing"):
            resolve_simulator(self.root / "absent", root=self.root, environment={}, required=False)

    def test_bundle_path_and_version_identity(self):
        bundle = self.root / "Pinned.app"
        executable = bundle / "Contents/MacOS/LTspice"
        executable.parent.mkdir(parents=True)
        executable.write_text("native binary placeholder")
        executable.chmod(0o755)
        (bundle / "Contents/Info.plist").write_bytes(plistlib.dumps({"CFBundleShortVersionString": "17.2.4"}))
        self.assertEqual(resolve_simulator(bundle, root=self.root, environment={}), executable)
        self.assertEqual(simulator_identity(executable)["version"], "17.2.4")

    def test_non_executable_file_is_rejected(self):
        self.executables[0].chmod(0o644)
        with self.assertRaises(FileNotFoundError):
            resolve_simulator(self.executables[0], root=self.root, environment={})

    def test_external_runtime_provenance(self):
        self.assertEqual(source_path(self.root / "local", root=self.root), "local")
        external = self.root.parent / "installed/LTspice"
        self.assertEqual(source_path(external, root=self.root), str(external))

    def test_relocated_source_manifest_keeps_hash_requirement(self):
        alias = self.root / "previous"
        alias.symlink_to("local")
        saved = {"previous": "original-content-hash"}
        self.assertTrue(source_manifest_matches(saved, {"local": "original-content-hash"}, self.root))
        self.assertFalse(source_manifest_matches(saved, {"local": "changed-content-hash"}, self.root))
        self.assertFalse(source_manifest_matches(saved, {}, self.root))


if __name__ == "__main__":
    unittest.main()
