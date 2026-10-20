"""Non-scientific mock-process tests for the Step 79 runner hardening."""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import run_vina_3L4W_alternative_ligand_pilot as runner  # noqa: E402


MOCK_PROGRAM = """\
from pathlib import Path
import sys
import time

behavior, output = sys.argv[1:3]
if behavior == "success":
    Path(output).write_text("MODEL 1\\nREMARK MOCK NON-SCIENTIFIC\\nENDMDL\\n", encoding="utf-8")
    raise SystemExit(0)
if behavior == "nonzero":
    print("mock failure", file=sys.stderr)
    raise SystemExit(7)
if behavior == "timeout":
    time.sleep(30)
raise SystemExit(9)
"""


class RunnerHardeningTests(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory(dir=ROOT)
        self.temp = Path(self._temporary.name)
        self.mock = self.temp / "mock_vina.py"
        self.mock.write_text(MOCK_PROGRAM, encoding="utf-8")
        self.receptor = self.temp / "receptor.pdbqt"
        self.ligand = self.temp / "qualified.pdbqt"
        self.start = self.temp / "start.pdbqt"
        self.receptor.write_text("RECEPTOR MOCK\n", encoding="utf-8")
        self.ligand.write_text("LIGAND MOCK\n", encoding="utf-8")
        self.start.write_text("START MOCK\n", encoding="utf-8")

    def tearDown(self) -> None:
        self._temporary.cleanup()

    def job(self, behavior: str = "success", seed: int = 1101) -> dict:
        folder = self.temp / "outputs" / "LIG02_neutral_D542minus"
        output = folder / f"seed_{seed}_docked.pdbqt"
        row = {
            "structure_id": "LIG02", "receptor_state": "Asp542-", "start_seed": str(seed),
            "receptor_pdbqt": str(self.receptor), "qualified_pdbqt": str(self.ligand),
            "start_pdbqt": str(self.start), "receptor_sha256": runner.base.digest(self.receptor),
            "qualified_pdbqt_sha256": runner.base.digest(self.ligand),
            "start_sha256": runner.base.digest(self.start), "num_modes": "20",
        }
        return {"row": row, "seed": seed, "folder": folder, "output": output,
                "cmd": [sys.executable, str(self.mock), behavior, str(output)],
                "control": self.temp / "control", "manifest_sha256": "a" * 64}

    def test_success_and_verified_reuse(self) -> None:
        job = self.job()
        result = runner.run_one(job, timeout_seconds=5)
        self.assertEqual(result["status"], "completed")
        classified = runner.classify_job(job)
        self.assertEqual(classified["classification"], "REUSE")
        self.assertIn("MOCK NON-SCIENTIFIC", job["output"].read_text(encoding="utf-8"))

    def test_nonzero_exit_is_finalized_and_not_reused(self) -> None:
        job = self.job("nonzero")
        result = runner.run_one(job, timeout_seconds=5)
        self.assertEqual(result["status"], "vina_failed")
        self.assertEqual(result["exit_code"], 7)
        self.assertEqual(runner.classify_job(job)["classification"], "MISSING")

    def test_timeout_kills_and_finalizes_mock_process(self) -> None:
        job = self.job("timeout")
        result = runner.run_one(job, timeout_seconds=0.2)
        self.assertEqual(result["status"], "timeout")
        record = json.loads(result["logs"][-1].read_text(encoding="utf-8"))
        self.assertEqual(record["exit_code"], "timeout_killed")
        self.assertLess(record["elapsed_seconds"], 10)

    def test_launch_without_final_record_is_blocked(self) -> None:
        job = self.job()
        launch = runner.launch_folder(job) / "seed_1101_attempt_1_launch.json"
        runner.atomic_json_new(launch, {"pid": 424242, "command": job["cmd"]})
        classified = runner.classify_job(job)
        self.assertEqual(classified["classification"], "BLOCKED")
        self.assertIn("no final record", classified["reason"])
        self.assertTrue(launch.is_file())

    def test_altered_output_hash_is_blocked(self) -> None:
        job = self.job()
        self.assertEqual(runner.run_one(job, timeout_seconds=5)["status"], "completed")
        with job["output"].open("a", encoding="utf-8") as stream:
            stream.write("REMARK ALTERED AFTER RECORD\n")
        classified = runner.classify_job(job)
        self.assertEqual(classified["classification"], "BLOCKED")
        self.assertIn("SHA256", classified["reason"])

    def test_invalid_attempt_json_is_blocked(self) -> None:
        job = self.job()
        job["folder"].mkdir(parents=True)
        attempt = job["folder"] / "seed_1101_attempt_1.json"
        attempt.write_text("{not valid JSON", encoding="utf-8")
        self.assertEqual(runner.classify_job(job)["classification"], "BLOCKED")

    def test_output_without_attempt_json_is_blocked_and_not_overwritten(self) -> None:
        job = self.job()
        job["folder"].mkdir(parents=True)
        original = b"MODEL 1\nREMARK EXISTING MOCK\nENDMDL\n"
        job["output"].write_bytes(original)

        def forbidden_popen(*args, **kwargs):
            raise AssertionError("mock process must not be started")

        result = runner.run_one(job, timeout_seconds=5, popen_factory=forbidden_popen)
        self.assertEqual(result["status"], "blocked")
        self.assertEqual(job["output"].read_bytes(), original)
        self.assertFalse(list(job["folder"].glob("*attempt*.json")))

    def test_atomic_write_never_overwrites(self) -> None:
        path = self.temp / "atomic.json"
        runner.atomic_write_new(path, b"first")
        with self.assertRaises(FileExistsError):
            runner.atomic_write_new(path, b"second")
        self.assertEqual(path.read_bytes(), b"first")
        self.assertFalse(list(self.temp.glob(".atomic.json.*.tmp")))

    def test_dry_run_counting(self) -> None:
        reuse = self.job("success", 1101)
        missing = self.job("success", 2202)
        blocked = self.job("success", 3303)
        self.assertEqual(runner.run_one(reuse, timeout_seconds=5)["status"], "completed")
        blocked["folder"].mkdir(parents=True, exist_ok=True)
        blocked["output"].write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
        items = [runner.classify_job(job) for job in (reuse, missing, blocked)]
        self.assertEqual(runner.classification_counts(items),
                         {"REUSE": 1, "MISSING": 1, "BLOCKED": 1})


if __name__ == "__main__":
    unittest.main()
