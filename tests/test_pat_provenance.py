"""Small fixtures for run identity, portable publication and refusal behavior."""
import argparse
import csv
import json
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflows/pat"))
import prepare_pat as pat


class ProvenanceTests(unittest.TestCase):
    def test_cross_node_device_requires_matching_content(self):
        record = pat.file_identity(self.bam, hash_content=True)
        record["device"] += 1
        pat.check_identity(self.bam, record)
        record["sha256"] = "0" * 64
        with self.assertRaises(pat.LineageInputError):
            pat.check_identity(self.bam, record)
        del record["sha256"]
        with self.assertRaises(pat.LineageInputError):
            pat.check_identity(self.bam, record)

    def test_cross_node_device_does_not_relax_inode_check(self):
        record = pat.file_identity(self.bam, hash_content=True)
        record["device"] += 1
        record["inode"] += 1
        with self.assertRaises(pat.LineageInputError):
            pat.check_identity(self.bam, record)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "E31"
        bamdir = self.source / "bismark/deduplicated"
        bamdir.mkdir(parents=True)
        self.bam = bamdir / "S1.deduplicated.sorted.bam"
        self.bam.write_bytes(b"source-bam")
        Path(f"{self.bam}.bai").write_bytes(b"index")
        self.bundle = self.root / "bundle"
        (self.bundle / "pat").mkdir(parents=True)
        self.manifest = {"S1": {"sample": "S1", "pat": "pat/S1.pat.gz",
                                "csi": "pat/S1.pat.gz.csi", "lbeta": "pat/S1.lbeta"}}
        for kind in ("pat", "csi", "lbeta"):
            (self.bundle / self.manifest["S1"][kind]).write_bytes(b"12345678")
        self.reference = {"cpg_bed_sha256": "expected", "cpg_count": 2}
        self.versions = {"wgbstools": "wgbstools version 0.3.0", "samtools": "samtools 1.22.1"}
        self.prov = {"schema_version": 2, "results_dir": str(self.source), "sample_order": ["S1"],
                     "genome": "hg38", "conversion_mode": "rebuild", "reference": self.reference.copy(),
                     "versions": self.versions.copy(), "bam2pat_defaults": pat.BAM2PAT_DEFAULTS.copy(),
                     "bam_inputs": [{"sample": "S1", "bam": str(self.bam), "bam_bytes": self.bam.stat().st_size,
                                     "bai": str(self.bam)+".bai", "bam_identity": pat.file_identity(self.bam, hash_content=True),
                                     "bai_identity": pat.file_identity(Path(f"{self.bam}.bai"), hash_content=True)}],
                     "pat_artifacts": [{**self.manifest["S1"], "mode": "rebuilt", "pat_bytes": 8, "lbeta_bytes": 8,
                                        "identities": {kind: pat.file_identity(self.bundle / self.manifest["S1"][kind], hash_content=True)
                                                       for kind in ("pat", "csi", "lbeta")}}]}
        self.write_bundle()
        mock_run = patch.object(pat, "run")
        self.run = mock_run.start(); self.run.return_value.stdout = "chr1\n"
        self.addCleanup(mock_run.stop)

    def write_bundle(self):
        (self.bundle / "pat_provenance.json").write_text(json.dumps(self.prov))
        (self.bundle / "PAT_COMPLETED").write_text("done\n")
        with (self.bundle / "pat_manifest.tsv").open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample", "pat", "csi", "lbeta"], delimiter="\t")
            writer.writeheader(); writer.writerows(self.manifest.values())

    def check(self, *, legacy=False):
        prov, manifest = pat.read_bundle(self.bundle, legacy=legacy)
        return pat.validate_reuse("S1", self.bam, self.bundle, prov, manifest, self.reference,
                                  self.versions, Path("tabix"), legacy=legacy)

    def test_valid_bundle(self):
        self.assertEqual(set(self.check()), {"pat", "csi", "lbeta"})

    def test_old_run_same_sample_is_rejected(self):
        for run in ("E20", "E21", "E27"):
            with self.subTest(run=run):
                self.prov["bam_inputs"][0]["bam"] = str(self.bam).replace("E31", run)
                self.write_bundle()
                with self.assertRaisesRegex(pat.LineageInputError, "Source BAM"):
                    self.check()

    def test_false_source_run_rejected(self):
        self.prov["results_dir"] = str(self.root / "E27"); self.write_bundle()
        with self.assertRaisesRegex(pat.LineageInputError, "Source methylseq run"):
            self.check()

    def test_wrong_dictionary_same_length_rejected(self):
        self.prov["reference"]["cpg_bed_sha256"] = "wrong"; self.write_bundle()
        with self.assertRaises(pat.LineageInputError): self.check()

    def test_changed_conversion_policy_rejected(self):
        self.prov["bam2pat_defaults"]["mapq"] = 20; self.write_bundle()
        with self.assertRaises(pat.LineageInputError): self.check()

    def test_changed_version_rejected(self):
        self.prov["versions"]["wgbstools"] = "different"; self.write_bundle()
        with self.assertRaises(pat.LineageInputError): self.check()

    def test_missing_and_unknown_samples_rejected(self):
        with self.assertRaisesRegex(pat.LineageInputError, "Missing BAM"):
            pat.discover_bams(self.source, [{"sample": "S1"}, {"sample": "missing"}])
        with self.assertRaisesRegex(pat.LineageInputError, "absent from metadata"):
            pat.discover_bams(self.source, [{"sample": "missing"}])

    def test_manifest_duplicate_rejected(self):
        with (self.bundle / "pat_manifest.tsv").open("a") as handle:
            handle.write("S1\tpat/S1.pat.gz\tpat/S1.pat.gz.csi\tpat/S1.lbeta\n")
        with self.assertRaisesRegex(pat.LineageInputError, "duplicate"):
            self.check()

    def test_missing_provenance_sample_rejected(self):
        self.prov["sample_order"] = ["S1", "missing"]; self.write_bundle()
        with self.assertRaises(pat.LineageInputError): self.check()

    def test_stale_manifest_rejected(self):
        (self.bundle / "pat/S1.pat.gz").unlink()
        with self.assertRaisesRegex(pat.LineageInputError, "Stale"):
            self.check()

    def test_unsafe_manifest_rejected(self):
        for value in ("../outside", "/absolute/S1.pat.gz"):
            self.manifest["S1"]["pat"] = value; self.write_bundle()
            with self.assertRaises(pat.LineageInputError): self.check()

    def test_changed_index_rejected(self):
        (self.bundle / "pat/S1.pat.gz.csi").write_bytes(b"87654321")
        with self.assertRaises(pat.LineageInputError): self.check()

    def test_bam_replacement_same_size_rejected(self):
        replacement = self.bam.with_suffix(".replacement")
        replacement.write_bytes(b"changedbam")
        replacement.replace(self.bam)
        with self.assertRaises(pat.LineageInputError): self.check()

    def test_publication_by_rename(self):
        target = self.root / "published"
        self.bundle.rename(target); self.bundle = target
        self.assertEqual(len(self.check()), 3)

    def test_legacy_requires_explicit_reconciliation(self):
        self.prov.pop("schema_version"); self.write_bundle()
        with self.assertRaisesRegex(pat.LineageInputError, "Legacy"):
            self.check()
        self.assertEqual(len(self.check(legacy=True)), 3)

    def test_incomplete_bundle_rejected(self):
        (self.bundle / "PAT_COMPLETED").unlink()
        with self.assertRaisesRegex(pat.LineageInputError, "Incomplete"):
            self.check()

    def test_complete_bundle_reuses_one_and_builds_only_missing(self):
        second = self.bam.with_name("S2.deduplicated.sorted.bam")
        second.write_bytes(b"secondbam")
        Path(f"{second}.bai").write_bytes(b"index2")
        metadata = self.root / "metadata.tsv"
        metadata.write_text("sample\tpatient\tcompartment\nS1\tP1\tBM\nS2\tP1\tPB\n")
        toolroot = self.root / "tools"
        executable = toolroot / "src/python/wgbs_tools.py"
        executable.parent.mkdir(parents=True)
        executable.write_text("tool"); executable.chmod(0o755)
        ref = toolroot / "references/hg38"; ref.mkdir(parents=True)
        (ref / "CpG.bed.gz").write_bytes(b"dictionary")
        (ref / "CpG.chrome.size").write_text("chr1\t29401795\n")
        self.prov["reference"] = {"cpg_bed_sha256": pat.HG38_CPG_SHA256, "cpg_count": 29401795}
        self.prov["bam_inputs"][0]["primary_mapped_reads"] = 100
        self.prov.pop("schema_version"); self.write_bundle()
        args = argparse.Namespace(wgbstools=executable, samtools=executable, tabix=executable,
            metadata=metadata, results_dir=self.source, allow_subset=False, reference_root=toolroot,
            genome="hg38", reference_manifest=None, conversion_mode="complete", reuse_sample=["S1"],
            reuse_dir=[], legacy_reuse_dir=[self.bundle], output_dir=self.root / "new",
            preflight_only=False, threads=1)
        conversion_samples = []
        def fake_run(command, **kwargs):
            if "bam2pat" in command:
                sample = Path(command[2]).stem; conversion_samples.append(sample)
                output = Path(command[command.index("--out_dir") + 1])
                for suffix in ("pat.gz", "pat.gz.csi", "lbeta"):
                    (output / f"{sample}.{suffix}").write_bytes(b"new-data")
            class Result:
                stdout = "chr1\n"
            result = Result()
            if command[-1] == "--version":
                result.stdout = self.versions["wgbstools" if fake_run.version_calls == 0 else "samtools"]
                fake_run.version_calls += 1
            return result
        fake_run.version_calls = 0
        real_sha = pat.sha256
        def fake_sha(path):
            return pat.HG38_CPG_SHA256 if path.name == "CpG.bed.gz" else real_sha(path)
        def tiny_pat(sample, directory, tabix, expected):
            return {kind: directory / f"{sample}.{suffix}" for kind, suffix in
                    (("pat", "pat.gz"), ("csi", "pat.gz.csi"), ("lbeta", "lbeta"))}
        def bam_check(sample, bam, samtools):
            return {"sample": sample, "bam": str(bam), "bam_bytes": bam.stat().st_size,
                    "bai": str(bam)+".bai", "primary_mapped_reads": 100}
        before = {path.name: path.read_bytes() for path in (self.bundle / "pat").iterdir()}
        with patch.object(pat, "run", side_effect=fake_run), patch.object(pat, "sha256", side_effect=fake_sha), \
             patch.object(pat, "validate_pat", side_effect=tiny_pat), patch.object(pat, "validate_bam", side_effect=bam_check):
            pat.prepare(args)
        self.assertEqual(conversion_samples, ["S2"])
        self.assertTrue((args.output_dir / "pat/S1.pat.gz").is_symlink())
        self.assertFalse((args.output_dir / "pat/S2.pat.gz").is_symlink())
        prov, manifest = pat.read_bundle(args.output_dir)
        self.assertEqual(list(manifest), ["S1", "S2"])
        self.assertEqual(len(prov["legacy_reconciliation"]), 1)
        self.assertIn("started_utc", prov); self.assertIn("completed_utc", prov)
        timing=json.loads((args.output_dir / "pat_execution.json").read_text())
        self.assertEqual(timing["status"], "completed")
        self.assertEqual(timing["started_utc"], prov["started_utc"])
        self.assertEqual(before, {path.name: path.read_bytes() for path in (self.bundle / "pat").iterdir()})
        fake_run.version_calls = 0
        with self.assertRaisesRegex(pat.LineageInputError, "Output must be absent"):
            with patch.object(pat, "run", side_effect=fake_run), patch.object(pat, "sha256", side_effect=fake_sha), \
                 patch.object(pat, "validate_pat", side_effect=tiny_pat):
                pat.prepare(args)

        args.output_dir = self.root / "all_new"
        args.reuse_sample = []; args.legacy_reuse_dir = []; args.reuse_dir = []
        conversion_samples.clear(); fake_run.version_calls = 0
        with patch.object(pat, "run", side_effect=fake_run), patch.object(pat, "sha256", side_effect=fake_sha), \
             patch.object(pat, "validate_pat", side_effect=tiny_pat), patch.object(pat, "validate_bam", side_effect=bam_check):
            pat.prepare(args)
        self.assertEqual(conversion_samples, ["S1", "S2"])
        self.assertFalse((args.output_dir / "pat/S1.pat.gz").is_symlink())

        args.output_dir = self.root / "failed"
        fake_run.version_calls = 0
        with patch.object(pat, "run", side_effect=fake_run), patch.object(pat, "sha256", side_effect=fake_sha), \
             patch.object(pat, "validate_bam", side_effect=pat.LineageInputError("fixture failure")):
            with self.assertRaisesRegex(pat.LineageInputError, "fixture failure"):
                pat.prepare(args)
        timing=json.loads((args.output_dir / "pat_execution.json").read_text())
        self.assertEqual(timing["status"], "failed_or_interrupted")
        self.assertIn("ended_utc", timing); self.assertNotIn("completed_utc", timing)
        self.assertFalse((args.output_dir / "PAT_COMPLETED").exists())

        args.output_dir = self.root / "failed_manifest"
        fake_run.version_calls = 0
        original_open = Path.open
        def fail_manifest(path, *open_args, **kwargs):
            if path == args.output_dir / "pat_manifest.tsv":
                raise OSError("fixture manifest failure")
            return original_open(path, *open_args, **kwargs)
        with patch.object(pat, "run", side_effect=fake_run), patch.object(pat, "sha256", side_effect=fake_sha), \
             patch.object(pat, "validate_pat", side_effect=tiny_pat), patch.object(pat, "validate_bam", side_effect=bam_check), \
             patch.object(Path, "open", fail_manifest):
            with self.assertRaisesRegex(OSError, "fixture manifest failure"):
                pat.prepare(args)
        prov=json.loads((args.output_dir / "pat_provenance.json").read_text())
        self.assertNotIn("completed_utc", prov)
        self.assertEqual(prov["status"], "failed_or_interrupted")

    def test_direct_output_names_and_internal_owner(self):
        with self.assertRaisesRegex(pat.LineageInputError, 'requires YYYY'):
            pat.validate_cli_output(self.root / 'undated')
        with self.assertRaises(ValueError):
            pat.validate_cli_output(self.root / '2026-02-30_12-00-invalid-v1')
        owner = self.root / '2026-11-02_13-00-fixture-v1'
        owner.mkdir()
        receipt = owner / '.run-ownership.json'
        receipt.write_text(json.dumps({'output': str(owner), 'status': 'running'}))
        pat.validate_cli_output(owner / 'prepared', owner)
        with self.assertRaisesRegex(pat.LineageInputError, 'does not belong'):
            pat.validate_cli_output(self.root / 'foreign', owner)
        receipt.write_text(json.dumps({'output': str(owner), 'status': 'completed'}))
        with self.assertRaisesRegex(pat.LineageInputError, 'immutable'):
            pat.validate_cli_output(owner / 'prepared', owner)

    def test_direct_output_parent_alias_is_rejected(self):
        alias = self.root / 'alias'
        alias.symlink_to(self.source, target_is_directory=True)
        with self.assertRaisesRegex(pat.LineageInputError, 'canonical'):
            pat.validate_cli_output(alias / '2026-11-02_13-00-example-v1')

    def test_symlink_never_replaces_existing_file(self):
        destination = self.root / "keep"
        destination.write_text("keep")
        with self.assertRaises(FileExistsError): pat.atomic_symlink(self.bam, destination)
        self.assertEqual(destination.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
