"""
test_domains_scan.py — offline unit tests for scripts/domains_scan.py

Mocks pfam_scan.scan_sequences so these tests don't require the actual
~2.2GB Pfam-A.hmm file or pyhmmer to be installed/exercised — they only
verify domains_scan.py's own cache-hit/miss and output-format logic.
"""

import domains_scan


def _write_fasta(path, name, seq):
    path.write_text(f">{name}\n{seq}\n")


class TestDomainsScanCacheMiss:

    def test_on_demand_scan_writes_output_and_cache(self, tmp_path, monkeypatch):
        ref_out_dir = tmp_path / "ref"
        ref_out_dir.mkdir()
        monkeypatch.setattr(domains_scan, "_load_reference_data_config",
                            lambda: {"out_dir": str(ref_out_dir)})
        monkeypatch.setattr(domains_scan, "scan_sequences",
                            lambda records: {"SNB-1": [("PF00957", "Synaptobrevin", 23, 103, 109.8)]})

        fasta_in = tmp_path / "snb1.fa"
        _write_fasta(fasta_in, "SNB-1", "MDAQ")
        outputfile = tmp_path / "out.txt"

        domains_scan.main(str(fasta_in), "", str(tmp_path), "", str(outputfile))

        assert outputfile.read_text() == "Pfam\t23\t103\tSynaptobrevin\n"
        cache_file = ref_out_dir / "pfam_scan_cache" / "SNB-1_domains.txt"
        assert cache_file.exists()
        assert cache_file.read_text() == outputfile.read_text()

    def test_zero_hits_writes_empty_file(self, tmp_path, monkeypatch):
        ref_out_dir = tmp_path / "ref"
        ref_out_dir.mkdir()
        monkeypatch.setattr(domains_scan, "_load_reference_data_config",
                            lambda: {"out_dir": str(ref_out_dir)})
        monkeypatch.setattr(domains_scan, "scan_sequences", lambda records: {"FOO": []})

        fasta_in = tmp_path / "foo.fa"
        _write_fasta(fasta_in, "FOO", "MDAQ")
        outputfile = tmp_path / "out.txt"

        domains_scan.main(str(fasta_in), "", str(tmp_path), "", str(outputfile))

        assert outputfile.read_text() == ""


class TestDomainsScanCacheHit:

    def test_precomputed_cache_file_is_reused_without_scanning(self, tmp_path, monkeypatch):
        ref_out_dir = tmp_path / "ref"
        cache_dir = ref_out_dir / "pfam_scan_cache"
        cache_dir.mkdir(parents=True)
        (cache_dir / "PEZO-1_domains.txt").write_text(
            "Pfam\t27\t769\tPiezo TM1-24\nPfam\t2113\t2436\tPiezo non-specific cation channel, cap domain\n"
        )
        monkeypatch.setattr(domains_scan, "_load_reference_data_config",
                            lambda: {"out_dir": str(ref_out_dir)})

        def _scan_sequences_should_not_be_called(records):
            raise AssertionError("scan_sequences should not run on a cache hit")
        monkeypatch.setattr(domains_scan, "scan_sequences", _scan_sequences_should_not_be_called)

        fasta_in = tmp_path / "pezo1.fa"
        _write_fasta(fasta_in, "PEZO-1", "M" * 2442)
        outputfile = tmp_path / "out.txt"

        domains_scan.main(str(fasta_in), "", str(tmp_path), "", str(outputfile))

        rows = outputfile.read_text().splitlines()
        assert len(rows) == 2
        assert rows[0] == "Pfam\t27\t769\tPiezo TM1-24"
