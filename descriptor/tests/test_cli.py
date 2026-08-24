from __future__ import annotations

from elf_anisotropy.cli import run


def test_missing_elfcar_returns_clear_error(capsys, tmp_path) -> None:
    status = run(["analyze", str(tmp_path / "ELFCAR")])
    captured = capsys.readouterr()
    assert status == 2
    assert "ELFCAR file not found" in captured.err
