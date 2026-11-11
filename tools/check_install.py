"""Build from an sdist after removing staging; test a clean installed wheel.

Run with a Python containing setuptools and wheel. Use --wheelhouse to test
offline with downloaded dependency wheels; otherwise pip resolves dependencies.
All generated artifacts and environments are temporary.
"""

from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import venv


def run(args, cwd, env):
    subprocess.run([str(arg) for arg in args], cwd=cwd, env=env, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wheelhouse', type=Path, help='offline directory containing dependency wheels')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"
    with tempfile.TemporaryDirectory(prefix="topostich-install-") as directory:
        work = Path(directory)
        staging = work / "source"
        staging.mkdir()
        for name in ("pyproject.toml", "README.md", "MANIFEST.in", "AGENTS.md", "CHANGELOG.md"):
            shutil.copy2(root / name, staging / name)
        for name in ("src", "tests", "tools", "docs", "schemas"):
            shutil.copytree(root / name, staging / name,
                            ignore=shutil.ignore_patterns("__pycache__", "*.egg-info"))
        archives = work / "archives"
        archives.mkdir()
        run([sys.executable, "-c", "from setuptools.build_meta import build_sdist; build_sdist('../archives')"], staging, env)
        sdist, = archives.glob("*.tar.gz")
        shutil.rmtree(staging)
        relocated = work / "relocated"
        relocated.mkdir()
        with tarfile.open(sdist) as archive:
            # Only extract our freshly built archive, with explicit safe paths
            # for Python 3.10 as well as newer extraction-filter implementations.
            for member in archive.getmembers():
                target = (relocated / member.name).resolve()
                if not target.is_relative_to(relocated.resolve()) or not (member.isfile() or member.isdir()):
                    raise RuntimeError("Unsafe sdist archive member")
            archive.extractall(relocated)
        extracted, = relocated.iterdir()
        run([sys.executable, "-c", "from setuptools.build_meta import build_wheel; build_wheel('../../archives')"], extracted, env)
        wheel, = archives.glob("*.whl")
        installed = work / "installed"
        venv.EnvBuilder(with_pip=True, symlinks=os.name != "nt").create(installed)
        bindir = installed / ("Scripts" if os.name == "nt" else "bin")
        python = bindir / ("python.exe" if os.name == "nt" else "python")
        command = bindir / ("topostich.exe" if os.name == "nt" else "topostich")
        install_args = [python, "-m", "pip", "install"]
        if args.wheelhouse:
            install_args.extend(['--no-index', '--find-links', args.wheelhouse.resolve()])
        run([*install_args, wheel], work, env)
        run([python, '-m', 'pip', 'check'], work, env)
        run([command, "--version"], work, env)
        run([command, "--help"], work, env)
        smoke_input = work / "smoke.pdb"
        smoke_input.write_text(
            "CRYST1   20.000   30.000   40.000  90.00  90.00  90.00 P 1           1\n"
            "ATOM      1  O5'  DG A   1       1.000   2.000   3.000  1.00 20.00           O\n"
            "ATOM      2  O3'  DG A   1       4.000   5.000   6.000  1.00 20.00           O\n"
            "END\n"
        )
        smoke_output = work / "smoke.cif"
        run([command, smoke_input, "-atom", "O3'", "-atom", "O5'", "-out", smoke_output, "--json"], work, env)
        if "_atom_site_fract_x" not in smoke_output.read_text():
            raise RuntimeError("Installed CLI did not produce core CIF")
        smoke_mmcif = work / 'smoke.mmcif'
        run([command, smoke_input, '-out', smoke_mmcif, '-mmCif'], work, env)
        run([command, smoke_mmcif, '-atom', "O5'", '-atom', "O3'", '-out', work/'native.cif'], work, env)
        run([command, smoke_mmcif, '-atom', '-all', '-out', work/'all.cif'], work, env)
        cgd_command = bindir / ('topostich-cgd.exe' if os.name == 'nt' else 'topostich-cgd')
        run([cgd_command, '--version'], work, env)
        run([cgd_command, smoke_mmcif, '--out', work/'graph.cgd',
             '--cmm', work/'graph.cmm', '--crystal-cmm', work/'copies.cmm',
             '--receipt', work/'graph.json'], work, env)
        run([python, '-c', "import json,xml.etree.ElementTree as E; "
             "r=json.load(open('graph.json')); assert r['node_count']==1 and r['edge_count']==0; "
             "assert len(E.parse('graph.cmm').getroot().findall('.//marker'))==1"], work, env)
        run([python, '-c', "from topostich.nucleotides import classify_atom; assert classify_atom('DZ', 'C1').role == 'base'"], work, env)
        run([python, "-m", "unittest", "discover", "-s", extracted / "tests", "-v"], work, env)
        print("PASS: sdist, relocated wheel build, clean install, installed CLI and complete installed suite")


if __name__ == "__main__":
    main()
