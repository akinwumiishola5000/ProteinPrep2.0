#!/usr/bin/env python3
"""
proteinprep.py

CLI tool to:
 - fetch PDB by ID (with retries)
 - clean protein: remove waters, remove heteroatoms (optionally keep named ligands)
 - optionally keep only specified chains (A,B,...)
 - optionally add hydrogens (OpenBabel)
 - optionally protonate (PDB2PQR + PROPKA)
 - convert final structure to both PDB and PDBQT
 - produce a JSON log file with a short report

Author: Ishola Abeeb Akinwumi (2.0)
"""
import os
import sys
import time
import json
import shutil
import subprocess
from typing import List, Optional
import requests
import typer

app = typer.Typer(add_completion=False)

# -------------------- Utilities --------------------
def which(program: str) -> Optional[str]:
    return shutil.which(program)

def obabel_available() -> bool:
    return which("obabel") is not None or which("babel") is not None

def pdb2pqr_available() -> bool:
    return which("pdb2pqr") is not None

def try_install_openbabel() -> bool:
    print("[INSTALL] Attempting to install OpenBabel (conda -> pip fallback)...")
    if which("conda"):
        try:
            subprocess.check_call(["conda", "install", "-c", "conda-forge", "openbabel", "-y"])
            return obabel_available()
        except Exception as e:
            print("[WARN] Conda install failed:", e)
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "openbabel-wheel"])
        return obabel_available()
    except Exception as e:
        print("[WARN] pip install openbabel-wheel failed:", e)
    return False

# -------------------- Download PDB --------------------
def download_pdb(pdb_id: str, out_path: str, retries: int = 3, timeout: int = 10):
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            resp = requests.get(url, timeout=timeout)
            resp.raise_for_status()
            with open(out_path, "wb") as fh:
                fh.write(resp.content)
            return
        except requests.RequestException as e:
            last_err = e
            print(f"[DOWNLOAD] Attempt {attempt}/{retries} failed: {e}")
            time.sleep(2)
    raise RuntimeError(f"Failed to download {pdb_id} after {retries} attempts: {last_err}")

# -------------------- Clean PDB --------------------
def clean_pdb(input_pdb: str,
              output_pdb: str,
              remove_waters: bool = True,
              remove_hetero: bool = True,
              keep_chains: Optional[List[str]] = None,
              keep_ligands: Optional[List[str]] = None) -> dict:
    keep_chains = [c.strip().upper() for c in (keep_chains or [])] if keep_chains else None
    keep_ligands = [k.strip().upper() for k in (keep_ligands or [])] if keep_ligands else []
    removed = {"waters": 0, "hetero_residues": 0, "skipped_chains": 0}
    wrote_any = False

    with open(input_pdb, "r") as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM  ", "HETATM", "TER", "END")):
                rec = line[0:6].strip()
                chain_id = line[21:22].strip() if len(line) >= 22 else ""
                resname = line[17:20].strip().upper() if len(line) >= 20 else ""
                if keep_chains and chain_id and chain_id.upper() not in keep_chains:
                    removed["skipped_chains"] += 1
                    continue
                if remove_waters and resname in ("HOH", "H2O", "WAT"):
                    removed["waters"] += 1
                    continue
                if remove_hetero and rec == "HETATM" and resname not in keep_ligands:
                    removed["hetero_residues"] += 1
                    continue
                fout.write(line)
                wrote_any = True
    if not wrote_any:
        raise RuntimeError("No ATOM/HETATM records written — check input or filters.")
    return removed

# -------------------- OpenBabel helpers --------------------
def add_hydrogens_with_obabel(input_file: str, output_file: str) -> dict:
    exe = which("obabel") or which("babel")
    if not exe:
        raise RuntimeError("OpenBabel CLI not found.")
    cmd = [exe, input_file, "-O", output_file, "-h"]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return {"cmd": " ".join(cmd), "rc": proc.returncode, "output": proc.stdout}

def convert_to_pdbqt_with_obabel(input_file: str, output_file: str) -> dict:
    exe = which("obabel") or which("babel")
    if not exe:
        raise RuntimeError("OpenBabel CLI not found.")
    cmd = [exe, input_file, "-O", output_file]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return {"cmd": " ".join(cmd), "rc": proc.returncode, "output": proc.stdout}

# -------------------- PDB2PQR helpers --------------------
def protonate_with_pdb2pqr(input_file: str, output_file: str, ph: float = 7.4) -> dict:
    if not pdb2pqr_available():
        raise RuntimeError("PDB2PQR CLI not found.")
    cmd = ["pdb2pqr", "--ff=PARSE", f"--with-ph={ph}", input_file, output_file]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return {"cmd": " ".join(cmd), "rc": proc.returncode, "output": proc.stdout}

# -------------------- Main workflow --------------------
@app.command()
def main(
    pdb: str = typer.Argument(..., help="PDB ID (4 chars) or local PDB path"),
    out_dir: str = typer.Option(".", help="Output directory"),
    remove_waters: bool = typer.Option(True, help="Remove water molecules"),
    remove_hetero: bool = typer.Option(True, help="Remove heteroatoms"),
    keep_chains: Optional[str] = typer.Option(None, help="Comma-separated chain IDs to keep, e.g. 'A,C'"),
    keep_ligands: Optional[str] = typer.Option(None, help="Comma-separated 3-letter ligand names to keep, e.g. NAD,HEM"),
    add_hydrogens: bool = typer.Option(False, help="Add hydrogens using OpenBabel"),
    protonate: bool = typer.Option(False, help="Protonate using PDB2PQR + PROPKA"),
    auto_pdbqt: bool = typer.Option(True, help="Convert final structure to PDBQT"),
    batch_file: Optional[str] = typer.Option(None, help="Optional: path to a newline-separated file of PDB IDs or paths"),
    ph: float = typer.Option(7.4, help="pH for protonation")
):
    os.makedirs(out_dir, exist_ok=True)
    targets = []

    if batch_file:
        if not os.path.exists(batch_file):
            typer.echo(f"[ERROR] Batch file not found: {batch_file}")
            raise typer.Exit(code=1)
        with open(batch_file, "r") as f:
            targets = [ln.strip() for ln in f if ln.strip()]
    else:
        targets = [pdb]

    keep_chain_list = [c.strip().upper() for c in keep_chains.split(",")] if keep_chains else None
    keep_lig_list = [k.strip().upper() for k in keep_ligands.split(",")] if keep_ligands else None

    obabel_ready = obabel_available()
    pdb2pqr_ready = pdb2pqr_available()

    if (add_hydrogens or auto_pdbqt) and not obabel_ready:
        typer.echo("[WARN] OpenBabel not found, attempting install...")
        try_install_openbabel()
        obabel_ready = obabel_available()

    reports = []

    for t in targets:
        try:
            # Fetch or use local
            if os.path.exists(t):
                fetched = t
                pdb_label = os.path.splitext(os.path.basename(t))[0]
            else:
                pdb_label = t.upper()
                fetched = os.path.join(out_dir, f"{pdb_label}.pdb")
                typer.echo(f"[FETCH] Downloading {pdb_label} ...")
                download_pdb(pdb_label, fetched)

            # Clean
            cleaned = os.path.join(out_dir, f"{pdb_label}_clean.pdb")
            report = {"input": t, "fetched": fetched, "cleaned": cleaned}
            removed = clean_pdb(fetched, cleaned, remove_waters, remove_hetero, keep_chain_list, keep_lig_list)
            report["removed"] = removed

            processed = cleaned

            # Protonation
            if protonate:
                if pdb2pqr_ready:
                    proton_file = os.path.join(out_dir, f"{pdb_label}_final.pdb")
                    typer.echo(f"[PDB2PQR] Protonating {processed} -> {proton_file} at pH {ph}")
                    pqr_res = protonate_with_pdb2pqr(processed, proton_file, ph)
                    report["protonate"] = pqr_res
                    if pqr_res["rc"] == 0:
                        processed = proton_file
                        typer.echo(f"[SUCCESS] Protonation complete: {proton_file}")
                    else:
                        typer.echo("[FAIL] Protonation failed.")
                else:
                    report["protonate_skipped"] = "PDB2PQR not available"

            # Hydrogen addition (if not protonated)
            if add_hydrogens and not protonate:
                if obabel_ready:
                    hydrog_file = os.path.join(out_dir, f"{pdb_label}_final.pdb")
                    typer.echo(f"[OBABEL] Adding hydrogens {processed} -> {hydrog_file}")
                    obres = add_hydrogens_with_obabel(processed, hydrog_file)
                    report["add_hydrogens"] = obres
                    if obres["rc"] == 0:
                        processed = hydrog_file
                        typer.echo(f"[SUCCESS] Hydrogen addition complete: {hydrog_file}")
                    else:
                        typer.echo("[FAIL] Hydrogen addition failed.")
                else:
                    report["add_hydrogens_skipped"] = "OpenBabel not available"

            # PDBQT conversion
            if auto_pdbqt:
                if obabel_ready:
                    pdbqt_file = os.path.join(out_dir, f"{pdb_label}_final.pdbqt")
                    typer.echo(f"[OBABEL] Converting {processed} -> {pdbqt_file}")
                    obres = convert_to_pdbqt_with_obabel(processed, pdbqt_file)
                    report["pdbqt"] = obres
                    if obres["rc"] == 0:
                        typer.echo(f"[SUCCESS] PDBQT conversion complete: {pdbqt_file}")
                    else:
                        typer.echo("[FAIL] PDBQT conversion failed.")
                else:
                    report["pdbqt_skipped"] = "OpenBabel not available"

            reports.append(report)
            typer.echo(f"[DONE] {pdb_label}")

        except Exception as e:
            typer.echo(f"[ERROR] processing {t}: {e}")
            reports.append({"input": t, "error": str(e)})

    logf = os.path.join(out_dir, "proteinprep_log.json")
    with open(logf, "w") as fh:
        json.dump(reports, fh, indent=2)
    typer.echo(f"[FINISHED] Reports: {logf}")

if __name__ == "__main__":
    app()
