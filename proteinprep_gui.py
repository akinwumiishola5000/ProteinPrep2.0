#!/usr/bin/env python3
"""
proteinprep_gui.py

PySimpleGUI wrapper for proteinprep.py CLI with live log streaming.
Supports:
 - hydrogen addition
 - protonation with custom pH (PDB2PQR/PROPKA)
 - chain selection
 - PDBQT conversion
 - Success/Fail messages
"""

import os
import sys
import subprocess
import threading
import PySimpleGUI as sg

def run_cmd_in_thread(args, window):
    """Run subprocess and stream output to GUI."""
    try:
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            window.write_event_value("-OUT-", line.rstrip("\n"))
        proc.wait()
        window.write_event_value("-DONE-", proc.returncode)
    except Exception as e:
        window.write_event_value("-OUT-", f"[ERROR] {e}")
        window.write_event_value("-DONE-", -1)

def build_args(pdb, outdir, remove_waters, remove_hetero, keep_chains,
               add_hydrogens, protonate, ph_value, auto_pdbqt, batch_file):
    """Construct CLI command from GUI inputs."""
    args = [sys.executable, os.path.join(os.path.dirname(__file__), "proteinprep.py")]
    if batch_file:
        args.extend(["--batch-file", batch_file])
    args.append(pdb)
    if outdir:
        args.extend(["--out-dir", outdir])
    if remove_waters is not None:
        args.append("--remove-waters" if remove_waters else "--no-remove-waters")
    if remove_hetero is not None:
        args.append("--remove-hetero" if remove_hetero else "--no-remove-hetero")
    if keep_chains:
        args.extend(["--keep-chains", keep_chains])
    if add_hydrogens:
        args.append("--add-hydrogens")
    if protonate:
        args.append("--protonate")
        if ph_value:
            args.extend(["--ph", str(ph_value)])
    if auto_pdbqt:
        args.append("--auto-pdbqt")
    return args

def main():
    sg.theme("LightGreen")
    layout = [
        [sg.Text("ProteinPrep GUI", font=("Helvetica", 16))],
        [sg.Text("PDB ID or local PDB file:"), sg.Input(key="-PDB-"), sg.FileBrowse(file_types=(("PDB Files","*.pdb"),), target="-PDB-")],
        [sg.Text("Or batch file (one ID/path per line):"), sg.Input(key="-BATCH-"), sg.FileBrowse(file_types=(("Text Files","*.txt"),), target="-BATCH-")],
        [sg.Checkbox("Remove waters", default=True, key="-WATERS-"),
         sg.Checkbox("Remove heteroatoms", default=True, key="-HETERO-")],
        [sg.Text("Keep chains (comma separated):"), sg.Input(key="-CHAINS-", size=(20,1))],
        [sg.Checkbox("Add hydrogens (OpenBabel)", key="-HYDROG-"),
         sg.Checkbox("Protonate (PDB2PQR/PROPKA)", key="-PROT-"),
         sg.Text("pH:"), sg.Input(default_text="7.4", size=(5,1), key="-PH-"),
         sg.Checkbox("Auto convert to PDBQT (OpenBabel)", key="-PDBQT-")],
        [sg.Text("Output folder:"), sg.Input(default_text=".", key="-OUT-"), sg.FolderBrowse(target="-OUT-")],
        [sg.Button("Run"), sg.Button("Exit")],
        [sg.Multiline(size=(100,20), key="-LOG-", autoscroll=True, disabled=True)]
    ]

    window = sg.Window("ProteinPrep GUI", layout, finalize=True)
    thread = None

    while True:
        event, values = window.read(timeout=100)
        if event in (sg.WIN_CLOSED, "Exit"):
            break

        if event == "Run":
            pdb = values["-PDB-"].strip()
            batch = values["-BATCH-"].strip()
            if not pdb and not batch:
                sg.popup_error("Please enter a PDB ID or select a batch file.")
                continue

            outdir = values["-OUT-"].strip() or "."
            ph_value = values["-PH-"].strip()
            try:
                ph_value = float(ph_value)
            except ValueError:
                sg.popup_error("Invalid pH value. Please enter a number.")
                continue

            args = build_args(
                pdb if pdb else batch,
                outdir,
                values["-WATERS-"],
                values["-HETERO-"],
                values["-CHAINS-"].strip(),
                values["-HYDROG-"],
                values["-PROT-"],
                ph_value if values["-PROT-"] else None,
                values["-PDBQT-"],
                batch if batch else None
            )

            window["-LOG-"].update("")
            window["Run"].update(disabled=True)
            thread = threading.Thread(target=run_cmd_in_thread, args=(args, window), daemon=True)
            thread.start()

        elif event == "-OUT-":
            window["-LOG-"].print(values[event])

        elif event == "-DONE-":
            rc = values[event]
            if rc == 0:
                window["-LOG-"].print("\n[Process finished SUCCESSFULLY ✅]")
            else:
                window["-LOG-"].print(f"\n[Process FAILED ❌ | Return code {rc}]")
            window["Run"].update(disabled=False)

    window.close()

if __name__ == "__main__":
    main()
