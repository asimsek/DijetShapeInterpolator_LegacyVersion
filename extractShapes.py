#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
from ROOT import TFile, gROOT
import numpy as np
import os
import sys
import glob
import re


def infer_type_from_group(group_name: str) -> str:
    g = (group_name or "").strip().lower()
    if g.startswith("rsgtoqq"):
        return "qq"
    if g.startswith("rsgtogg"):
        return "gg"
    if g.startswith("qstar"):
        return "qg"
    # Safe default if nothing matches
    return "qq"


def parse_group_names(list_path):
    names, seen = [], set()
    with open(list_path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.endswith(":") and ("/" not in line):
                name = line[:-1].strip()
                if name and name not in seen:
                    names.append(name)
                    seen.add(name)
    return names

def resolve_input_file(base_dir, group):
    combined = os.path.join(base_dir, group, "combined")
    preferred = os.path.join(combined, f"InputShapes_{group}.root")
    if os.path.isfile(preferred):
        return preferred

    pattern = os.path.join(combined, f"InputShapes_{group}_*.root")
    matches = sorted(glob.glob(pattern))
    if matches:
        return matches[0]

    return None

_mass_re = re.compile(r"_M(\d+)(?:_|$)")

def extract_shapes_from_file(input_file, tdir="", model="qq", debug=False):
    tf = TFile.Open(input_file)
    if not tf or tf.IsZombie():
        sys.stderr.write(f"ERROR: Could not open ROOT file: {input_file}\n")
        return {}, []

    directory = tf if tdir == "" else tf.Get(tdir)
    if not directory:
        sys.stderr.write(f"ERROR: Could not find directory '{tdir}' in file: {input_file}\n")
        tf.Close()
        return {}, []

    keys = directory.GetListOfKeys()
    if not keys:
        sys.stderr.write(f"ERROR: No keys found in: {input_file}\n")
        tf.Close()
        return {}, []

    shapes = {}
    binxcenters = []
    nEntries = keys.GetEntries()

    for i in range(nEntries):
        key = keys.At(i)
        hName = key.GetName()
        obj = directory.Get(hName)

        # Keep only TH1
        if not obj or not obj.InheritsFrom("TH1"):
            if debug:
                print(f"[debug] Skipping non-TH1 object: {hName}")
            continue

        name_lower = hName.lower()
        if not (name_lower.startswith(f"h_{model}_") or name_lower.startswith(f"h_mjj_ratio_{model}_")):
            if debug:
                print(f"[debug] Skip (model mismatch): {hName}")
            continue

        # Parse mass from histogram name
        m = _mass_re.search(hName)
        if not m:
            if debug:
                print(f"[debug] Could not parse mass from histogram name: {hName}")
            continue

        mass = int(m.group(1))

        if debug:
            print(f"[debug] Extracting shapes for m = {mass} GeV from {hName}")

        # Bin contents (normalize)
        bincontents = [obj.GetBinContent(b) for b in range(1, obj.GetNbinsX() + 1)]

        if not binxcenters:
            binxcenters = [obj.GetBinCenter(b) for b in range(1, obj.GetNbinsX() + 1)]

        norm = np.array(bincontents, dtype=float)
        total = np.sum(norm)
        if total <= 0:
            if debug:
                print(f"[debug] Histogram {hName} has non-positive integral; skipping.")
            continue
        norm /= total

        shapes[mass] = norm.tolist()

    tf.Close()
    return shapes, binxcenters


def write_python_module(out_path, src_root_path, tdir, shapes, binxcenters):
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("# Extracted shapes:\n")
        f.write(f"# Input file: {src_root_path}\n")
        if tdir:
            f.write(f"# Directory: {tdir}\n")
        f.write("\nshapes = {\n")
        for key in sorted(shapes.keys()):
            f.write(f"  {key}: {shapes[key]},\n")
        f.write("}\n\n")
        f.write(f"binxcenters = {binxcenters}\n")



def main():
    usage = (
        "Examples:\n"
        "  python3 extractShapes.py -l inputLists/inputShapes_RSGToQQ_kMpl01.txt -b signalShapes\n"
        "  python3 extractShapes.py -l inputLists/inputShapes_RSGToQQ_kMpl01.txt -b signalShapes --dir mySubDir\n"
    )

    parser = ArgumentParser(
        description="Extract resonance shapes per group by reading group names from a dictionary txt file.",
        epilog=usage
    )
    parser.add_argument("-l", "--listFile", dest="list_file", required=True,
                        help="Dictionary text file (only group names are read).")
    parser.add_argument("-b", "--baseDir", dest="base_dir", default="signalShapes",
                        help="Base directory that contains per-group subfolders (default: %(default)s)")
    parser.add_argument("--dir", dest="dir", default="",
                        help="Path to TDirectory containing histograms (optional; default: top-level)")
    parser.add_argument("--debug", dest="debug", default=False, action="store_true",
                        help="Debug printout")

    args = parser.parse_args()

    if not os.path.isfile(args.list_file):
        sys.stderr.write(f"ERROR: List file not found: {args.list_file}\n")
        sys.exit(1)
    if not os.path.isdir(args.base_dir):
        sys.stderr.write(f"ERROR: Base directory not found: {args.base_dir}\n")
        sys.exit(1)

    gROOT.SetBatch(True)

    groups = parse_group_names(args.list_file)
    if not groups:
        sys.stderr.write(f"ERROR: No group names found in: {args.list_file}\n")
        sys.exit(1)

    print(f"\033[1;34m[INFO]\033[0m Found {len(groups)} group(s) in {args.list_file}")

    ok, fail = 0, 0

    for group in groups:
        model = infer_type_from_group(group)
        root_path = resolve_input_file(args.base_dir, group)
        if not root_path:
            print(f"\033[1;33m[WARN]\033[0m Missing InputShapes for group '{group}' under {args.base_dir}/{group}/combined/")
            fail += 1
            continue

        print(f"\033[1;31mProcessing group ->\033[0m {group}")
        shapes, binx = extract_shapes_from_file(root_path, tdir=args.dir, model=model, debug=args.debug)
        if not shapes:
            print(f"\033[1;33m[WARN]\033[0m No shapes extracted for group: {group}")
            fail += 1
            continue

        out_path = root_path.replace(".root", ".py")
        try:
            write_python_module(out_path, root_path, args.dir, shapes, binx)
            print(f"\033[1;33m  Wrote:\033[0m {out_path}")
            ok += 1
        except Exception as e:
            sys.stderr.write(f"ERROR: Failed to write output file '{out_path}': {e}\n")
            fail += 1

        print("")
    
    print(f"\033[1;34m[SUMMARY]\033[0m Groups success: {ok}, failed: {fail}")
    if ok == 0 and fail > 0:
        sys.exit(2)

if __name__ == '__main__':
    main()

