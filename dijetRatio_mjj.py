#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import ROOT as rt
from optparse import OptionParser


base_cut = "abs(deltaETAjj)<1.1 && abs(etaWJ_j1)<2.5 && abs(etaWJ_j2)<2.5 && pTWJ_j1>60 && pTWJ_j2>30 && IdTight_j1 && IdTight_j2"

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


def parse_input_list(list_path: str):
    groups = {}
    current_group = None

    with open(list_path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            # Group header: e.g. "RSGToQQ_kMpl01_Run3Summer22EE:"
            if line.endswith(":") and ("/" not in line):
                current_group = line[:-1].strip()
                groups[current_group] = {}
                continue

            # Items: "1000: /eos/.../file.root"
            if ":" in line and current_group is not None:
                left, right = line.split(":", 1)
                mass_str = left.strip()
                path = right.strip()
                if not mass_str.isdigit():
                    continue
                mass = int(mass_str)
                groups[current_group][mass] = path

    return groups


def make_output_names(base_outdir: str, group_name: str, input_root_path: str):
    out_dir = os.path.join(base_outdir, group_name)
    os.makedirs(out_dir, exist_ok=True)

    parent_dir = os.path.basename(os.path.dirname(input_root_path))
    tag = parent_dir.split("MiniAOD")[0]  # drop MiniAOD* suffix if present
    out_file = os.path.join(out_dir, f"{tag}_divide.root")

    return out_dir, out_file, tag


def process_one(entry_path: str, mass: int, sample_type: str, out_file: str):
    f = rt.TFile.Open(entry_path)
    if not f or f.IsZombie():
        print(f"\033[1;31m[ERROR]\033[0m Could not open: {entry_path}")
        return False

    tree = f.Get("rootTupleTree/tree")
    if not tree:
        print(f"\033[1;31m[ERROR]\033[0m Missing TTree in: {entry_path}")
        f.Close()
        return False

    hist_name = f"h_mjj_ratio_{sample_type}_M{mass}"
    h = rt.TH1D(hist_name, hist_name, 75, 0.0, 1.5)

    expr = f"mjj/{float(mass):.1f}"
    nfilled = tree.Project(hist_name, expr, base_cut)

    if nfilled <= 0:
        print(f"\033[1;33m[WARN]\033[0m No entries passed selection for mass={mass} in {entry_path}")

    out = rt.TFile.Open(out_file, "RECREATE")
    if not out or out.IsZombie():
        print(f"\033[1;31m[ERROR]\033[0m Could not create output file: {out_file}")
        f.Close()
        return False

    h.Write()
    out.Close()
    f.Close()
    return True


# python3 dijetRatio_mjj.py -l inputLists/inputShapes_RSGToQQ_kMpl01.txt -o SignalShapes

def main():
    parser = OptionParser()
    parser.add_option("-l", "--listFile", dest="listFile", type="string",
                      help="Path to input list (dictionary) txt file", metavar="INPUT_LIST")
    parser.add_option("-o", "--outdir", dest="outDir", default="./signalShapes", type="string",
                      help="Base output directory (default: %(default)s)")

    (options, args) = parser.parse_args()

    if not options.listFile:
        parser.error("Please provide -l/--listFile pointing to your input list txt file.")

    list_path = options.listFile
    out_base = options.outDir

    if not os.path.isfile(list_path):
        print(f"\033[1;31m[ERROR]\033[0m Input list file not found: {list_path}")
        sys.exit(1)

    # ROOT in batch mode
    rt.gROOT.SetBatch(True)

    groups = parse_input_list(list_path)
    if not groups:
        print(f"\033[1;31m[ERROR]\033[0m No groups/masses parsed from: {list_path}")
        sys.exit(1)

    print(f"\033[1;34m[INFO]\033[0m Parsed {len(groups)} group(s) from {list_path}")

    total_ok, total_fail = 0, 0

    for group_name, mass_map in groups.items():
        if not mass_map:
            print(f"\033[1;33m[WARN]\033[0m Group '{group_name}' has no entries, skipping.")
            continue

        # Determine TYPE (qq/gg/qg)
        typ = infer_type_from_group(group_name)

        ###########
        print(f"\033[1;31mProcessing ->\033[0m {group_name}")
        print("\033[1;33m              Masses:\033[0m ", end="", flush=True)
        printed_any = False
        for mass in sorted(mass_map.keys()):
            in_path = mass_map[mass]
            out_dir, out_file, tag = make_output_names(out_base, group_name, in_path)

            # live, inline printing
            if printed_any:
                print(", ", end="", flush=True)
            print(mass, end="", flush=True)

            ok = process_one(in_path, mass, typ, out_file)
            if ok:
                total_ok += 1
            else:
                total_fail += 1

            printed_any = True

        print ("")
        ##########

    print(f"\033[1;34m[SUMMARY]\033[0m Success: {total_ok}, Failed: {total_fail}")
    if total_fail > 0:
        sys.exit(2)


if __name__ == "__main__":
    main()

