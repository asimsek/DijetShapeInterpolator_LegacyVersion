#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from optparse import OptionParser
import ROOT as rt
import os
import sys
import time


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


def find_root_files_top_level(inputdir):
    root_files = []
    try:
        for name in os.listdir(inputdir):
            full = os.path.join(inputdir, name)
            if os.path.isfile(full) and name.endswith(".root"):
                root_files.append(full)
    except OSError as e:
        sys.stderr.write("ERROR: Could not list directory %s: %s\n" % (inputdir, str(e)))
        return []
    root_files.sort()
    return root_files


def parse_mass_from_filename(path):
    try:
        base = os.path.basename(path)
        base = base.replace("-", "_")
        part = base.split("_M_")[-1]
        mass_str = part.split("_")[0]
        return int(mass_str)
    except Exception:
        return None


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-l','--listFile', dest="listFile", type="string",
                      help="Path to dictionary txt file (groups only will be read)")
    parser.add_option('-b','--baseDir', dest="baseDir", default="signalShapes", type="string",
                      help="Base folder containing per-group subfolders (default: signalShapes)")
    parser.add_option('-t','--type', dest="type", default="nom", type="string",
                      help="type of result (e.g., nom, JESup, JESdown, ...)")

    (options, args) = parser.parse_args()

    if not options.listFile:
        parser.error("Please provide the input dictionary file with -l/--listFile")

    if not os.path.isfile(options.listFile):
        sys.stderr.write("ERROR: Input list file does not exist: %s\n" % options.listFile)
        sys.exit(1)

    if not os.path.isdir(options.baseDir):
        sys.stderr.write("ERROR: Base directory does not exist: %s\n" % options.baseDir)
        sys.exit(1)

    TYPE = '' if options.type.lower() == 'nom' else options.type.upper()

    # ROOT in batch mode
    rt.gROOT.SetBatch(True)

    groups = parse_group_names(options.listFile)
    if not groups:
        sys.stderr.write("ERROR: No group names found in: %s\n" % options.listFile)
        sys.exit(1)

    print("\033[1;34m[INFO]\033[0m Found %d group(s) in %s" % (len(groups), options.listFile))

    total_groups_ok = 0
    total_groups_fail = 0

    for group_name in groups:
        group_dir = os.path.join(options.baseDir, group_name)
        mdl = infer_type_from_group(group_name)

        if not os.path.isdir(group_dir):
            print("\033[1;33m[WARN]\033[0m Missing group directory (skip): %s" % group_dir)
            total_groups_fail += 1
            continue

        # Collect ROOT files at top-level of the group folder
        files = find_root_files_top_level(group_dir)
        if not files:
            print("\033[1;33m[WARN]\033[0m No .root files in %s (skip)" % group_dir)
            total_groups_fail += 1
            continue

        # Output dir inside the group
        outDir = os.path.join(group_dir, "combined")
        try:
            os.makedirs(outDir)
            print("\033[1;33m - Created directory:\033[0m %s" % outDir)
        except OSError:
            pass

        # Stable base name for output (use the group name)
        outFileName = group_name.split("MiniAOD")[0] or os.path.splitext(group_name)[0]

        # Load histograms
        histos = []
        mass_points = []

        # Live, same-line listing of masses as we process
        print("\033[1;31mProcessing group ->\033[0m %s" % (group_name))
        print("\033[1;33m  Masses:\033[0m ", end="", flush=True)
        first_print = True

        for fpath in files:
            mass = parse_mass_from_filename(fpath)
            if mass is None:
                # not a mass file; quietly skip
                continue

            # inline mass print
            if not first_print:
                print(", ", end="", flush=True)
            print(mass, end="", flush=True)
            first_print = False
            #time.sleep(1)

            tfileIn = rt.TFile.Open(fpath)
            if not tfileIn or tfileIn.IsZombie():
                print("\n  [skip] Could not open ROOT file:", fpath)
                continue

            hname = 'h_mjj_ratio_%s_M%i' % (mdl, mass)
            h = tfileIn.Get(hname)
            if not h:
                tfileIn.Close()
                continue

            # Rename for output
            h.SetName('h_%s_%s_%s_M%i_WJ' % (mdl, outFileName, TYPE, mass))
            h.SetTitle('Mjj_WJ/M')
            h.SetDirectory(0)
            histos.append(h)
            mass_points.append(mass)
            tfileIn.Close()

        print()  # newline after inline masses

        if not histos:
            print("\033[1;33m[WARN]\033[0m No histograms collected for group: %s" % group_name)
            total_groups_fail += 1
            continue

        # Compose output file path
        if TYPE == '':
            out_path = os.path.join(outDir, 'InputShapes_%s.root' % outFileName)
        else:
            out_path = os.path.join(outDir, 'InputShapes_%s_%s.root' % (outFileName, TYPE))

        tfileOut = rt.TFile.Open(out_path, 'RECREATE')
        if not tfileOut or tfileOut.IsZombie():
            sys.stderr.write("ERROR: Could not create output file: %s\n" % out_path)
            total_groups_fail += 1
            continue

        tfileOut.cd()
        for h in histos:
            h.Write()
        tfileOut.Close()

        print("\033[1;31m  Wrote\033[0m %d \033[1;31mhistograms to:\033[0m %s" % (len(histos), out_path))
        print ("")
        total_groups_ok += 1

    print("\033[1;34m[SUMMARY]\033[0m Groups success: %d, failed: %d" % (total_groups_ok, total_groups_fail))
    if total_groups_fail > 0 and total_groups_ok == 0:
        sys.exit(2)




