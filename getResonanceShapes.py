#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, glob, importlib.util
from argparse import ArgumentParser
from array import array
from ROOT import TFile, TH1D, Math, gROOT
import numpy as np


def infer_type_from_group(group_name: str) -> str:
    g = (group_name or "").strip().lower()
    if g.startswith("rsgtoqq"):
        return "qq"
    if g.startswith("rsgtogg"):
        return "gg"
    if g.startswith("qstar"):
        return "qg"
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


def resolve_shapes_py(base_dir, group):
    combined = os.path.join(base_dir, group, "combined")
    preferred = os.path.join(combined, f"InputShapes_{group}.py")
    if os.path.isfile(preferred):
        return preferred
    matches = sorted(glob.glob(os.path.join(combined, f"InputShapes_{group}_*.py")))
    return matches[0] if matches else None


def load_shapes_module(module_path):
    modname = f"input_shapes_{abs(hash(module_path))}"
    spec = importlib.util.spec_from_file_location(modname, module_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    # quick sanity check
    if not hasattr(mod, "shapes") or not hasattr(mod, "binxcenters"):
        raise RuntimeError(f"{module_path} does not define 'shapes' and 'binxcenters'")

    return mod

class ShapeStorage:
    def __init__(self, shapes, binxcenters):
        self.shapes = shapes
        self.binxcenters = binxcenters

        if len(self.shapes) < 2:
            print("** ERROR: ** Need at least 2 input shapes, %i provided. Aborting." % (len(self.shapes)))
            sys.exit(1)

        nbins = [len(self.binxcenters)]
        for key in self.shapes.keys():
            norm = sum(self.shapes[key])
            if abs(norm - 1.0) > 0.01:
                print("** ERROR: ** Input shape for m =", key, "GeV not normalized to unity. Aborting.")
                sys.exit(3)
            nbins.append(len(self.shapes[key]))
        if len(set(nbins)) > 1:
            print("** ERROR: ** Bin counts differ among inputs/bin centers. Aborting.")
            sys.exit(2)


def LineShapePDF(shapes, mass, histo):
    x = shapes.binxcenters
    if mass in shapes.shapes:
        y = np.array(shapes.shapes[mass])
    else:
        input_masses = sorted(shapes.shapes.keys())
        min_mass = input_masses[0]
        max_mass = input_masses[-1]
        if mass < min_mass:
            print("** WARNING: ** Extrapolating below lowest input mass.")
            ml, mh = input_masses[0], input_masses[1]
        elif mass > max_mass:
            print("** WARNING: ** Extrapolating above highest input mass.")
            ml, mh = input_masses[-2], input_masses[-1]
        else:
            ml = max([m for m in input_masses if m < mass])
            mh = min([m for m in input_masses if m > mass])

        yl = np.array(shapes.shapes[ml])
        yh = np.array(shapes.shapes[mh])
        y = ((yh - yl) / float(mh - ml)) * float(mass - ml) + yl

    interpolator = Math.Interpolator(len(x))
    interpolator.SetData(len(x), array('d', x), array('d', y.tolist()))

    for i in range(1, histo.GetNbinsX() + 1):
        xcenter = histo.GetBinCenter(i) / float(mass)
        if xcenter > shapes.binxcenters[0] and xcenter < shapes.binxcenters[-1]:
            xlow = histo.GetXaxis().GetBinLowEdge(i) / float(mass)
            if xlow < shapes.binxcenters[0]: xlow = shapes.binxcenters[0]
            xhigh = histo.GetXaxis().GetBinUpEdge(i) / float(mass)
            if xhigh > shapes.binxcenters[-1]: xhigh = shapes.binxcenters[-1]
            integral = interpolator.Integ(xlow, xhigh)
            histo.SetBinContent(i, (integral if integral >= 0. else 0.))
        else:
            histo.SetBinContent(i, 0.)
    integral = histo.Integral()
    if integral > 0:
        histo.Scale(1. / integral)



def main():
    usage = (
        "Examples:\n"
        "  python3 getResonanceShapes.py -l inputLists/inputShapes_RSGToQQ_kMpl01.txt -b signalShapes --step 100\n"
        "  python3 getResonanceShapes.py -l inputLists/inputShapes_RSGToQQ_kMpl01.txt -b signalShapes --step 50 --fineBinning\n"
    )

    parser = ArgumentParser(
        description="Interpolate resonance shapes per group by reading group names from a dictionary txt file.",
        epilog=usage
    )

    parser.add_argument("-l", "--listFile", dest="list_file", required=True,
                        help="Dictionary text file (only group names are read).")
    parser.add_argument("-b", "--baseDir", dest="base_dir", default="signalShapes",
                        help="Base directory containing per-group subfolders (default: %(default)s)")
    parser.add_argument("--step", type=int, required=True,
                        help="Mass interval in GeV (e.g. 100). Range is [min_mass, max_mass] per group.")
    parser.add_argument("--fineBinning", dest="fineBinning", default=False, action="store_true",
                        help="Use fine, 1-GeV binning")
    parser.add_argument("--storePDF", dest="storePDF", default=False, action="store_true",
                        help="Also store a 1-GeV-binned PDF")
    parser.add_argument("--storeCDF", dest="storeCDF", default=False, action="store_true",
                        help="Also store a 1-GeV-binned CDF")
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

    # Standard dijet mass binning
    binBoundaries = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325,
        354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687,
        1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509,
        4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 
        10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]

    groups = parse_group_names(args.list_file)
    if not groups:
        sys.stderr.write(f"ERROR: No group names found in: {args.list_file}\n")
        sys.exit(1)

    print(f"\033[1;34m[INFO]\033[0m Found {len(groups)} group(s) in {args.list_file}")

    # Output directory
    outDir = "interpolatedResonanceShapes"
    os.makedirs(outDir, exist_ok=True)

    ok, fail = 0, 0

    for group in groups:
        model = infer_type_from_group(group)
        shapes_py = resolve_shapes_py(args.base_dir, group)
        if not shapes_py:
            print(f"\033[1;33m[WARN]\033[0m Missing InputShapes .py for group '{group}' under {args.base_dir}/{group}/combined/")
            fail += 1
            continue

        try:
            mod = load_shapes_module(shapes_py)
        except Exception as e:
            sys.stderr.write(f"ERROR: Failed to load shapes from '{shapes_py}': {e}\n")
            fail += 1
            continue

        # Validate and prepare shapes
        try:
            # ensure integer keys
            shapes_dict = {int(k): v for k, v in mod.shapes.items()}
            storage = ShapeStorage(shapes_dict, mod.binxcenters)
        except Exception as e:
            sys.stderr.write(f"ERROR: Invalid shapes content in '{shapes_py}': {e}\n")
            fail += 1
            continue

        # Mass list per group: from group’s min..max -- in steps of args.step
        mass_keys = sorted(storage.shapes.keys())
        min_mass, max_mass = mass_keys[0], mass_keys[-1]
        masses = list(range(min_mass, max_mass + 1, args.step))
        if masses[-1] != max_mass:
            masses.append(max_mass)  # ensure inclusion of the endpoint

        # Build output file name mirroring the input file
        output_file = os.path.basename(shapes_py).replace(".py", ".root").replace("InputShapes", "ResonanceShapes")
        output_path = os.path.join(outDir, output_file)

        # Create output ROOT
        output = TFile(output_path, "RECREATE")
        if not output or output.IsZombie():
            sys.stderr.write(f"ERROR: Could not create output: {output_path}\n")
            fail += 1
            continue

        # Live, same-line print of masses
        print(f"\033[1;31mProcessing group ->\033[0m {group}")
        print("\033[1;33m  Masses:\033[0m ", end="", flush=True)
        first = True

        for mass in masses:
            if not first:
                print(", ", end="", flush=True)
            print(mass, end="", flush=True)
            first = False

            histname = f"h_{model}_{int(mass)}"
            if args.fineBinning:
                h_shape = TH1D(histname, f"{model} Resonance Shape", 14000, 0, 14000)
            else:
                h_shape = TH1D(histname, f"{model} Resonance Shape", len(binBoundaries)-1, array('d', binBoundaries))
            h_shape.SetXTitle("Dijet Mass [GeV]")
            h_shape.SetYTitle("Probability")

            LineShapePDF(storage, mass, h_shape)

            output.cd()
            h_shape.Write()

            if args.storePDF or args.storeCDF:
                h_pdf = TH1D(histname + "_pdf", f"{model} Resonance Shape PDF", 14000, 0, 14000)
                h_cdf = TH1D(histname + "_cdf", f"{model} Resonance Shape CDF", 14000, 0, 14000)

                # PDF: uniformly spread each coarse bin content into 1 GeV bins
                for i in range(1, h_shape.GetNbinsX() + 1):
                    bin_min = h_pdf.GetXaxis().FindBin(h_shape.GetXaxis().GetBinLowEdge(i) + 0.5)
                    bin_max = h_pdf.GetXaxis().FindBin(h_shape.GetXaxis().GetBinUpEdge(i) - 0.5)
                    width = max(1, bin_max - bin_min + 1)
                    val = h_shape.GetBinContent(i) / float(width)
                    for b in range(bin_min, bin_max + 1):
                        h_pdf.SetBinContent(b, val)

                # CDF: cumulative of the PDF
                for i in range(1, h_cdf.GetNbinsX() + 1):
                    bin_min = h_pdf.GetXaxis().FindBin(h_cdf.GetXaxis().GetBinLowEdge(i) + 0.5)
                    bin_max = h_pdf.GetXaxis().FindBin(h_cdf.GetXaxis().GetBinUpEdge(i) - 0.5)
                    curr = 0.0
                    for b in range(bin_min, bin_max + 1):
                        curr += h_pdf.GetBinContent(b)
                    prev = h_cdf.GetBinContent(i - 1)
                    h_cdf.SetBinContent(i, prev + curr)

                output.cd()
                if args.storePDF: h_pdf.Write()
                if args.storeCDF: h_cdf.Write()

        print()  # newline after inline masses
        output.Close()
        print(f"\033[1;33m  Wrote:\033[0m {output_path}\n")
        ok += 1

    print(f"\033[1;34m[SUMMARY]\033[0m Groups success: {ok}, failed: {fail}")
    if ok == 0 and fail > 0:
        sys.exit(2)

if __name__ == '__main__':
    main()

