# Legacy Interpolator for the Dijet Resonance Shapes (Run 3)

These scripts are designed to run within the CMS software environment (CMSSW) on lxplus.

---

## Installation

1. Log in to lxplus and set up the CMSSW environment:

   ```bash
   cmsrel CMSSW_15_0_13
   cd CMSSW_15_0_13/src
   cmsenv
   ```

2. Clone this repository:

   ```bash
   git clone --recursive https://github.com/asimsek/DijetShapeInterpolator_LegacyVersion DijetShapeInterpolator
   cd DijetShapeInterpolator
   ```


---

## Usage

### 1) Convert mjj distributions into a uniform ratio format

```bash
python3 dijetRatio_mjj.py -l inputLists/inputShapes_RSGToQQ_kMpl01.txt -o SignalShapes
```

### 2) Merge all formatted shapes into one ROOT per group

```bash
python3 addSignalShapes.py -l inputLists/inputShapes_RSGToQQ_kMpl01.txt -b SignalShapes -t nom
```

### 3) Extract shapes from the merged ROOT into a small Python module

```bash
python3 extractShapes.py -l inputLists/inputShapes_RSGToQQ_kMpl01.txt -b SignalShapes
```

### 4) Interpolate shapes on a mass grid

```bash
python3 getResonanceShapes.py -l inputLists/inputShapes_RSGToQQ_kMpl01.txt -b SignalShapes --step 100 --fineBinning
```

---



