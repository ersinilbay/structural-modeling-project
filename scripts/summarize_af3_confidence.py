#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path
from statistics import mean


def main():
    parser = argparse.ArgumentParser(
        description="Summarize AlphaFold 3 confidence metrics across seeds and samples."
    )
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()

    files = list(args.directory.glob(
        "**/seed-*_sample-*/*_summary_confidences.json"
    ))

    def key(path):
        m = re.search(r"seed-(\d+)_sample-(\d+)", str(path))
        return tuple(map(int, m.groups())) if m else (999, 999)

    files = sorted(files, key=key)

    if not files:
        raise SystemExit("No summary confidence files found.")

    iptms = []

    for path in files:
        m = re.search(r"seed-(\d+)_sample-(\d+)", str(path))
        seed, sample = m.groups()

        data = json.loads(path.read_text())
        iptm = data.get("iptm")
        iptms.append(iptm)

        print(
            f"seed={seed} sample={sample} "
            f"ranking={data.get('ranking_score')} "
            f"iptm={iptm} "
            f"ptm={data.get('ptm')} "
            f"chain_ptm={data.get('chain_ptm')}"
        )

    print()
    print("models:", len(files))
    print("mean ipTM:", round(mean(iptms), 3))
    print("ipTM range:", min(iptms), "-", max(iptms))


if __name__ == "__main__":
    main()
