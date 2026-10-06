#!/usr/bin/env python3

import argparse
from collections import Counter
from pathlib import Path

from Bio.PDB import MMCIFParser, NeighborSearch


def main():
    parser = argparse.ArgumentParser(
        description="Find residue-level contacts shared across predicted complexes."
    )
    parser.add_argument("directory", type=Path)
    parser.add_argument("--chain-a", default="A")
    parser.add_argument("--chain-b", default="B")
    parser.add_argument("--cutoff", type=float, default=4.0)
    args = parser.parse_args()

    files = sorted(args.directory.glob("*.cif"))

    if not files:
        raise SystemExit("No CIF files found.")

    parser_cif = MMCIFParser(QUIET=True)
    contacts_a = []
    contacts_b = []

    for path in files:
        structure = parser_cif.get_structure(path.stem, str(path))
        model = structure[0]

        atoms_a = list(model[args.chain_a].get_atoms())
        atoms_b = list(model[args.chain_b].get_atoms())

        search = NeighborSearch(atoms_a + atoms_b)
        atom_pairs = search.search_all(args.cutoff, level="A")

        residue_pairs = set()

        for atom1, atom2 in atom_pairs:
            r1 = atom1.get_parent()
            r2 = atom2.get_parent()
            c1 = r1.get_parent().id
            c2 = r2.get_parent().id

            if c1 == args.chain_a and c2 == args.chain_b:
                residue_pairs.add((r1.id[1], r2.id[1]))
            elif c1 == args.chain_b and c2 == args.chain_a:
                residue_pairs.add((r2.id[1], r1.id[1]))

        a = {x for x, _ in residue_pairs}
        b = {y for _, y in residue_pairs}

        contacts_a.append(a)
        contacts_b.append(b)

        print(
            f"{path.name}: "
            f"{len(a)} chain-A residues, "
            f"{len(b)} chain-B residues, "
            f"{len(residue_pairs)} residue pairs"
        )

    count_a = Counter(x for model in contacts_a for x in model)
    count_b = Counter(x for model in contacts_b for x in model)

    n = len(files)

    print("\nChain A consensus")
    for count in range(n, 0, -1):
        residues = sorted(r for r, c in count_a.items() if c == count)
        if residues:
            print(f"{count}/{n}: {','.join(map(str, residues))}")

    print("\nChain B consensus")
    for count in range(n, 0, -1):
        residues = sorted(r for r, c in count_b.items() if c == count)
        if residues:
            print(f"{count}/{n}: {','.join(map(str, residues))}")


if __name__ == "__main__":
    main()
