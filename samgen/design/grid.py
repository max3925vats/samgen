"""Two-component placement from a 0/1 pattern file.

The pattern is a whitespace-separated integer grid; integers map to component
keys via `mapping` (e.g. {0: "base", 1: "ligand"}). Out-of-range sites fall
back to the 0 label (padding with the base component).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List


@dataclass
class Grid:
    cells: List[List[int]]
    mapping: Dict[int, str]

    @classmethod
    def from_file(cls, path: str, mapping: Dict) -> "Grid":
        cells = []
        with open(path) as fh:
            for line in fh:
                row = [int(v) for v in line.split()]
                if row:
                    cells.append(row)
        # config YAML keys arrive as str/int; normalise to int keys
        mapping = {int(k): v for k, v in mapping.items()}
        return cls(cells=cells, mapping=mapping)

    def label(self, row: int, col: int) -> str:
        if row < len(self.cells) and col < len(self.cells[row]):
            code = self.cells[row][col]
        else:
            code = 0  # pad out-of-range sites with the base component
        return self.mapping[code]
