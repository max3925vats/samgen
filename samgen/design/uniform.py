"""Single-component: every site gets the same molecule."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Uniform:
    component: str

    def label(self, row: int, col: int) -> str:
        return self.component
