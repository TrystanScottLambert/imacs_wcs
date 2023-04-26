"""Telescope and instrument constants."""

from dataclasses import dataclass

@dataclass
class Constants:
    """All of the constants that are needed."""

    iroa_d: float = -46.25
    offx: int = 4240
    offy: int = 4240
    nx_num: int = 8800
    ny_num: int = 8800
