"""Telescope and instrument constants."""

from dataclasses import dataclass

@dataclass
class Constants:
    """All of the constants that are needed."""
    mscale: float = 67.18
    pixscale: float = 0.202*66.66667/mscale
    iroa_d: float = -46.25
    offx: float = 4240
    offy: float = 4240
    nx_num: float = 8800
    ny_num: float = 8800
