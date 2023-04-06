"""Chip class which stores all the chip constants."""

from dataclasses import dataclass

@dataclass
class Chip:
    """Chip dependent information"""
    sxz: int
    syz: int
    x0_SITe: float
    y0_SITe: float
    x0_E2V: float
    y0_E2V: float
    xz: int
    yz: int


# Arrays from iwcs.cl
xz = [0,0,0,0,2049,2049,2049,2049]
yz = [4097,4097,4097,4097,0,0,0,0]
sxz = [1,1,1,1,-1,-1,-1,-1]
syz = [-1,-1,-1,-1,1,1,1,1]
x0_SITe = [-61.555,-30.070,1.435,32.910,-30.115,-61.610,32.900,1.362]
y0_SITe = [0.250,0.168,0.185,0.135,-61.912,-61.565,-62.070,-62.020]
x0_E2V = [-61.524,-29.681,2.190,34.043,-29.694,-61.545,34.041,2.159]
y0_E2V = [2.014,1.993,2.0,1.964,-59.982,-59.965,-60.018,-59.978]

chips = {
    f'c{i+1}': Chip(
    sxz[i], syz[i], x0_SITe[i], y0_SITe[i], x0_E2V[i], y0_E2V[i], xz[i], yz[i]
    ) for i in range(len(sxz))
    }
