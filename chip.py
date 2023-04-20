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
    angular_offset: float
    pa_offset: float

# Arrays from iwcs.cl
xz = [0,0,0,0,2049,2049,2049,2049]
yz = [4097,4097,4097,4097,0,0,0,0]
sxz = [1,1,1,1,-1,-1,-1,-1]
syz = [-1,-1,-1,-1,1,1,1,1]
x0_SITe = [-61.555,-30.070,1.435,32.910,-30.115,-61.610,32.900,1.362]
y0_SITe = [0.250,0.168,0.185,0.135,-61.912,-61.565,-62.070,-62.020]
x0_E2V = [-61.524,-29.681,2.190,34.043,-29.694,-61.545,34.041,2.159]
y0_E2V = [2.014,1.993,2.0,1.964,-59.982,-59.965,-60.018,-59.978]

#Determined these values using test WCS object
angular_offsets = [0.11439789554022421, 0.08568281915726744, 0.0672471279094002, 0.12289945235938936, 0.06388453802234773, 0.13484887345132082, 0.11316830018295365, 0.08997800570772167] # Degrees
pa_offsets = [-2.0539729107877136, -2.7068906459672206, 2.7913314336967727, 2.104049211290143, -0.7075454619615309, -0.8661395815991457, 1.087876384556838, 0.37682784738123765] # radians

chips = {
    f'c{i+1}': Chip(
    sxz[i], syz[i], x0_SITe[i], y0_SITe[i], x0_E2V[i], y0_E2V[i], xz[i], yz[i], angular_offsets[i], pa_offsets[i]
    ) for i in range(len(sxz))
    }
