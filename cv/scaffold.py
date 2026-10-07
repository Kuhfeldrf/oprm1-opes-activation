"""Shared CV definitions (human P35372 numbering; GPCRdb oprm_human segments).

Alignment scaffold = helical cores of TM1-TM5, excluding the cytoplasmic end of
TM5 (moves with TM6 on activation), ICL3, and ECL2 (9MQH/9MQI/9MQJ gap 225-226).
"""
SCAFFOLD = (list(range(72, 97)) + list(range(105, 132)) + list(range(140, 171))
            + list(range(184, 207)) + list(range(228, 256)))
NPXXYA = list(range(334, 340))
CV1_PAIR = (167, 281)          # R3.50 CA - T6.34 CA
HBOND_PAIR = ((167, "CZ"), (281, "OG1"))
