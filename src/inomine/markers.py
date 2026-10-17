"""Inovirus marker HMM definitions and shared naming helpers."""

import re

RAW_MARKERS = {
    "zot": ["core_Zot_IPR008900"],
    "rstA": ["core_RstA_IPR003491"],
    "rstB": ["core_RstB_IPR010008", "CTX-_RstB_IPR009712"],
    "ace": ["core_Ace_IPR019670"],
    "g3p": ["CTX-_G3P_Self-built", "G2P-2", "G3P_IPR032014"],
    "g8p": ["CTX-_G8P_IPR008020"],
    "ctx": ["CTX+_CtxA_IPR001144", "CTX+_CtxB_IPR001835"],
    "g2p": ["G2P_IPR006516"],
    "g5p": ["G5P_IPR003512"],
    "shared": ["core_RstR_IPR010982", "CTX-I_RstR_IPR021096",
               "CTX-I_XafT_Self-built", "core_RstA_IPR008766"],
}

MARKER_KIND = {name: kind for kind, names in RAW_MARKERS.items() for name in names}

INOVIRUS_SIGNAL = {"zot", "rstA", "rstB", "ace", "g3p", "g8p", "ctx", "g2p", "g5p"}


def gene_name(marker):
    m = re.sub(r"^(core|dif|int|CTX-I|CTX[-+])_", "", marker)
    m = re.sub(r"_(IPR\d+|Self-?built|PHA\d+)$", "", m)
    return {"G2P-2": "G3P"}.get(m, m)
