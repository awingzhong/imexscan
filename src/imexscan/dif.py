"""Vibrio dif-like att motif scanner.

Python port of find_Vibrio_dif_motif.sh (gawk three-round regex scan).
Outputs the same CSV columns and match records as the original script.
"""

import csv
import gzip
import re

ROUNDS = [
    ("full-length",
     r"[ATCG][TC][TG][ATCG][AT][AC]CATAA[ATCG]{6,12}TAAT[AG][TC]G[ATC]A[TCG][TC]",
     r"[AG][AGC]T[GAT]C[AG][TC]ATTA[ATCG]{6,12}TTATG[GT][AT][ATCG][CA][GA][ATCG]"),
    ("left-armed",
     r"[ATCG][TC][TG][ATCG][AT][AC]CATAA[TC][AG][TC]A[TC]A",
     r"T[AG]T[AG][TC][AG]TTATG[GT][AT][ATCG][CA][GA][ATCG]"),
    ("right-armed",
     r"[TC][AG][TC]A[TC]ATAAT[AG][TC]G[ATC]A[TCG][TC]",
     r"[AG][AGC]T[GAT]C[AG][TC]ATTAT[AG]T[AG][TC][AG]"),
]

CSV_HEADER = ["Sequence_ID", "Match_Type", "Start_Position", "End_Position",
              "Match_Length", "Match_Sequence", "Type_of_DIF-like_att"]

_COMPILED = [(name, re.compile(fwd, re.IGNORECASE), re.compile(rev, re.IGNORECASE))
             for name, fwd, rev in ROUNDS]


def open_maybe_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _scan(seq, regex, blocked):
    """Non-overlapping scan; hits overlapping a blocked range are skipped.

    Mirrors the awk script: the next search starts after each hit, and only
    round-1 (full-length) hits are recorded as blockers.
    """
    hits = []
    pos = 0
    while True:
        m = regex.search(seq, pos)
        if not m:
            break
        s, e = m.start(), m.end()
        if not any(s < be and e > bs for bs, be in blocked):
            hits.append((s, e))
        pos = e
    return hits


def scan(path):
    """Return match records [seq_id, type, start, end, len, seq, round_name] (1-based)."""
    records = []
    seq_id, seq_lines = None, []

    def flush():
        if not seq_lines:
            return
        seq = "".join(seq_lines)
        blocked = []
        for round_i, (round_name, fwd, rev) in enumerate(_COMPILED):
            for regex, match_type in ((fwd, "Forward"), (rev, "Reverse_Complement")):
                for s, e in _scan(seq, regex, blocked):
                    records.append([seq_id if seq_id is not None else "",
                                    match_type, s + 1, e, e - s, seq[s:e], round_name])
                    if round_i == 0:
                        blocked.append((s, e))
        seq_lines.clear()

    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                seq_id = line[1:].rstrip()
            else:
                seq_lines.append(line.strip())
        flush()
    return records


def write_csv(records, out_path):
    with open(out_path, "w", newline="") as fh:
        fh.write(",".join(CSV_HEADER) + "\n")
        writer = csv.writer(fh, quoting=csv.QUOTE_ALL, lineterminator="\n")
        writer.writerows(records)
