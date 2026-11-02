#!/usr/bin/env python3

import csv
import os
from pathlib import Path
from collections import defaultdict

PROJECT = Path(
    "/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project"
)

SIM = Path(os.environ["SIM"])
SIM1_MAN = Path(os.environ["SIM1_MAN"])

HARM = (
    PROJECT
    / "08_harmonized_calls"
    / "simulation1"
    / "BICseq2"
)

OUT = (
    PROJECT
    / "09_metrics"
    / "simulation1"
    / "BICseq2"
)

OUT.mkdir(parents=True, exist_ok=True)

CALLS_FILE = HARM / "BICseq2_calls_log2_0.2.tsv"
TASKS_FILE = SIM1_MAN / "bicseq2_tasks.tsv"

TRUTH_DIR = (
    SIM
    / "healthy_cfDNA_spikein"
    / "truth"
)

THRESHOLDS = [0.10, 0.30, 0.50]
PRIMARY = 0.50


# ============================================================
# General helpers
# ============================================================

def norm_chr(x):
    x = str(x).strip()
    if x.lower().startswith("chr"):
        x = x[3:]
    return x


def norm_type(x):
    x = str(x).strip().upper()

    if x in {"GAIN", "DUP", "DUPLICATION", "AMP", "AMPLIFICATION"}:
        return "GAIN"

    if x in {"LOSS", "DEL", "DELETION"}:
        return "LOSS"

    return x


def safe_div(a, b):
    return a / b if b else 0.0


def calc_f1(p, r):
    return (
        2 * p * r / (p + r)
        if (p + r) > 0
        else 0.0
    )


def mean(values):
    values = [
        float(x)
        for x in values
        if x not in ("", None)
    ]

    return (
        sum(values) / len(values)
        if values
        else ""
    )


def overlap_bp(a1, a2, b1, b2):
    return max(
        0,
        min(a2, b2) - max(a1, b1)
    )


def reciprocal_overlap(a1, a2, b1, b2):

    ov = overlap_bp(a1, a2, b1, b2)

    len_a = max(1, a2 - a1)
    len_b = max(1, b2 - b1)

    return (
        ov,
        ov / len_a,
        ov / len_b,
    )


# ============================================================
# Read the 684 expected BIC-seq2 configurations
# ============================================================

with TASKS_FILE.open() as fh:
    tasks = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )

if len(tasks) != 684:
    raise SystemExit(
        f"ERROR: expected 684 configurations, found {len(tasks)}"
    )

print(f"Configurations found: {len(tasks)}")


# ============================================================
# Read exact Simulation-1 truth
#
# Exact format:
#
# #chrom start end cnv_type copy_number
# size_class cnv_id cnv_size
#
# 30 truth CNVs per replicate.
# ============================================================

truth_by_rep = {}

for rep in range(1, 6):

    path = TRUTH_DIR / f"rep{rep}.truth.tsv"

    if not path.exists():
        raise SystemExit(
            f"ERROR: truth file missing: {path}"
        )

    events = []

    with path.open() as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in reader:

            # Header is "#chrom", so use it explicitly.
            chrom = norm_chr(row["#chrom"])

            start = int(row["start"])
            end = int(row["end"])

            cnv_type = norm_type(
                row["cnv_type"]
            )

            copy_number = int(
                row["copy_number"]
            )

            size_class = (
                row["size_class"]
                .strip()
                .lower()
            )

            cnv_id = row["cnv_id"].strip()

            cnv_size = int(
                row["cnv_size"]
            )

            events.append({
                "replicate": rep,
                "event_id": cnv_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "type": cnv_type,
                "copy_number": copy_number,
                "size_class": size_class,
                "size_bp": cnv_size,
                "truth_file": str(path),
            })

    if len(events) != 30:
        raise SystemExit(
            f"ERROR: rep{rep}: expected 30 truth CNVs, "
            f"found {len(events)}"
        )

    gains = sum(
        e["type"] == "GAIN"
        for e in events
    )

    losses = sum(
        e["type"] == "LOSS"
        for e in events
    )

    print(
        f"Truth rep{rep}: "
        f"{len(events)} events "
        f"({gains} gain, {losses} loss)"
    )

    truth_by_rep[rep] = events


# ============================================================
# Read harmonized BIC-seq2 calls
# ============================================================

if not CALLS_FILE.exists():
    raise SystemExit(
        f"ERROR: harmonized calls missing:\n{CALLS_FILE}"
    )

calls_by_config = defaultdict(list)

with CALLS_FILE.open() as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        sample = row["sample_id"]
        bin_bp = str(row["bin_bp"])

        state = norm_type(
            row["state_log2_0.2"]
        )

        if state not in {"GAIN", "LOSS"}:
            continue

        start = int(
            float(row["start_0based"])
        )

        end = int(
            float(row["end"])
        )

        calls_by_config[
            (sample, bin_bp)
        ].append({
            "chrom": norm_chr(row["chrom"]),
            "start": start,
            "end": end,
            "type": state,
            "size_bp": end - start,
            "log2_ratio": row.get(
                "log2_ratio",
                ""
            ),
            "source_file": row.get(
                "source_file",
                ""
            ),
        })


# ============================================================
# One-to-one reciprocal-overlap matching
#
# Requirements:
#   same chromosome
#   same CNV type
#   reciprocal overlap >= threshold
#
# Candidate ranking:
#   highest minimum reciprocal overlap
#   then highest absolute overlap
# ============================================================

def match_events(truth, calls, threshold):

    candidates = []

    for ti, t in enumerate(truth):

        for ci, c in enumerate(calls):

            if t["chrom"] != c["chrom"]:
                continue

            if t["type"] != c["type"]:
                continue

            ov, truth_frac, call_frac = (
                reciprocal_overlap(
                    t["start"],
                    t["end"],
                    c["start"],
                    c["end"],
                )
            )

            rec = min(
                truth_frac,
                call_frac
            )

            if rec >= threshold:

                candidates.append({
                    "truth_index": ti,
                    "call_index": ci,
                    "overlap_bp": ov,
                    "truth_overlap_fraction":
                        truth_frac,
                    "call_overlap_fraction":
                        call_frac,
                    "minimum_reciprocal_overlap":
                        rec,
                })

    candidates.sort(
        key=lambda x: (
            x["minimum_reciprocal_overlap"],
            x["overlap_bp"],
        ),
        reverse=True,
    )

    used_truth = set()
    used_calls = set()

    matches = []

    for m in candidates:

        ti = m["truth_index"]
        ci = m["call_index"]

        if (
            ti in used_truth
            or ci in used_calls
        ):
            continue

        used_truth.add(ti)
        used_calls.add(ci)

        t = truth[ti]
        c = calls[ci]

        start_error = abs(
            c["start"] - t["start"]
        )

        end_error = abs(
            c["end"] - t["end"]
        )

        size_error = abs(
            c["size_bp"] - t["size_bp"]
        )

        m = dict(m)

        m.update({
            "start_error_bp":
                start_error,

            "end_error_bp":
                end_error,

            "mean_breakpoint_error_bp":
                (start_error + end_error) / 2,

            "size_error_bp":
                size_error,

            "size_error_fraction":
                safe_div(
                    size_error,
                    t["size_bp"]
                ),
        })

        matches.append(m)

    return (
        matches,
        used_truth,
        used_calls,
    )


# ============================================================
# Output containers
# ============================================================

config_rows = []

truth_rows = []
call_rows = []
match_rows = []

negative_rows = []


# ============================================================
# Score every expected configuration
# ============================================================

for task in tasks:

    sample = task["sample_id"]
    bin_bp = str(task["bin_bp"])

    background = task["background"]
    depth = task["target_depth"]
    sample_class = task["class"].lower()

    rep_raw = task["replicate"]

    if sample_class == "positive":
        replicate = int(rep_raw)
        truth = truth_by_rep[replicate]
    else:
        replicate = rep_raw
        truth = []

    calls = calls_by_config.get(
        (sample, bin_bp),
        []
    )

    threshold_info = {}

    threshold_truth_maps = {}
    threshold_call_maps = {}

    for threshold in THRESHOLDS:

        matches, used_t, used_c = (
            match_events(
                truth,
                calls,
                threshold,
            )
        )

        tp = len(matches)
        fp = len(calls) - tp
        fn = len(truth) - tp

        precision = safe_div(
            tp,
            tp + fp
        )

        recall = safe_div(
            tp,
            tp + fn
        )

        threshold_info[threshold] = {
            "matches": matches,
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "precision": precision,
            "recall": recall,
            "F1": calc_f1(
                precision,
                recall
            ),
            "FDR": safe_div(
                fp,
                tp + fp
            ),
        }

        threshold_truth_maps[threshold] = {
            m["truth_index"]
            for m in matches
        }

        threshold_call_maps[threshold] = {
            m["call_index"]
            for m in matches
        }


    primary = threshold_info[PRIMARY]
    primary_matches = primary["matches"]

    config_rows.append({
        "tool": "BICseq2",

        "sample_id": sample,
        "background": background,
        "target_depth": depth,
        "replicate": replicate,
        "class": sample_class,
        "bin_bp": int(bin_bp),

        "truth_events": len(truth),
        "called_events": len(calls),

        "TP_50": primary["TP"],
        "FP_50": primary["FP"],
        "FN_50": primary["FN"],
        "precision_50":
            primary["precision"],
        "recall_50":
            primary["recall"],
        "sensitivity_50":
            primary["recall"],
        "F1_50":
            primary["F1"],
        "FDR_50":
            primary["FDR"],

        "TP_30":
            threshold_info[0.30]["TP"],
        "FP_30":
            threshold_info[0.30]["FP"],
        "FN_30":
            threshold_info[0.30]["FN"],
        "precision_30":
            threshold_info[0.30]["precision"],
        "recall_30":
            threshold_info[0.30]["recall"],
        "F1_30":
            threshold_info[0.30]["F1"],

        "TP_10":
            threshold_info[0.10]["TP"],
        "FP_10":
            threshold_info[0.10]["FP"],
        "FN_10":
            threshold_info[0.10]["FN"],
        "precision_10":
            threshold_info[0.10]["precision"],
        "recall_10":
            threshold_info[0.10]["recall"],
        "F1_10":
            threshold_info[0.10]["F1"],

        "mean_breakpoint_error_bp_50":
            mean([
                m[
                    "mean_breakpoint_error_bp"
                ]
                for m in primary_matches
            ]),

        "mean_size_error_bp_50":
            mean([
                m["size_error_bp"]
                for m in primary_matches
            ]),

        "mean_size_error_fraction_50":
            mean([
                m["size_error_fraction"]
                for m in primary_matches
            ]),

        "negative_no_call":
            (
                1
                if sample_class != "positive"
                and len(calls) == 0
                else 0
            ),
    })


    # --------------------------------------------------------
    # Truth-event table
    # --------------------------------------------------------

    primary_by_truth = {
        m["truth_index"]: m
        for m in primary_matches
    }

    for ti, t in enumerate(truth):

        m = primary_by_truth.get(ti)

        truth_rows.append({
            "tool": "BICseq2",
            "sample_id": sample,
            "background": background,
            "target_depth": depth,
            "replicate": replicate,
            "bin_bp": int(bin_bp),

            "truth_event_id":
                t["event_id"],

            "chrom":
                t["chrom"],

            "start":
                t["start"],

            "end":
                t["end"],

            "cnv_type":
                t["type"],

            "copy_number":
                t["copy_number"],

            "size_class":
                t["size_class"],

            "cnv_size_bp":
                t["size_bp"],

            "detected_10":
                int(
                    ti
                    in threshold_truth_maps[0.10]
                ),

            "detected_30":
                int(
                    ti
                    in threshold_truth_maps[0.30]
                ),

            "detected_50":
                int(
                    ti
                    in threshold_truth_maps[0.50]
                ),

            "status_50":
                (
                    "MATCHED"
                    if m
                    else "MISSED"
                ),

            "overlap_bp_50":
                (
                    m["overlap_bp"]
                    if m
                    else 0
                ),

            "truth_overlap_fraction_50":
                (
                    m[
                        "truth_overlap_fraction"
                    ]
                    if m
                    else 0
                ),

            "call_overlap_fraction_50":
                (
                    m[
                        "call_overlap_fraction"
                    ]
                    if m
                    else 0
                ),

            "minimum_reciprocal_overlap_50":
                (
                    m[
                        "minimum_reciprocal_overlap"
                    ]
                    if m
                    else 0
                ),

            "breakpoint_error_bp_50":
                (
                    m[
                        "mean_breakpoint_error_bp"
                    ]
                    if m
                    else ""
                ),

            "size_error_bp_50":
                (
                    m["size_error_bp"]
                    if m
                    else ""
                ),
        })


    # --------------------------------------------------------
    # Called-event table
    # --------------------------------------------------------

    primary_by_call = {
        m["call_index"]: m
        for m in primary_matches
    }

    for ci, c in enumerate(calls):

        m = primary_by_call.get(ci)

        call_rows.append({
            "tool": "BICseq2",

            "sample_id": sample,
            "background": background,
            "target_depth": depth,
            "replicate": replicate,
            "class": sample_class,
            "bin_bp": int(bin_bp),

            "chrom": c["chrom"],
            "start": c["start"],
            "end": c["end"],
            "cnv_type": c["type"],
            "call_size_bp": c["size_bp"],
            "log2_ratio": c["log2_ratio"],

            "matched_10":
                int(
                    ci
                    in threshold_call_maps[0.10]
                ),

            "matched_30":
                int(
                    ci
                    in threshold_call_maps[0.30]
                ),

            "matched_50":
                int(
                    ci
                    in threshold_call_maps[0.50]
                ),

            "status_50":
                (
                    "TP"
                    if m
                    else "FP"
                ),

            "source_file":
                c["source_file"],
        })


    # --------------------------------------------------------
    # Detailed 50% matches
    # --------------------------------------------------------

    for m in primary_matches:

        t = truth[
            m["truth_index"]
        ]

        c = calls[
            m["call_index"]
        ]

        match_rows.append({
            "tool": "BICseq2",

            "sample_id": sample,
            "background": background,
            "target_depth": depth,
            "replicate": replicate,
            "bin_bp": int(bin_bp),

            "truth_event_id":
                t["event_id"],

            "truth_chrom":
                t["chrom"],

            "truth_start":
                t["start"],

            "truth_end":
                t["end"],

            "truth_cnv_type":
                t["type"],

            "truth_copy_number":
                t["copy_number"],

            "truth_size_class":
                t["size_class"],

            "truth_size_bp":
                t["size_bp"],

            "call_chrom":
                c["chrom"],

            "call_start":
                c["start"],

            "call_end":
                c["end"],

            "call_cnv_type":
                c["type"],

            "call_size_bp":
                c["size_bp"],

            "call_log2_ratio":
                c["log2_ratio"],

            "overlap_bp":
                m["overlap_bp"],

            "truth_overlap_fraction":
                m[
                    "truth_overlap_fraction"
                ],

            "call_overlap_fraction":
                m[
                    "call_overlap_fraction"
                ],

            "minimum_reciprocal_overlap":
                m[
                    "minimum_reciprocal_overlap"
                ],

            "start_error_bp":
                m["start_error_bp"],

            "end_error_bp":
                m["end_error_bp"],

            "mean_breakpoint_error_bp":
                m[
                    "mean_breakpoint_error_bp"
                ],

            "size_error_bp":
                m["size_error_bp"],

            "size_error_fraction":
                m["size_error_fraction"],
        })


    # --------------------------------------------------------
    # Negative controls
    # --------------------------------------------------------

    if sample_class != "positive":

        negative_rows.append({
            "tool": "BICseq2",

            "sample_id": sample,
            "background": background,
            "target_depth": depth,
            "replicate": replicate,
            "bin_bp": int(bin_bp),

            "false_calls":
                len(calls),

            "false_altered_bp":
                sum(
                    c["size_bp"]
                    for c in calls
                ),

            "no_false_call":
                int(
                    len(calls) == 0
                ),
        })


# ============================================================
# Writers
# ============================================================

def write_tsv(path, rows):

    if not rows:
        path.write_text("")
        return

    with path.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=list(
                rows[0].keys()
            ),
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)


write_tsv(
    OUT / "BICseq2_configuration_metrics.tsv",
    config_rows,
)

write_tsv(
    OUT / "BICseq2_event_matches_50pct.tsv",
    match_rows,
)

write_tsv(
    OUT / "BICseq2_truth_event_detection.tsv",
    truth_rows,
)

write_tsv(
    OUT / "BICseq2_called_events.tsv",
    call_rows,
)

write_tsv(
    OUT / "BICseq2_negative_controls.tsv",
    negative_rows,
)


# Explicit matched/unmatched tables

write_tsv(
    OUT / "BICseq2_truth_matched_50pct.tsv",
    [
        r
        for r in truth_rows
        if r["detected_50"] == 1
    ],
)

write_tsv(
    OUT / "BICseq2_truth_unmatched_50pct.tsv",
    [
        r
        for r in truth_rows
        if r["detected_50"] == 0
    ],
)

write_tsv(
    OUT / "BICseq2_calls_matched_50pct.tsv",
    [
        r
        for r in call_rows
        if r["matched_50"] == 1
    ],
)

write_tsv(
    OUT / "BICseq2_calls_unmatched_50pct.tsv",
    [
        r
        for r in call_rows
        if r["matched_50"] == 0
    ],
)


# ============================================================
# Configuration-level summaries
# ============================================================

def summarize_configs(
    group_cols,
    filename,
):

    groups = defaultdict(list)

    for row in config_rows:

        key = tuple(
            row[c]
            for c in group_cols
        )

        groups[key].append(row)

    output = []

    for key, rows in groups.items():

        tp = sum(
            int(r["TP_50"])
            for r in rows
        )

        fp = sum(
            int(r["FP_50"])
            for r in rows
        )

        fn = sum(
            int(r["FN_50"])
            for r in rows
        )

        p = safe_div(
            tp,
            tp + fp
        )

        rec = safe_div(
            tp,
            tp + fn
        )

        out = {
            c: v
            for c, v in zip(
                group_cols,
                key
            )
        }

        out.update({
            "n_configurations":
                len(rows),

            "TP_50": tp,
            "FP_50": fp,
            "FN_50": fn,

            "precision_50": p,
            "recall_50": rec,
            "sensitivity_50": rec,
            "F1_50":
                calc_f1(p, rec),

            "FDR_50":
                safe_div(
                    fp,
                    tp + fp
                ),

            "mean_breakpoint_error_bp_50":
                mean([
                    r[
                        "mean_breakpoint_error_bp_50"
                    ]
                    for r in rows
                ]),

            "mean_size_error_bp_50":
                mean([
                    r[
                        "mean_size_error_bp_50"
                    ]
                    for r in rows
                ]),
        })

        output.append(out)

    write_tsv(
        OUT / filename,
        sorted(
            output,
            key=lambda x: tuple(
                str(x[c])
                for c in group_cols
            )
        ),
    )


summarize_configs(
    ["target_depth"],
    "BICseq2_by_depth.tsv",
)

summarize_configs(
    ["bin_bp"],
    "BICseq2_by_bin.tsv",
)

summarize_configs(
    ["target_depth", "bin_bp"],
    "BICseq2_by_depth_bin.tsv",
)

summarize_configs(
    ["background"],
    "BICseq2_by_background.tsv",
)

summarize_configs(
    ["replicate"],
    "BICseq2_by_replicate.tsv",
)


# ============================================================
# Truth-event recall summaries
# ============================================================

def summarize_truth(
    group_cols,
    filename,
):

    groups = defaultdict(list)

    for row in truth_rows:

        key = tuple(
            row[c]
            for c in group_cols
        )

        groups[key].append(row)

    output = []

    for key, rows in groups.items():

        n = len(rows)

        d10 = sum(
            int(r["detected_10"])
            for r in rows
        )

        d30 = sum(
            int(r["detected_30"])
            for r in rows
        )

        d50 = sum(
            int(r["detected_50"])
            for r in rows
        )

        out = {
            c: v
            for c, v in zip(
                group_cols,
                key
            )
        }

        out.update({
            "truth_event_instances": n,

            "detected_10": d10,
            "missed_10": n - d10,
            "recall_10":
                safe_div(d10, n),

            "detected_30": d30,
            "missed_30": n - d30,
            "recall_30":
                safe_div(d30, n),

            "detected_50": d50,
            "missed_50": n - d50,
            "recall_50":
                safe_div(d50, n),
        })

        output.append(out)

    write_tsv(
        OUT / filename,
        sorted(
            output,
            key=lambda x: tuple(
                str(x[c])
                for c in group_cols
            )
        ),
    )


summarize_truth(
    ["cnv_size_bp"],
    "BICseq2_by_CNV_size_bp.tsv",
)

summarize_truth(
    ["size_class"],
    "BICseq2_by_CNV_size_class.tsv",
)

summarize_truth(
    ["cnv_type"],
    "BICseq2_by_CNV_type.tsv",
)

summarize_truth(
    ["target_depth", "cnv_size_bp"],
    "BICseq2_by_depth_CNV_size.tsv",
)

summarize_truth(
    ["bin_bp", "cnv_size_bp"],
    "BICseq2_by_bin_CNV_size.tsv",
)

summarize_truth(
    ["target_depth", "cnv_type"],
    "BICseq2_by_depth_CNV_type.tsv",
)

summarize_truth(
    ["bin_bp", "cnv_type"],
    "BICseq2_by_bin_CNV_type.tsv",
)


# ============================================================
# Negative-control summaries
# ============================================================

def summarize_negative(
    group_cols,
    filename,
):

    groups = defaultdict(list)

    for row in negative_rows:

        key = tuple(
            row[c]
            for c in group_cols
        )

        groups[key].append(row)

    output = []

    for key, rows in groups.items():

        out = {
            c: v
            for c, v in zip(
                group_cols,
                key
            )
        }

        out.update({
            "n_negative_configurations":
                len(rows),

            "total_false_calls":
                sum(
                    int(r["false_calls"])
                    for r in rows
                ),

            "mean_false_calls":
                safe_div(
                    sum(
                        int(r["false_calls"])
                        for r in rows
                    ),
                    len(rows),
                ),

            "total_false_altered_bp":
                sum(
                    int(r["false_altered_bp"])
                    for r in rows
                ),

            "negative_configs_with_no_calls":
                sum(
                    int(r["no_false_call"])
                    for r in rows
                ),

            "negative_no_call_rate":
                safe_div(
                    sum(
                        int(r["no_false_call"])
                        for r in rows
                    ),
                    len(rows),
                ),
        })

        output.append(out)

    write_tsv(
        OUT / filename,
        sorted(
            output,
            key=lambda x: tuple(
                str(x[c])
                for c in group_cols
            )
        ),
    )


summarize_negative(
    ["target_depth"],
    "BICseq2_negative_by_depth.tsv",
)

summarize_negative(
    ["bin_bp"],
    "BICseq2_negative_by_bin.tsv",
)

summarize_negative(
    ["target_depth", "bin_bp"],
    "BICseq2_negative_by_depth_bin.tsv",
)


# ============================================================
# Final validation / summary
# ============================================================

positive_configs = sum(
    r["class"] == "positive"
    for r in config_rows
)

negative_configs = sum(
    r["class"] != "positive"
    for r in config_rows
)

matched_truth = sum(
    int(r["detected_50"])
    for r in truth_rows
)

missed_truth = (
    len(truth_rows)
    - matched_truth
)

matched_calls = sum(
    int(r["matched_50"])
    for r in call_rows
)

unmatched_calls = (
    len(call_rows)
    - matched_calls
)


summary = f"""BIC-seq2 Simulation-1 scoring summary

Expected configurations: 684
Observed configurations: {len(config_rows)}

Positive configurations: {positive_configs}
Negative configurations: {negative_configs}

Truth event instances: {len(truth_rows)}
Expected truth event instances: {95 * 6 * 30}

Called CNV events: {len(call_rows)}

Primary reciprocal-overlap threshold: 50%
Sensitivity thresholds also scored: 10%, 30%

Matched truth events (50%): {matched_truth}
Missed truth events (50%): {missed_truth}

Matched BIC-seq2 calls (50%): {matched_calls}
Unmatched BIC-seq2 calls / FP (50%): {unmatched_calls}

Detailed output directory:
{OUT}
"""

(
    OUT
    / "BICseq2_SCORING_SUMMARY.txt"
).write_text(summary)


print()
print("============================================")
print("BIC-seq2 Simulation-1 scoring COMPLETE")
print("============================================")
print(summary)
