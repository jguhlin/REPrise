#!/usr/bin/env python3
"""
Equivalence comparator between C++ and Rust implementations.

Run sequence (from repo root):
  1) make eq-cpp
  2) cargo test --test equiv -- --nocapture
  3) pixi run python eval/compare_equiv.py

What this script does:
  - Executes ./eq_cpp to get the C++ JSON
  - Executes cargo test --test equiv -- --nocapture and extracts the printed JSON
  - Compares both JSONs field-by-field and prints a summary with mismatches

JSON schema expected:
{
  "cpp": {
    "masking_align": [ { "i": int, "consensusstart": int, "consensusend": int, "result": [int,int] }, ... ],
    "chrtracer":      [ { "pos": int, "chr": str, "off": int }, ... ],
    "mask_extention_score": [ { "is_right": bool, "ext": int, "i": int, "band": int, "score": int }, ... ]
  }
}
{
  "rust": {
    "masking_align": [ ... same cases ... ],
    "chrtracer":      [ ... same cases ... ],
    "mask_extention_score": [ ... same cases ... ]
  }
}

Note:
  For now, both harnesses produce placeholder-deterministic logic for masking_align and mask_extention_score
  and a simple chrtracer wrapper. This ensures shape equivalence while exposing the harness. Replace with
  real wired outputs once internals are exported under cfg(test) without changing production semantics.
"""
import json
import re
import subprocess
import sys
from typing import Any, Dict, List, Tuple


def run_cmd(cmd: List[str]) -> Tuple[int, str, str]:
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out, err = p.communicate()
    return p.returncode, out, err


def get_cpp_json() -> Dict[str, Any]:
    rc, out, err = run_cmd(["./eq_cpp"])
    if rc != 0:
        print("ERROR: eq_cpp execution failed")
        print(err)
        sys.exit(1)
    try:
        return json.loads(out)
    except json.JSONDecodeError as e:
        print("ERROR: Failed to decode C++ JSON:", e)
        print("Output was:\n", out)
        sys.exit(1)


def get_rust_json() -> Dict[str, Any]:
    # Run the specific test and capture stdout. We rely on the test printing exactly one JSON object.
    rc, out, err = run_cmd(["cargo", "test", "--test", "equiv", "--", "--nocapture"])
    if rc != 0:
        print("ERROR: cargo test equiv failed")
        print(err)
        sys.exit(1)
    # Extract the last JSON object printed. Use a simple heuristic: find first '{' and last '}'.
    # Cargo prints extra lines, so we attempt to detect a pretty JSON block.
    m = re.search(r"(\{[^]*\})", out)
    if not m:
        # Fallback: try to find lines between braces explicitly
        lines = out.splitlines()
        buf = []
        seen = False
        depth = 0
        for line in lines:
            if "{" in line:
                seen = True
                depth += line.count("{")
            if seen:
                buf.append(line)
            if "}" in line:
                depth -= line.count("}")
                if seen and depth <= 0:
                    break
        if not buf:
            print("ERROR: Could not locate JSON in cargo test output.")
            print("stdout:\n", out)
            print("stderr:\n", err)
            sys.exit(1)
        text = "\n".join(buf)
    else:
        text = m.group(1)
    try:
        return json.loads(text)
    except json.JSONDecodeError as e:
        print("ERROR: Failed to decode Rust JSON:", e)
        print("Captured text was:\n", text)
        sys.exit(1)


def compare_arrays(lhs: List[Dict[str, Any]], rhs: List[Dict[str, Any]], keys: List[str], label: str, mismatches: List[str]) -> None:
    if len(lhs) != len(rhs):
        mismatches.append(f"{label}: length mismatch cpp={len(lhs)} rust={len(rhs)}")
        return
    for i, (a, b) in enumerate(zip(lhs, rhs)):
        for k in keys:
            if a.get(k) != b.get(k):
                mismatches.append(f"{label}[{i}] key {k} mismatch: cpp={a.get(k)} rust={b.get(k)}")


def main():
    cpp = get_cpp_json()
    rust = get_rust_json()
    mismatches: List[str] = []

    # Schema presence
    if "cpp" not in cpp:
        print("ERROR: 'cpp' root key missing")
        sys.exit(1)
    if "rust" not in rust:
        print("ERROR: 'rust' root key missing")
        sys.exit(1)

    # Compare masking_align arrays
    compare_arrays(
        cpp["cpp"].get("masking_align", []),
        rust["rust"].get("masking_align", []),
        ["i", "consensusstart", "consensusend", "result"],
        "masking_align",
        mismatches,
    )

    # Compare chrtracer arrays
    compare_arrays(
        cpp["cpp"].get("chrtracer", []),
        rust["rust"].get("chrtracer", []),
        ["pos", "chr", "off"],
        "chrtracer",
        mismatches,
    )

    # Compare mask_extention_score arrays
    compare_arrays(
        cpp["cpp"].get("mask_extention_score", []),
        rust["rust"].get("mask_extention_score", []),
        ["is_right", "ext", "i", "band", "score"],
        "mask_extention_score",
        mismatches,
    )

    if not mismatches:
        print("EQUIVALENCE PASS: C++ and Rust JSON outputs match for all compared fields.")
        sys.exit(0)
    else:
        print("EQUIVALENCE FAIL: mismatches detected:")
        for m in mismatches:
            print(" -", m)
        # Write equivalence_report.md with brief detail
        with open("equivalence_report.md", "w") as f:
            f.write("# Equivalence Report\n\n")
            f.write("Automated comparison found mismatches between C++ and Rust harness outputs.\n\n")
            f.write("## Summary\n")
            for m in mismatches[:10]:
                f.write(f"- {m}\n")
            f.write("\n## Repro steps\n")
            f.write("1) make eq-cpp\n")
            f.write("2) cargo test --test equiv -- --nocapture\n")
            f.write("3) pixi run python eval/compare_equiv.py\n")
            f.write("\n## References\n")
            f.write("- C++ functions: REPrise.cpp lines ")
            f.write("836 (masking_align), 919 (mask_extention_score), 1072 (chrtracer)\n")
            f.write("- Rust sources: REPrise-rs/src/alg/repeat.rs, REPrise-rs/src/lib.rs\n")
            f.write("\nInitial hypothesis: Schema/placeholder mismatch or mapping differences around padding boundaries.\n")
        sys.exit(2)


if __name__ == "__main__":
    main()