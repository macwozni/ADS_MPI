#!/usr/bin/env python3
"""Process launcher stub with one explicit, transient control-plane failure."""

from __future__ import annotations

import json
import os
import sys


FAILURE_PREFIX = "ADS_BENCHMARK_LAUNCHER_FAILURE "


def main() -> int:
    try:
        separator = sys.argv.index("--")
    except ValueError:
        print("missing payload separator", file=sys.stderr)
        return 2
    payload = sys.argv[separator + 1 :]
    if not payload:
        print("missing payload", file=sys.stderr)
        return 2

    invocation_log = "fake_launcher_invocations.log"
    try:
        with open(invocation_log, "r", encoding="utf-8") as stream:
            invocation = sum(1 for _ in stream) + 1
    except FileNotFoundError:
        invocation = 1
    with open(invocation_log, "a", encoding="utf-8") as stream:
        stream.write(f"{os.getpid()}\n")

    if invocation == 1:
        print(
            FAILURE_PREFIX
            + json.dumps(
                {
                    "kind": "mpi",
                    "schema_version": 1,
                    "source": "fake-control-plane-launcher",
                },
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        return 86

    os.execvpe(payload[0], payload, os.environ)
    return 127


if __name__ == "__main__":
    raise SystemExit(main())
