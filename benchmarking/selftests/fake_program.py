#!/usr/bin/env python3
"""Tiny external process used to exercise the generic benchmark executor."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=(
            "success",
            "nonzero",
            "fail-second",
            "sleep",
            "child",
            "detached",
            "bad",
        ),
    )
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--expected-cwd", required=True)
    parser.add_argument("--write-samples", choices=("0", "1"), required=True)
    options = parser.parse_args()
    invocation_log = "fake_invocations.log"
    try:
        with open(invocation_log, "r", encoding="utf-8") as stream:
            invocation = sum(1 for _ in stream) + 1
    except FileNotFoundError:
        invocation = 1
    with open(invocation_log, "a", encoding="utf-8") as stream:
        stream.write(f"{os.getpid()}\n")

    print(f"fake stdout mode={options.mode} invocation={invocation}")
    print("fake stderr", file=sys.stderr)
    if options.mode == "sleep":
        time.sleep(5)
    if options.mode == "child":
        read_fd, write_fd = os.pipe()
        subprocess.Popen(
            [
                sys.executable,
                "-c",
                "import os,signal,sys,time; "
                "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
                "fd=int(sys.argv[1]); os.write(fd,b'R'); os.close(fd); time.sleep(5)",
                str(write_fd),
            ],
            pass_fds=(write_fd,),
        )
        os.close(write_fd)
        ready = os.read(read_fd, 1)
        os.close(read_fd)
        if ready != b"R":
            return 8
        time.sleep(5)
    if options.mode == "detached":
        subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(5)"],
            start_new_session=True,
        )
        time.sleep(5)
    if options.mode == "nonzero":
        return 7
    if options.mode == "fail-second" and invocation == 2:
        return 7
    if options.mode == "bad":
        print("not a tagged result")
        return 0
    working_directory = os.path.realpath(os.getcwd())
    if working_directory != os.path.realpath(options.expected_cwd):
        print("fake program received the wrong working directory", file=sys.stderr)
        return 9
    if options.write_samples == "1":
        with open("field_samples.csv", "w", encoding="utf-8") as stream:
            stream.write(f"x,y,z,value\n0,0,0,{invocation}\n")
    print(
        "ADS_BENCHMARK_RESULT "
        + json.dumps(
            {
                "checksum": "fake-ok",
                "invocation": invocation,
                "physical_step_wall_seconds": round(
                    0.20 + invocation * 0.01, 2
                ),
                "steps": options.steps,
                "working_directory": working_directory,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
