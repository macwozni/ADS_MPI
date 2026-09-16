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
        "--mode", choices=("success", "nonzero", "sleep", "child", "bad")
    )
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--expected-cwd", required=True)
    options = parser.parse_args()
    print(f"fake stdout mode={options.mode}")
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
    if options.mode == "nonzero":
        return 7
    if options.mode == "bad":
        print("not a tagged result")
        return 0
    working_directory = os.path.realpath(os.getcwd())
    if working_directory != os.path.realpath(options.expected_cwd):
        print("fake program received the wrong working directory", file=sys.stderr)
        return 9
    print(
        "ADS_BENCHMARK_RESULT "
        + json.dumps(
            {
                "checksum": "fake-ok",
                "steps": options.steps,
                "working_directory": working_directory,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
