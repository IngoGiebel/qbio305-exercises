#! /usr/bin/env python3

from os import chdir
from os.path import dirname
from pathlib import Path
from subprocess import run
import sys


def main():
    
    script_path = Path(dirname(Path(sys.argv[0]).resolve()))
    classpath = script_path / "current"
    
    run([
            "java",
            "-cp", classpath,
            "jgi.SplitPairsAndSingles", "rp"
        ] +
        sys.argv[1:])


if __name__ == "__main__":
    main()
