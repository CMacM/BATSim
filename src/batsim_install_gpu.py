import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path


CUPY_SPECS = {
    12: "cupy-cuda12x>=13.6,<14",
    13: "cupy-cuda13x>=13.6,<14",
}
NUMPY_SPEC = "numpy>=1.26,<2.0"


def _parse_cuda_major(version):
    if version is None:
        return None
    match = re.search(r"(\d+)(?:\.\d+)?", str(version))
    if match is None:
        return None
    return int(match.group(1))


def _run_command(args):
    try:
        return subprocess.run(
            args,
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return None


def _detect_from_env():
    for name in ("BAT_SIM_CUDA_VERSION", "CUDA_VERSION"):
        major = _parse_cuda_major(os.environ.get(name))
        if major is not None:
            return major, name
    return None, None


def _detect_from_nvcc():
    result = _run_command(["nvcc", "--version"])
    if result is None or result.returncode != 0:
        return None, None
    match = re.search(r"release\s+(\d+)(?:\.\d+)?", result.stdout)
    if match is None:
        return None, None
    return int(match.group(1)), "nvcc"


def _detect_from_nvidia_smi():
    result = _run_command(["nvidia-smi"])
    if result is None or result.returncode != 0:
        return None, None
    match = re.search(r"CUDA Version:\s*(\d+)(?:\.\d+)?", result.stdout)
    if match is None:
        return None, None
    return int(match.group(1)), "nvidia-smi"


def _detect_from_cuda_files():
    cuda_root = Path(os.environ.get("CUDA_HOME") or os.environ.get("CUDA_PATH") or "/usr/local/cuda")
    version_json = cuda_root / "version.json"
    if version_json.exists():
        try:
            data = json.loads(version_json.read_text())
            major = _parse_cuda_major(data.get("cuda", {}).get("version"))
            if major is not None:
                return major, str(version_json)
        except (OSError, json.JSONDecodeError):
            pass

    version_txt = cuda_root / "version.txt"
    if version_txt.exists():
        try:
            major = _parse_cuda_major(version_txt.read_text())
            if major is not None:
                return major, str(version_txt)
        except OSError:
            pass

    return None, None


def detect_cuda_major():
    for detector in (
        _detect_from_env,
        _detect_from_nvcc,
        _detect_from_nvidia_smi,
        _detect_from_cuda_files,
    ):
        major, source = detector()
        if major is not None:
            return major, source
    return None, None


def cupy_spec_for_cuda(cuda_major):
    return CUPY_SPECS.get(cuda_major)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Install the CuPy wheel matching the detected CUDA major version."
    )
    parser.add_argument(
        "--cuda",
        choices=("12", "13"),
        help="Override CUDA detection with an explicit CUDA major version.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the pip command without running it.",
    )
    args = parser.parse_args(argv)

    if args.cuda:
        cuda_major = int(args.cuda)
        source = "--cuda"
    else:
        cuda_major, source = detect_cuda_major()

    spec = cupy_spec_for_cuda(cuda_major)
    if spec is None:
        supported = ", ".join(str(version) for version in sorted(CUPY_SPECS))
        detected = "unknown" if cuda_major is None else str(cuda_major)
        print(
            f"Could not select a CuPy package for CUDA {detected}. "
            f"Supported CUDA major versions: {supported}.",
            file=sys.stderr,
        )
        print(
            "Set BAT_SIM_CUDA_VERSION=12 or BAT_SIM_CUDA_VERSION=13, "
            "or pass --cuda 12/13 to override detection.",
            file=sys.stderr,
        )
        return 1

    command = [sys.executable, "-m", "pip", "install", NUMPY_SPEC, spec]
    print(f"Detected CUDA {cuda_major} from {source}; selected {spec}")
    print("Running: " + shlex.join(command))

    if args.dry_run:
        return 0

    return subprocess.call(command)


if __name__ == "__main__":
    raise SystemExit(main())
