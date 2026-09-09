"""Environment setup and verification utilities."""

import os
import subprocess
import sys

def check_colabfold():
    """Verify ColabFold is installed and accessible."""
    try:
        result = subprocess.run(
            ["colabfold_batch", "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            return True, result.stdout.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    return False, "Not found or not accessible"

def check_pandas():
    """Verify pandas is installed."""
    try:
        import pandas
        return True, f"pandas {pandas.__version__}"
    except ImportError:
        return False, "Not installed"

def check_environment():
    """Run all environment checks."""
    checks = [
        ("ColabFold", check_colabfold),
        ("pandas", check_pandas),
    ]

    results = {}
    for name, check_fn in checks:
        ok, msg = check_fn()
        results[name] = (ok, msg)
        status = "✓" if ok else "✗"
        print(f"{status} {name}: {msg}")

    return all(ok for ok, _ in results.values())
