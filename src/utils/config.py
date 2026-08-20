"""
Parse configuration files with a simple key: value format.
Lines starting with '#' are treated as comments.
"""

from pathlib import Path
from typing import Dict

def load_config(file_path: Path) -> Dict[str, str]:
    """
    Read a configuration file and return a dictionary of parameters.
    Expected format:
        step:    identification
        i/input:    reads.fastq
        o/output:    results/
    """
    config = {}
    with open(file_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if ':' not in line:
                raise ValueError(f"Invalid config line {line_num}: missing ':' separator")
            key, value = line.split(':', 1)
            key = key.strip()
            value = value.strip()
            config[key] = value
    return config
