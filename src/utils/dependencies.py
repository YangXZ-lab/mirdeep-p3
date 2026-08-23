"""
External software dependency checks.
"""

import shutil
import subprocess
from pathlib import Path
from typing import Dict, Optional, List, Tuple

def check_external_tools(required: Dict[str, Optional[str]]) -> Tuple[bool, List[str]]:
    """
    Check availability of command-line tools by attempting to execute them.
    required: dict mapping tool name (used as default command) -> custom path or None.
    If custom_path is provided:
        - If it contains a directory separator, treat as a literal file path.
        - Otherwise, treat as a command name to be found in PATH.
    Returns (all_ok, missing_list) and prints status.
    """
    missing = []
    print("External dependency check:")
    for tool, custom_path in required.items():
        command = None

        if custom_path:
            p = Path(custom_path)
            # Check if it looks like a plain name (no directory component)
            if p.name == custom_path and p.parent == Path('.'):
                # It's a bare name; find in PATH
                found = shutil.which(custom_path)
                if found:
                    command = custom_path
                else:
                    print(f"  {tool} ✘ ('{custom_path}' not found in PATH)")
                    missing.append(tool)
                    continue
            else:
                # It's a specific path
                if p.is_file() and (p.stat().st_mode & 0o111):
                    command = str(p)
                else:
                    print(f"  {tool} ✘ (specified path not executable: {custom_path})")
                    missing.append(tool)
                    continue
        else:
            # No custom path; search tool name in PATH
            if shutil.which(tool):
                command = tool
            else:
                print(f"  {tool} ✘ (not found in PATH)")
                missing.append(tool)
                continue

        # Test execution (try --version, -version, or empty args)
        success = False
        for test_arg in ['--version', '-version', '']:
            try:
                cmd = [command] + ([test_arg] if test_arg else [])
                proc = subprocess.run(cmd, stdout=subprocess.DEVNULL,
                                      stderr=subprocess.DEVNULL, timeout=20)
                if proc.returncode == 0:
                    success = True
                    break
            except (subprocess.TimeoutExpired, FileNotFoundError, PermissionError, OSError):
                success = False
                break

        if success:
            print(f"  {tool} ✔ ({command})")
        else:
            print(f"  {tool} ✘ (unable to execute, check installation)")
            missing.append(tool)
    return len(missing) == 0, missing
