"""
Custom help display for combined steps.
This will eventually show a unified help screen describing parameters
shared across steps and those unique to each step.
"""

def print_combined_step_help(step_str: str):
    """
    Display combined help for multiple steps.
    Currently shows which steps are included and a summary of shared parameters.
    """
    steps = [s.strip() for s in step_str.split(',')]
    step_name_map = {
        '1': 'identification',
        'identification': 'identification',
        '2': 'annotation',
        'annotation': 'annotation',
        '3': 'analysis',
        'analysis': 'analysis',
    }
    resolved = [step_name_map.get(s, s) for s in steps]

    print(f"Steps: {', '.join(resolved)}\n")
    print("Combined help for multiple steps:")
    print("(Full combined help will be implemented later.)")
    print("Shared parameters: --config, --threads, etc.")
    print("Use mirdeep-p3.py <step> -h to see step-specific parameters.\n")
