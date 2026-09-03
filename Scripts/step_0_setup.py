#!/usr/bin/env python3

"""
step_0a_setup.py - Script to automatically set up work environment

INSTRUCTIONS:

    Run this script in the terminal as `python3 step_0a_setup.py`
"""

import os
import config as c

os.makedirs(c.ORIGINAL_DATA, exist_ok=True)
os.makedirs(c.DATADIR, exist_ok=True)
os.makedirs(c.METADATA, exist_ok=True)
os.makedirs(c.ANALYSIS_DATA, exist_ok=True)
os.makedirs(c.IMPORTABLE_DATA, exist_ok=True)

print("Work environment setup concluded")
