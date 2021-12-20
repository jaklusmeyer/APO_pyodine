# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 10:42:24 2021

Some functions to return information about the current git branch and
commit hash.

@author: pheeren
"""


import subprocess

def get_git_revision_hash() -> str:
    return subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()

def get_git_revision_short_hash() -> str:
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()

def get_git_branch_name() -> str:
    return subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'], universal_newlines=True).strip()