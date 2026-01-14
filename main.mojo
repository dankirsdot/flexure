"""
Flexure Amplifier Analysis - Main Entry Point

This is the main entry point for the flexure amplifier library.
It runs example analyses demonstrating the library capabilities.

The library implements the unified two-port dynamic stiffness model from:
"A general two-port dynamic stiffness model and static/dynamic comparison
for three bridge-type flexure displacement amplifiers"
by Mingxiang Ling (2019, Mechanical Systems and Signal Processing)

Usage:
    mojo run main.mojo
    
Outputs:
    - Console: Single-point analysis results
    - fig4_rhombic.png: Rhombic amplifier sweep (Figure 4)
    - fig5_parallel.png: Parallel amplifier sweep (Figure 5)
    - fig6_aligned.png: Aligned amplifier sweep (Figure 6)
"""

from examples import run_all_examples


fn main():
    """Main entry point - run all examples."""
    run_all_examples()