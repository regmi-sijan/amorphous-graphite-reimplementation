# Contributing

This is an independent reproduction study inspired by the cited paper. Changes
should make the workflow easier to reproduce and its scientific assumptions
clearer. Preserve attribution and distinguish published findings from newly
computed results.

## Reporting a problem

Include the command, configuration, Python/platform details, relevant dependency
versions, and the smallest useful log excerpt. State expected and actual behavior.
For scientific discrepancies, provide units, geometry provenance, convergence
settings and the cited comparison. Avoid uploading large model files or full
wavefunction/trajectory dumps to an issue; share an appropriate data reference.

## Proposing changes

Explain the scientific or usability problem, the change, and relevant validation.
Run `python -m unittest -v test_structure.py test_workflow.py` for affected analysis
code. For changes to engines, forces, stress or transport, include a suitably
small actual calculation or an analytic check that can reveal a wrong result.
Documentation-only edits do not require rerunning expensive simulations.

Do not replace the force model silently, prescribe layers in random initial
structures, or change formulas just to match a target curve. Document parameter
choices, definitions, method substitutions and failed convergence explicitly.
Keep generated runs and external models outside the source repository.

The maintainer has chosen to leave the original project code unlicensed for now.
No contributor licensing agreement or automatic rights grant is established by
this guide; clarify contribution terms with the maintainer before contributing.
