# VorPy Agent Operating Rules

## Routine work

- Routine source-code edits, documentation updates, and focused unit tests within
  this repository are allowed.
- Do not automatically commit or push changes.

## Approval required

Obtain explicit user approval before:

- deleting an existing file;
- modifying files outside this repository;
- installing, upgrading, or removing dependencies;
- changing system configuration;
- staging, committing, pushing, pulling, switching branches, merging, rebasing,
  or deleting a Git branch;
- running expensive molecular benchmarks, full-system `p53tet` runs, or large
  network exports; or
- modifying scientific reference values or numerical tolerances.

Never execute destructive Git operations, including `git reset --hard`,
`git clean -fd`, or a force-push, without explicit authorization.

## Protect in-progress work

- Preserve all existing tracked and untracked changes.
- Never overwrite, discard, or revert another agent's work.

## Scientific work

- Small synthetic geometry tests, EDTA validation, and bounded cambrin
  diagnostics are allowed.
- Scientific changes must document assumptions, units, numerical methods, and
  validation evidence.
- Coordinate changes to shared interfaces among FORGE, GAUSS, ATLAS, and VECTOR
  with the affected owners or agents before making incompatible changes.

## Codex Context and Token Efficiency

### Objective

Minimize unnecessary model requests, repeated context processing,
and excessive session growth while maintaining scientific correctness,
test coverage, and reproducibility.

Token efficiency must never justify weakening scientific validation,
silently changing numerical tolerances, or omitting necessary tests.

### 1. Bounded Development Tasks

- Work on one clearly defined milestone at a time.
- Do not expand the assigned task into unrelated refactoring.
- Do not perform a comprehensive repository audit unless explicitly requested.
- Prefer localized modifications over architectural rewrites.
- Stop when the assigned acceptance criteria are satisfied.
- Report additional opportunities rather than implementing them without approval.

### 2. Repository Exploration

- Read AGENTS.md once at the beginning of the session.
- Inspect only files relevant to the assigned milestone.
- Prefer targeted searches:

  `rg -n "symbol_name" vorpy/src/`

- Prefer bounded file reads:

  `sed -n '120,220p' path/to/file.py`

- Avoid repeatedly reading unchanged files.
- Do not recursively print large directories or entire source trees.
- Avoid loading generated artifacts, archives, or large molecular datasets
  into the conversation unless necessary.
- Use existing documentation and previous milestone summaries to avoid
  repeating completed investigations.

### 3. Tool Output Management

- Keep terminal output concise.
- Prefer focused pytest commands over running the full test suite repeatedly.
- Use `pytest -q` for routine validation.
- Capture verbose output in temporary files when appropriate.
- Inspect only relevant error messages and traceback sections.
- Avoid printing complete logs, large CSV files, or extensive profiling traces.
- Do not repeatedly display the same successful test output.

### 4. Context Growth Management

- Monitor session context growth when usage information is available.
- Prefer completing a milestone in a fresh session with a concise handoff
  rather than extending an already large conversation.
- When context becomes unusually large, stop at a safe checkpoint.
- Do not continue an extended investigation solely to exhaust the
  remaining context window.
- Avoid unnecessary summaries of information already present.
- Do not repeatedly reproduce large code blocks or documentation.
- Preserve important findings in concise project documentation rather
  than relying on an indefinitely growing conversation.

### 5. Agent Coordination

- Each agent owns its assigned development area.
- Do not spawn additional agents unless explicitly authorized.
- Do not duplicate another agent's investigation.
- Do not modify another agent's files without coordination.
- Read the relevant handoff documentation before beginning work.
- If blocked by another agent's unfinished work, report the dependency
  and stop rather than independently rebuilding that work.

### 6. Testing and Validation

- Run the smallest test set that adequately validates the modification.
- Expand testing only when failures or affected dependencies justify it.
- Do not repeatedly run expensive molecular benchmarks.
- Prefer bounded synthetic systems and EDTA regression cases.
- Preserve existing numerical tolerances unless changes are explicitly
  justified and approved.
- Never substitute passing tests for mathematical or scientific validity.
- Report tests that were not run and explain why.

### 7. Performance Investigations

- Profile before optimizing.
- Identify the dominant bottleneck before modifying implementation.
- Implement one measurable optimization per milestone.
- Compare before and after using identical inputs and settings.
- Preserve geometry, topology, ownership, and scientific status semantics.
- Stop if optimization requires substantial architectural changes
  outside the assigned milestone.

### 8. Session Completion and Handoff

At the end of each milestone, provide a concise report containing:

1. Completed objective.
2. Files modified.
3. Tests executed and results.
4. Scientific or implementation limitations.
5. Remaining dependencies.
6. Recommended next action.

Keep the final report under 300 words unless additional detail is requested.

When approaching context or usage limits, provide a checkpoint containing:

- Current task and completion status.
- Modified files.
- Verified results.
- Unresolved issues.
- Exact next step.

Do not begin another milestone without authorization.

### 9. Repository Safety

- Preserve concurrent agent changes.
- Do not discard uncommitted work or untracked research artifacts.
- Do not change branches, stage files, commit, or push without approval.
- Do not change dependencies or project configuration without approval.
- Do not sacrifice scientific correctness to reduce token consumption.


### Python Environment

- The verified VorPy development environment is `.venv` at the repository root.
- Always execute tests using `.venv/bin/python -m pytest`.
- Do not use bare `pytest` or `python -m pytest` unless the active interpreter
  has been explicitly verified.
- Do not assume Conda `base` contains compatible project dependencies.
- Do not install, upgrade, or change dependencies without approval.
- If the project interpreter is unavailable, report the problem and stop.