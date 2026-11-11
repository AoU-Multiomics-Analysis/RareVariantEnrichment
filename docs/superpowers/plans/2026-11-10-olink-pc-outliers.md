# Olink PC outlier implementation plan

**Goal:** Run the existing outlier method at several PC counts in AoU Workbench and return one residual Z-score matrix per NPX type and PC count.

**Spec:** `docs/superpowers/specs/2026-11-10-olink-pc-outliers-design.md`

## Tasks

- [ ] Implement and test the staged R CLI. Compare prefix projection with direct QR regression. Cover PC=0, both NPX types, independent missingness, constants, duplicate keys, localization errors, and quoted paths.
- [ ] Implement the self-contained WDL 1.0 with typed files, task-local list creation, resource inputs, per-type preparation, per-PC adjustment, and summary collection. Validate with miniwdl and Cromwell's parser.
- [ ] Add the pinned micromamba image and Actions build/smoke/publish job. Keep image rebuilds restricted to bundled production code and build definitions. Run the complete synthetic workflow in Actions.
- [ ] Add an AST regression against workflow-scope file-writing functions and command-rendering tests that replace cloud File paths with localized files.
- [ ] Document current AoU registry, Google Batch, storage, and launch considerations. Add example inputs and Dockstore registration.
- [ ] Run required tests and independent review. Push the code and follow GitHub checks to completion. Report actual validation limits and image availability.

## Review focus

- Matrix orientation, string IDs, and correspondence between Z scores and outlier calls.
- Unavailable PC rank must stop the workflow instead of silently changing the grid.
- File values must survive until task localization; no URI hidden in an argument wrapper.
- Shared preparation must preserve the original analysis for every PC prefix.
- No controlled data, credentials, or computed participant results in the public repository or CI.
