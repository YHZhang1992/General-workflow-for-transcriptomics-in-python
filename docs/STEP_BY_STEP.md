# Step-by-step guide: Bulk transcriptomics

## Before execution

- Approve the analysis population, endpoint/contrast, units, missing-data rule, multiplicity strategy, and random seed.
- Validate identifiers before joins and confirm the independent experimental unit.
- Record input, configuration, reference, and code checksums.

## 1. Validate sample sheet and reference release

**Inputs:** verified artifacts from the previous step, with schema and checksum recorded.
**Action:** perform only this stage using frozen configuration and a clean run directory.
**Output:** a machine-readable table/object plus log entries sufficient to trace row/sample attrition.
**Acceptance:** command succeeds; required columns and counts are present; warnings and exclusions are reviewed.
## 2. Run sequencing/count QC when starting from FASTQ

**Inputs:** verified artifacts from the previous step, with schema and checksum recorded.
**Action:** perform only this stage using frozen configuration and a clean run directory.
**Output:** a machine-readable table/object plus log entries sufficient to trace row/sample attrition.
**Acceptance:** command succeeds; required columns and counts are present; warnings and exclusions are reviewed.
## 3. Normalize and model differential expression

**Inputs:** verified artifacts from the previous step, with schema and checksum recorded.
**Action:** perform only this stage using frozen configuration and a clean run directory.
**Output:** a machine-readable table/object plus log entries sufficient to trace row/sample attrition.
**Acceptance:** command succeeds; required columns and counts are present; warnings and exclusions are reviewed.
## 4. Run pathway scoring/enrichment

**Inputs:** verified artifacts from the previous step, with schema and checksum recorded.
**Action:** perform only this stage using frozen configuration and a clean run directory.
**Output:** a machine-readable table/object plus log entries sufficient to trace row/sample attrition.
**Acceptance:** command succeeds; required columns and counts are present; warnings and exclusions are reviewed.
## 5. Export tables, QC figures, and report

**Inputs:** verified artifacts from the previous step, with schema and checksum recorded.
**Action:** perform only this stage using frozen configuration and a clean run directory.
**Output:** a machine-readable table/object plus log entries sufficient to trace row/sample attrition.
**Acceptance:** command succeeds; required columns and counts are present; warnings and exclusions are reviewed.

## Final review

Perform structural, numerical, and scientific review. Archive the run manifest, environment/session information, logs, tables, figure source data, and reviewer disposition together.
