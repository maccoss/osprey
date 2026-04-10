# Release Notes

This directory contains per-version release notes for Osprey.

## Versioning Scheme

Osprey uses a `YY.feature.patch` versioning convention:

- **YY**: Two-digit year (e.g., `26` for 2026)
- **feature**: Incremented for each release containing new features
- **patch**: Incremented for bug-fix-only releases within the same feature version

Examples: `26.1.0` (first feature release of 2026), `26.1.1` (patch), `26.2.0` (second feature release).

The workspace version in `Cargo.toml` is updated only at release time, not during development.

## File Format

Each release gets one file: `RELEASE_NOTES_v{version}.md`. During development, the unreleased draft lives in `RELEASE_NOTES_next.md` and gets renamed at release time.

```text
release-notes/
  README.md                      # this file
  RELEASE_NOTES_next.md          # working draft for the next release
  RELEASE_NOTES_v26.1.0.md
  RELEASE_NOTES_v26.1.1.md
  RELEASE_NOTES_v26.2.0.md
```

## Writing Release Notes

### During Development

Maintain a draft release notes file for the next planned version (e.g., `RELEASE_NOTES_v26.2.0.md`). Append entries as features and fixes land on the development branch. This file is a working draft until the release is finalized.

### Content Structure

Each release notes file should use this structure:

```markdown
# Osprey v{version} Release Notes

One-sentence summary of the release.

## New Features

- Feature descriptions grouped by area (e.g., FDR Control, Reconciliation)
- Focus on what changed from the user's perspective, not implementation details

## Bug Fixes

- Description of the bug and its impact
- What was fixed

## Performance

- Performance improvements with context (e.g., "Reduced memory from 35 GB to 5 GB for 240-file experiments")

## Breaking Changes

- Any changes that require user action (config format changes, removed options, etc.)
- Omit this section if there are no breaking changes
```

Sections can be omitted if empty. For major releases with many changes, subsections within each category are fine. For patch releases, a flat list is sufficient.

### Style

- Write in past tense ("Added", "Fixed", "Removed")
- Lead with user impact, not implementation details
- Include specific numbers where relevant (e.g., memory reduction, precursor counts)
- Reference config options by their CLI flag or YAML key name

## Release Process

1. Finalize `RELEASE_NOTES_next.md` on the development branch
2. Rename it: `git mv release-notes/RELEASE_NOTES_next.md release-notes/RELEASE_NOTES_v{version}.md`
3. Update the title heading inside the file to match the version
4. Create a fresh empty `RELEASE_NOTES_next.md` for the following release
5. Update `version` in the workspace `Cargo.toml` to match the release
6. Commit the version bump and renames
7. Merge to `main`
8. Tag: `git tag v{version}`
9. Push: `git push origin main --tags`
10. CI builds and publishes the release artifacts
