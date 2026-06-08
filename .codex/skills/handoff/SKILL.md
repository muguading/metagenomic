---
name: handoff
description: Compact the current conversation into a handoff document for another agent to pick up. Use when the user explicitly asks for this skill or the workflow it describes. Upstream marks this skill as in progress.
---

> Codex conversion note: Upstream marks this skill as in progress; expect to refine it before relying on it heavily.


Write a handoff document summarising the current conversation so a fresh agent can continue the work. Save it to a path produced by `mktemp -t handoff-XXXXXX.md` (read the file before you write to it).

Suggest the skills to be used, if any, by the next session.

Do not duplicate content already captured in other artifacts (PRDs, plans, ADRs, issues, commits, diffs). Reference them by path or URL instead.

If the user passed arguments, treat them as a description of what the next session will focus on and tailor the doc accordingly.
