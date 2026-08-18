# ChatGPT Collaboration Workflow

This repository is designed to remain understandable across ChatGPT
conversations.

## Source of Truth

GitHub is the source of truth.

Conversation history and temporary ChatGPT storage are not persistent project
state and must not be required to continue development.

## Normal Bootstrap

A new conversation should read, in order:

1. `README.md`
2. `AGENTS.md`
3. `docs/PROJECT_STATE.md`
4. `docs/STAGES.md`

Then inspect the specific code relevant to the requested task.

## When ChatGPT Cannot Push Directly

Use the one-click update workflow:

1. ChatGPT prepares an update package containing a Windows launcher and the
   replacement/patch content.
2. The user double-clicks the launcher.
3. The launcher:
   - checks for Git;
   - clones the repository if needed;
   - updates `main`;
   - creates or recreates a work branch;
   - applies the prepared change;
   - commits;
   - pushes the branch.
4. ChatGPT reads that pushed branch from GitHub before continuing.

The user should not need to manually type a series of Git commands for routine
updates.

## Persistence Rule

If a decision, failure, experiment, workaround, or limitation would matter in a
fresh conversation, record it in `docs/PROJECT_STATE.md` before considering the
work complete.
