# ADR 0008: Archive agent sessions as sanitized replays and publish only the HTML

- Status: accepted
- Date: 2026-10-15

## Context

The user's `/publish-agent-trace` workflow archives the AI coding sessions that built the project.

## Decision

- `agent-traces/<slug>/` holds:
  - a sanitized `session.jsonl`: email redacted, home path shown as `~`, internal API-request and cache records dropped;
  - `replay.html`, generated with the ASA fork of claude-replay;
  - a README.
- Only the HTML replay is copied to `docs/agent-traces/` for Pages. Raw JSONL is never published there.
- Failed attempts stay in the trace.

## Alternatives

- Publishing the raw JSONL. Too much risk of leaking sensitive content.
- Not publishing at all. The user asked for publication.

## Consequences

Replays are public. Every new trace needs a sensitivity review before it is committed. Known residuals in the current trace: a directory listing of `~/git/asa` and the name of an API gateway in an error message.
