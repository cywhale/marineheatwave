# `tools/p2/probe.py`: fixed-request GetMap probe (Phase 2, Stage 0)

A small, dependency-free (standard library only) tool for timing a fixed list of
WMS/HTTP requests against one or two endpoints, with an optional paired A/B
design. It is a measurement tool only. It changes no server setting or client,
and it never purges or seeds a cache. The requests themselves are ordinary
requests, though: one that misses a cache (nginx, GeoWebCache) may leave a new
entry there, and that side effect must be part of the batch's approval.

Stage 0 covers only the tool and its local tests (`tests/test_p2_probe.py`,
which use a local fake HTTP server). Running it against a real GeoServer is a
separate, separately approved measurement batch.

## Commands

From the `API/` directory, with the project's test or tools venv:

```bash
python -m tools.p2.probe plan path/to/manifest.json
python -m tools.p2.probe run path/to/manifest.json --out path/to/new-output-dir
```

- `plan` validates the manifest and prints every request it would send, in
  order, with the full URL. It makes no network connection.
- `run` executes the batch. Options:
  - `--max-requests N` sets the batch cap and overrides the manifest's
    `max_requests`.
  - `--save-bodies` also writes every response body to `OUT/bodies/`.

Exit codes:

| Code | Meaning |
|---|---|
| 0 | complete |
| 2 | stopped because a request failed |
| 3 | invalid manifest, or the planned requests exceed the cap; nothing was sent |
| 4 | the output directory is not empty |

## What the tool does and does not do

- **Batch cap.** The cap is always explicit: `max_requests` in the manifest or
  `--max-requests`. If neither is given, the batch does not start. There is no
  built-in limit. If the planned sequence is larger than the cap, the batch does
  not start; it is never cut short to fit.
- **URLs.** Each URL is the target's `base_url`, plus the case's raw `query`,
  plus the arm's raw `extra_query` (when present), joined with `?`/`&`, exactly
  as written. Nothing is re-encoded or added: no cache-busting parameter.
- **No proxies.** Proxy environment variables are ignored; connections go
  directly to the target.
- **One attempt.** Each request is sent once. The first failure stops the batch,
  and the failure is written to the record. Failures are:
  - a connection error, refused connection or TLS error;
  - a timeout (`timeout_s` covers the whole request, including the body);
  - a status outside `expect.status`;
  - a wrong `Content-Type`;
  - a body that is not a PNG (when `expect.png` is on);
  - a body shorter than its `Content-Length`.
- **No cache operations.** The tool never purges or seeds a cache. A request
  that misses a cache may add an entry there; that is a property of the target,
  and must be part of the batch's approval.
- **Secrets.**
  - Manifests with credentials are rejected:
    - a user name or password in a URL;
    - a credential-like query parameter in `base_url`, a case `query` or an
      arm's `extra_query`: names such as `key`, `auth`, `authkey` (GeoServer's
      authkey module), `apikey`/`api_key`, `sig`/`signature`, `password`,
      `session`/`sessionid`/`sid`, or any name containing `token`, `secret`,
      `password`, `signature` or `session`;
    - `Authorization`, `Proxy-Authorization`, `Cookie`, `X-Api-Key` or
      `X-Auth-Token` request headers.
  - `Set-Cookie` and other cookie or authorization response headers are written
    as `[redacted]`.
  - URL-valued response headers (`Location`, `Content-Location`, `Refresh`,
    `Link`) are kept, but with URL user info and the values of credential-like
    query parameters replaced by `[redacted]`. The same redaction is applied to
    the recorded request URL and to the URLs `plan` prints.

## Timings

All times are in seconds, measured with `time.perf_counter()` from the start of
each request:

| Field | Meaning |
|---|---|
| `connect` | name resolution (DNS) + TCP connect + TLS handshake for https, finished |
| `first_byte` | the socket became readable after the request was sent (first response byte) |
| `headers` | status line and headers parsed |
| `total` | complete body received |

Each request opens a new connection (no keep-alive), so the times are those
of independent requests, including name resolution and connection setup. A
browser that reuses connections sees less of that overhead; interpret
browser-perceived latency separately.

## Design

- **Single arm** (`arms` has only `A`): for each case, `warmup` requests, then
  `repeats` timed requests.
- **Paired** (`A` and `B`): for each case, the sequence is:
  1. `warmup` rounds of A then B, not timed;
  2. `a_reference` A-only requests;
  3. `pairs` pairs in balanced order: AB, BA, BA, AB, AB, BA, ...;
  4. `a_reference` A-only requests again.
- **Arms.** Each arm is a target plus an optional extra query string. So A and B
  can be two endpoints with the same request, or one endpoint with two variants
  (for example two `STYLES`).
- **A drift.** This is the relative change of A over the batch, in percent:
  - when `a_reference` > 0: the A-only block after the pairs compared with the
    block before;
  - otherwise: the second half of A's paired observations compared with the
    first half.
- **Sample size.** How many warm-ups, pairs and references to use, and when a
  result is too noisy to read, is decided per batch, not by the tool.

## Manifest format (`schema: mhw-p2-probe/1`)

See [`example_manifest.json`](example_manifest.json). The manifest is validated completely before
anything is sent, so `plan` reports every input error:
- unknown fields;
- wrong types;
- invalid URLs or ports;
- a non-finite or non-positive timeout;
- a malformed `expect`;
- a case id that is not a safe file name. It uses a `localhost`
placeholder; replace it with a real target only in a manifest kept outside the
repository.

| Field | Required | Meaning |
|---|---|---|
| `schema` | yes | `"mhw-p2-probe/1"` |
| `batch` | yes | a label for this batch |
| `max_requests` | yes, unless `--max-requests` is given | approved request cap for this batch |
| `timeout_s` | yes | per-request limit, whole request including the body |
| `tls_verify` | no (default `true`) | `false` allows a self-signed certificate on a test target; the setting is recorded |
| `request_headers` | no | extra request headers; credentials are rejected |
| `targets` | yes | name → `{"base_url": "http(s)://host[:port]/path"}` |
| `arms` | yes | `{"A": {...}}` or `{"A": {...}, "B": {...}}`; each is `{"target": name, "extra_query": "K=V&..."}` |
| `cases` | yes | list of `{"id", "query", "expect"}`; `id` matches `[A-Za-z0-9][A-Za-z0-9._-]{0,63}`; `expect` defaults to `{"status": [200], "content_type": "image/png", "png": true}` (`content_type` may be `null`) |
| `design` | yes | `warmup`, then `pairs` (paired) or `repeats` (single arm), and optionally `a_reference` (paired only) |

## Output (`--out DIR`)

- **`requests.jsonl`:** one line per request sent, in order. Fields:
  - identity: `seq`, `case`, `arm`, `target`, `phase` (`warmup`, `a_ref_before`,
    `pair`, `a_ref_after` or `repeat`), `pair`, `position`, `started_at` (UTC),
    `url`;
  - timing: `t`;
  - response: `status`, `reason`, `headers` (all response headers in order,
    secrets redacted), `body_bytes`, `sha256`, `png_size`;
  - result: `ok`, `failure`.
- **`summary.json`:** the batch result. Fields:
  - `status` (`complete` or `stopped`) and `stop` (the failed request and why);
  - `planned`, `sent`, `ok`, `failed`, `not_sent`;
  - the manifest's SHA-256;
  - per case:
    - per arm: median, min and max of `total`, `headers` and `first_byte`; the
      distinct SHA-256 values and body sizes; the cache-related headers seen;
    - for paired batches: every pair with B − A, B/A and both SHA-256 values
      (`same_bytes`); the paired difference statistics; `a_drift`; and the byte
      comparison counted pair by pair (`pairs_compared`, `pairs_same_bytes`,
      `all_pairs_same_bytes`). Digests are compared within each pair, never as
      sets: A returning X, Y while B returns Y, X counts as two mismatches.
- **`bodies/`:** written only with `--save-bodies`.

Real manifests (with internal endpoints) and run outputs must not be committed.
Keep them outside the repository.
