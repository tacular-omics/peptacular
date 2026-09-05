# Manual agent evaluation

Use `agent_scenarios.json` with a local MCP client after installing this checkout.
These fixtures describe expected workflows. They do not claim a measured agent
success rate or require a provider API key in the automated test suite.

Register Peptacular using the command in `docs/mcp.rst`, then give the client one
scenario prompt at a time. All scenarios supply small synthetic sequences
inline. No workspace setup or input files are needed.

Record the client version, model, tool calls, arguments, result correctness,
retries, and whether it invented unsupported operations. Check each scenario's
listed assertions. A successful protocol connection alone does not establish
successful agent planning.

Confirm that agents pass actual annotations between calls, preserve any needed
source associations in their own context, and distinguish per-field diagnostics
from truncated calculations. After a truncated result, the agent should narrow
its request or adjust its bounded settings. It should not try to poll a job or
retrieve a stored result.

Restarting the server should require no cleanup or recovery. A subsequent call
supplies its inputs again and runs normally.
