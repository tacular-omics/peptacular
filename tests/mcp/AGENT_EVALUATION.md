# Manual agent evaluation

Use `agent_scenarios.json` with a local MCP client after installing this checkout.
These fixtures describe expected workflows. They do not claim a measured agent
success rate or require a provider API key in the automated test suite.

Run scenarios in a temporary workspace. For the FASTA scenario create this input:

```text
>duplicate first record
AKPEPTIDERAAK
>duplicate second record
MPEPTIDERAAK
```

Register Peptacular using the command in `docs/mcp.rst`, then give the client one
scenario prompt at a time. Record the client version, model, tool calls, arguments,
result correctness, retries, and whether it invented unsupported operations.
Check each scenario's listed assertions. A successful protocol connection alone
does not establish successful agent planning.

For cancellation, submit a deliberately large bounded job, cancel it, and verify
a terminal state with `get_job`. Check that the server still handles a small
inspection request afterward. For restart recovery, restart the client and use
`list_workspace` to locate retained completed results. Earlier live jobs should
report unavailable state, with no claim that they continued after disconnect.

Never use production research inputs for these setup checks. All examples above
are small synthetic sequences.
