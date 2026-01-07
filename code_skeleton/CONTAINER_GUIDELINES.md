# Container Development Guidelines

## Dockerfile Standards

### 1. No ENTRYPOINT
**Do NOT use `ENTRYPOINT` in your Dockerfiles.**
- **Reason**: Pipeline task specifications often need to override the command (e.g., to run `/bin/bash -c` for argument parsing or debugging). `ENTRYPOINT` forces its command to run and typically treats the pipeline's command as arguments, causing unexpected failures (like the `cp -c` error).
- **Alternative**: Use the pipeline task specification (e.g., in `task-spec.json`) to define the full command.

### 2. Script Locations
- Place executable scripts in `/usr/local/bin/` to ensure they are in the default `$PATH`.
- Example: `COPY run_script.sh /usr/local/bin/`
