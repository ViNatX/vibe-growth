# vibe-growth

AI agent platform for optimizing *Vibrio Natriegens* cell culture media, built for Track A of the Monomer Bio / Elnora AI Science Hackathon. The agent iteratively designs, executes, and analyzes experiments on a robotic workcell.

**Full agent context (Notion MCP, Monomer Bio MCP, workcell details):** see `AGENT_SETUP.md`

---

## Environment

This project uses [pixi](https://pixi.sh) for dependency management.

```bash
pixi install        # install all dependencies
pixi run <command>  # run a command in the pixi environment
```

**Before running any command**, always ensure the environment is fully set up:

1. `git pull` — pull the latest changes from teammates first
2. `git submodule update --init` — pull in the Elnora CLI submodule
3. `pixi install` — install all dependencies

Do not skip these steps or assume they have already been run. Teammates should pull often to stay in sync.

---

## Elnora CLI Setup (run once per machine)

The Elnora CLI plugin is checked in as a submodule at `./elnora-cli` and the marketplace manifest is at `.claude-plugin/marketplace.json`.

```bash
git submodule update --init

claude plugin marketplace add /absolute/path/to/vibe-growth

pixi run elnora auth login
```

---

## Team Collaboration via Notion

This is a shared project. At the start of every session, Claude should:

1. **Ask who you are** — so activity can be attributed correctly (e.g. "What's your name?").
2. **Check the vibe-growth Notion project page** — to see what has already been done and what teammates are currently working on, before starting any new work.
3. **Update the vibe-growth Notion project page after completing any analysis** — log what was done, by whom, and any key findings or next steps.

This keeps the team in sync across Elnora, Monomer, and local workflows.

> Example: Raul is setting up MCP connections and environments. After finishing, log it on the Notion page so teammates know the setup status.

---

## Notion MCP Setup (run once per machine)

`.mcp.json` is gitignored (contains secrets). The Notion token is also stored as a GitHub secret (`NOTION_TOKEN`) for CI use, but since GitHub secrets are write-only you'll need to get the value directly from a teammate. Then run:

```bash
claude mcp add notionApi -e NOTION_TOKEN=<your_token> --scope project -- npx -y @notionhq/notion-mcp-server
```
