# Agent Setup Instructions

This is a hackathon project for **Track A** of the Monomer Bio / Elnora AI Science Hackathon. The goal is autonomous experimentation: growing *Vibrio Natriegens* on Monomer's robotic workcell, with AI agents designing, executing, and analyzing experiments iteratively.

Elnora is already set up. Two more MCPs need to be installed before proceeding.

---

## 1. Notion MCP

Gives the agent access to the hackathon's Notion workspace (Track A docs, reagent database, setup instructions).

**Install:**
```bash
claude mcp add notionApi -e NOTION_TOKEN=PASTE_TOKEN_HERE -- npx -y @notionhq/notion-mcp-server
```

> The Notion integration secret (`ntn_...`) will be provided by the user — wait for it before running this command.

Once installed, you can read pages like:
- Track A Experiment and Software Setup
- Track A Materials / Reagent Database
- Elnora Setup Instructions

---

## 2. Monomer Bio MCP

Gives the agent control over the robotic workcell: submitting Workflows, triggering Routines (liquid handling, measurements), and reading experiment data via the Culture Monitor.

**Reference:** https://github.com/monomerbio/monomer-bio-hackathon-march14-15/

The specific MCP server URL, auth token, and `claude mcp add` command are in that GitHub repo. Check the README or setup docs there. The Monomer team will also provide credentials on the day.

Once connected, the agent can:
- Submit a **Workflow** (Python) to run a series of **Routines** on the workcell
- Read absorbance / growth curve data from the **Culture Monitor**
- Iterate: design → execute → analyze → redesign

---

## Context: What Track A Is

- **Organism:** *Vibrio Natriegens*
- **Goal:** Iteratively find the optimal media composition for growth
- **Workcell instruments:** Opentrons Flex liquid handler, Liconic incubators (4°C + 37°C), Tecan Infinite platereader, KX-2 robotic arm
- **Plates:** 24-well deep well (8.5 mL/well) for reagent stocks; 2 mL/well cell culture plates
- **Loop:** Design workflow → Monomer staff reviews & loads → executes (~1.5 hrs) → analyze absorbance data → propose next design

## Reagents Available

| Component | Stock Concentration |
|-----------|-------------------|
| MOPS (pH 7) | 1 M |
| Ammonium Sulfate | 1 M |
| Iron(II) sulfate heptahydrate | 1 M |

(Full reagent database is in Notion)

---

## Elnora (already set up)

Elnora CLI is installed. Use it to generate and iterate on protocols:

```bash
elnora auth status   # verify auth
elnora --compact projects list
elnora --compact tasks create --project <ID> --title "Media optimization" --message "Generate a Vibrio Natriegens growth media optimization protocol"
```

Note: global flags like `--compact` must go **before** the subcommand.
