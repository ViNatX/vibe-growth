"""
Post seeding density experiment results as Notion tables.
Replaces previously posted bullet-list blocks with proper table blocks.
"""

import json
import requests

# Load token from .mcp.json
with open(".mcp.json") as f:
    mcp = json.load(f)
TOKEN = mcp["mcpServers"]["notionApi"]["env"]["NOTION_TOKEN"]

PAGE_ID = "7de9fbd45d864039891965fad9230937"
HEADERS = {
    "Authorization": f"Bearer {TOKEN}",
    "Notion-Version": "2022-06-28",
    "Content-Type": "application/json",
}


def rt(text):
    """Build a rich_text array from a plain string."""
    return [{"type": "text", "text": {"content": str(text)}}]


def table_row(*values):
    return {
        "object": "block",
        "type": "table_row",
        "table_row": {"cells": [rt(v) for v in values]},
    }


def delete_block(block_id):
    r = requests.delete(f"https://api.notion.com/v1/blocks/{block_id}", headers=HEADERS)
    if r.status_code not in (200, 404):
        print(f"  Warning: failed to delete {block_id}: {r.status_code}")


def append_blocks(block_id, children):
    r = requests.patch(
        f"https://api.notion.com/v1/blocks/{block_id}/children",
        headers=HEADERS,
        json={"children": children},
    )
    if not r.ok:
        print(f"  Error {r.status_code}: {r.text}")
    r.raise_for_status()
    return r.json()


# ── ViNatX rows ───────────────────────────────────────────────────────────────
VINATX_ROWS = [
    table_row("Wells", "Media (µL)", "Cell Stock (µL)", "Cell %", "OD600 Rep1", "OD600 Rep2", "OD600 Rep3"),
    table_row("A2–A4", 100, 100, "50%", 0.87, 0.79, 0.75),
    table_row("B2–B4", 140,  60, "30%", 0.46, 0.40, 0.39),
    table_row("C2–C4", 160,  40, "20%", 0.30, 0.28, 0.27),
    table_row("D2–D4", 170,  30, "15%", 0.22, 0.19, 0.18),
    table_row("E2–E4", 180,  20, "10%", 0.14, 0.15, 0.14),
    table_row("F2–F4", 190,  10,  "5%", 0.12, 0.12, 0.12),
]

# ── CellAI rows — 1 row per condition, reps as columns ────────────────────────
# Layout: cols 1–10 = conditions (50%→5%), rows A/B/C = reps
CELLAI_ROWS = [
    table_row("Cell %", "Media (µL)", "Cell Stock (µL)", "OD600 Rep A", "OD600 Rep B", "OD600 Rep C"),
    table_row("50%", 100, 100, 0.9491, 0.4867, 0.4278),
    table_row("45%", 110,  90, 0.4257, 0.4271, 0.4595),
    table_row("40%", 120,  80, 0.4674, 0.5071, 0.5098),
    table_row("35%", 130,  70, 0.5001, 0.5344, 0.5386),
    table_row("30%", 140,  60, 0.5278, 0.5586, 0.5492),
    table_row("25%", 150,  50, 0.5130, 0.5091, 0.5035),
    table_row("20%", 160,  40, 0.4743, 0.4833, 0.4783),
    table_row("15%", 170,  30, 0.4587, 0.4181, 0.1261),
    table_row("10%", 180,  20, 0.1111, 0.1327, 0.1253),
    table_row( "5%", 190,  10, 0.1197, 0.1193, 0.1207),
]

# ── Delete old CellAI table and its heading, then replace ────────────────────
CELLAI_HEADING_ID = "323f9baa-cde3-8192-941c-e535c6ef47a0"
CELLAI_TABLE_ID   = "323f9baa-cde3-8188-83fe-e70cfdcd762f"
VINATX_TABLE_ID   = "323f9baa-cde3-8163-a96e-c8a7d51b9600"  # anchor: insert after

print("Deleting old CellAI heading and table...")
delete_block(CELLAI_TABLE_ID)
delete_block(CELLAI_HEADING_ID)
print("  Deleted.")

print("Posting new CellAI heading + table after ViNatX table...")
result = append_blocks(PAGE_ID, [
    {
        "object": "block",
        "type": "heading_3",
        "heading_3": {"rich_text": rt("CellAI — first OD read at 2:23 PM PT")},
    },
    {
        "object": "block",
        "type": "table",
        "table": {
            "table_width": 6,
            "has_column_header": True,
            "has_row_header": False,
            "children": CELLAI_ROWS,
        },
    },
])

print(f"Created {len(result['results'])} blocks.")
print("Done!")
