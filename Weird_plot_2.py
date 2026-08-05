import plotly.graph_objects as go
import numpy as np

# -------------------------
# Input data
# -------------------------
cell_types = ["SynTI", "SynTII", "GlyT1"]

counts_13_5 = {
    "SynTI": 5235,
    "SynTII": 3184,
    "GlyT1": 9292
}

counts_16_5 = {
    "SynTI": 7818,
    "SynTII": 7016,
    "GlyT1": 6602
}

# -------------------------
# Layout settings
# -------------------------
x_left = 0.15
x_right = 0.85
bar_width = 0.08
gap = 0.04

colors = {
    "GlyT1": "rgba(0, 114, 178, 0.85)",
    "SynTI":  "rgba(213, 94, 0, 0.85)",
    "SynTII":  "rgba(204, 121, 167, 0.85)"
}

ribbon_colors = {
    "GlyT1": "rgba(0, 114, 178, 0.35)",
    "SynTI":  "rgba(213, 94, 0, 0.35)",
    "SynTII":  "rgba(204, 121, 167, 0.35)"
}

# -------------------------
# Normalize heights separately
# so 16.5 can be bigger overall
# -------------------------
left_total = sum(counts_13_5.values())
right_total = sum(counts_16_5.values())

left_scale = 0.7 / left_total
right_scale = 0.9 / right_total

# -------------------------
# Compute stacked positions
# -------------------------
left_positions = {}
right_positions = {}

y_top_left = 0.9
for ct in cell_types:
    h = counts_13_5[ct] * left_scale
    left_positions[ct] = (y_top_left - h, y_top_left)
    y_top_left -= h + gap

y_top_right = 0.95
for ct in cell_types:
    h = counts_16_5[ct] * right_scale
    right_positions[ct] = (y_top_right - h, y_top_right)
    y_top_right -= h + gap

# -------------------------
# Helper to make smooth ribbon
# -------------------------
def ribbon_path(x0, x1, y0_bottom, y0_top, y1_bottom, y1_top, curve=0.22):
    c0 = x0 + curve * (x1 - x0)
    c1 = x1 - curve * (x1 - x0)

    return (
        f"M {x0},{y0_top} "
        f"C {c0},{y0_top} {c1},{y1_top} {x1},{y1_top} "
        f"L {x1},{y1_bottom} "
        f"C {c1},{y1_bottom} {c0},{y0_bottom} {x0},{y0_bottom} Z"
    )

# -------------------------
# Build figure
# -------------------------
fig = go.Figure()

# Add ribbons first
for ct in cell_types:
    y0b, y0t = left_positions[ct]
    y1b, y1t = right_positions[ct]

    fig.add_shape(
        type="path",
        path=ribbon_path(
            x_left + bar_width / 2,
            x_right - bar_width / 2,
            y0b, y0t, y1b, y1t
        ),
        fillcolor=ribbon_colors[ct],
        line=dict(width=0)
    )

# Add left and right bars
for ct in cell_types:
    # left block
    y0, y1 = left_positions[ct]
    fig.add_shape(
        type="rect",
        x0=x_left - bar_width / 2,
        x1=x_left + bar_width / 2,
        y0=y0,
        y1=y1,
        fillcolor=colors[ct],
        line=dict(width=0)
    )

    # right block
    y0, y1 = right_positions[ct]
    fig.add_shape(
        type="rect",
        x0=x_right - bar_width / 2,
        x1=x_right + bar_width / 2,
        y0=y0,
        y1=y1,
        fillcolor=colors[ct],
        line=dict(width=0)
    )

# Labels
for ct in cell_types:
    yl = sum(left_positions[ct]) / 2
    yr = sum(right_positions[ct]) / 2

    fig.add_annotation(
        x=x_left - 0.09, y=yl,
        text=f"{ct}<br>n={counts_13_5[ct]}",
        showarrow=False,
        xanchor="right",
        font=dict(size=24)
    )

    fig.add_annotation(
        x=x_right + 0.09, y=yr,
        text=f"{ct}<br>n={counts_16_5[ct]}",
        showarrow=False,
        xanchor="left",
        font=dict(size=24)
    )

# Column titles
fig.add_annotation(x=x_left, y=1.03, text="<b>E13.5</b>", showarrow=False, font=dict(size=16))
fig.add_annotation(x=x_right, y=1.03, text="<b>E16.5</b>", showarrow=False, font=dict(size=16))

# Layout cleanup
fig.update_xaxes(visible=False, range=[0, 1])
fig.update_yaxes(visible=False, range=[0, 1.08])
fig.update_layout(
    title="Cell Proliferation, Sankey-Style",
    width=1000,
    height=700,
    plot_bgcolor="white",
    paper_bgcolor="white",
    margin=dict(l=40, r=40, t=70, b=30)
)

fig.show()
