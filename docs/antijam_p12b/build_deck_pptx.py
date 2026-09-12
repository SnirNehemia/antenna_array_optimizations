"""Build the Mode C anti-jam briefing deck (.pptx) from the content below.

The authored source of truth for the briefing is ``modec_deck.html`` in this
folder; this script produces the PowerPoint rendering of the same story so the
deck can be handed to people who want to edit it. Content is duplicated here on
purpose -- the HTML is a 22-slide reading deck, this is an 18-slide speaking
deck -- so when a number changes, change it in both.

    python docs/antijam_p12b/build_deck_pptx.py

Writes ``results/antijam/p12b_modec_campaign/ModeC_AntiJam_Briefing.pptx``.
Chart images live in ``docs/antijam_p12b/deck_assets/``.

Not part of the optimization pipeline: the ``src/`` tree stays frozen at
Milestone 1, so this lives with the document it builds. Requires python-pptx.

Part of: Antenna Array Pattern Optimization Tool -- anti-jam milestone.
"""

import os

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

# ── design system, lifted from the HTML deck's CSS variables ──────────
INK = RGBColor(0x10, 0x15, 0x1C)
INK2 = RGBColor(0x3D, 0x47, 0x56)
INK3 = RGBColor(0x6B, 0x76, 0x88)
ACCENT = RGBColor(0x1C, 0x5C, 0xAB)
PANEL = RGBColor(0xF2, 0xF5, 0xF8)
PASS = RGBColor(0x00, 0x83, 0x00)
FAIL = RGBColor(0xC0, 0x39, 0x2F)
DARK_BG = RGBColor(0x10, 0x15, 0x1C)
DARK_PANEL = RGBColor(0x1B, 0x25, 0x34)
DARK_INK = RGBColor(0xEE, 0xF2, 0xF7)
DARK_INK2 = RGBColor(0xA9, 0xB6, 0xC6)
DARK_ACCENT = RGBColor(0x4E, 0x93, 0xE0)

SANS, SERIF, MONO = "Calibri", "Cambria", "Consolas"

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
ASSETS = os.path.join(HERE, "deck_assets")
OUT = os.path.join(REPO, "results", "antijam", "p12b_modec_campaign",
                   "ModeC_AntiJam_Briefing.pptx")


# ── primitives ───────────────────────────────────────────────────────
def textbox(slide, left, top, width, height):
    box = slide.shapes.add_textbox(Inches(left), Inches(top),
                                   Inches(width), Inches(height))
    tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    return tf


def para(tf, first=False, space_after=6, space_before=0, align=None):
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    p.space_after = Pt(space_after)
    p.space_before = Pt(space_before)
    if align is not None:
        p.alignment = align
    return p


def run(p, text, size, font=SANS, bold=False, color=INK2):
    r = p.add_run()
    r.text = text
    r.font.size = Pt(size)
    r.font.name = font
    r.font.bold = bold
    r.font.color.rgb = color
    return r


def rect(slide, left, top, width, height, fill):
    from pptx.enum.shapes import MSO_SHAPE
    sh = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(left), Inches(top),
                                Inches(width), Inches(height))
    sh.fill.solid()
    sh.fill.fore_color.rgb = fill
    sh.line.fill.background()
    sh.shadow.inherit = False
    return sh


def blank(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])


def chrome(slide, number, kicker, title):
    """Page number, eyebrow and title -- every content slide carries these."""
    tf = textbox(slide, 12.08, 0.5, 0.6, 0.3)
    run(para(tf, first=True), "%02d" % number, 10, MONO, False, INK3)
    tf = textbox(slide, 0.72, 0.52, 9.0, 0.3)
    run(para(tf, first=True), kicker.upper(), 11, SANS, True, ACCENT)
    tf = textbox(slide, 0.72, 0.92, 11.89, 1.05)
    run(para(tf, first=True), title, 34, SERIF, True, INK)


def bullets(slide, left, top, width, height, items, size=15.5, gap=10):
    """items: (lead, rest) -- lead is bold and dark, rest is body colour."""
    tf = textbox(slide, left, top, width, height)
    for i, (lead, rest) in enumerate(items):
        p = para(tf, first=(i == 0), space_after=gap)
        if lead:
            run(p, lead + " ", size, SANS, True, INK)
        if rest:
            run(p, rest, size, SANS, False, INK2)
    return tf


def prose(slide, left, top, width, height, paragraphs, size=15, gap=10):
    tf = textbox(slide, left, top, width, height)
    for i, chunks in enumerate(paragraphs):
        p = para(tf, first=(i == 0), space_after=gap)
        if isinstance(chunks, str):
            chunks = [(chunks, False)]
        for text, bold in chunks:
            run(p, text, size, SANS, bold, INK if bold else INK2)
    return tf


def footnote(slide, text, top=6.35):
    tf = textbox(slide, 0.72, top, 11.4, 0.5)
    run(para(tf, first=True), text, 12.5, SANS, False, INK3)


def statbox(slide, left, top, width, height, value, label, value_color=INK,
            fill=PANEL, label_color=INK3, value_size=33):
    rect(slide, left, top, width, height, fill)
    tf = textbox(slide, left + 0.28, top + 0.26, width - 0.56, height - 0.5)
    run(para(tf, first=True, space_after=2), value, value_size, SERIF, True, value_color)
    run(para(tf, space_after=0), label, 11.5, SANS, False, label_color)


def table(slide, left, top, width, rows, col_widths, header_h=0.34, row_h=0.40,
          size=14, header_size=10.5, left_cols=()):
    """rows[0] is the header. A cell is a str, or (text, colour, bold)."""
    n_r, n_c = len(rows), len(rows[0])
    height = header_h + row_h * (n_r - 1)
    gf = slide.shapes.add_table(n_r, n_c, Inches(left), Inches(top),
                                Inches(width), Inches(height))
    tbl = gf.table
    tbl.first_row = False
    tbl.horz_banding = True
    for ci, w in enumerate(col_widths):
        tbl.columns[ci].width = Inches(w)
    for ri, row in enumerate(rows):
        tbl.rows[ri].height = Inches(header_h if ri == 0 else row_h)
        for ci, cell in enumerate(row):
            text, color, bold = (cell, None, None) if isinstance(cell, str) else cell
            c = tbl.cell(ri, ci)
            c.fill.solid()
            c.fill.fore_color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
            c.margin_left = c.margin_right = Inches(0.08)
            c.margin_top = c.margin_bottom = Inches(0.02)
            p = c.text_frame.paragraphs[0]
            p.alignment = (PP_ALIGN.LEFT if ci == 0 or ci in left_cols
                           else PP_ALIGN.RIGHT)
            if ri == 0:
                run(p, text, header_size, SANS, True, INK3)
            else:
                run(p, text, size, SANS,
                    bold if bold is not None else (ci > 0),
                    color if color is not None else INK)
    return gf


def picture(slide, name, left, top, width, height):
    slide.shapes.add_picture(os.path.join(ASSETS, name), Inches(left),
                             Inches(top), Inches(width), Inches(height))


def caption(slide, left, top, width, height, text):
    tf = textbox(slide, left, top, width, height)
    run(para(tf, first=True), text, 11.5, SANS, False, INK3)


# ── slides ───────────────────────────────────────────────────────────
def slide_title(prs):
    s = blank(prs)
    rect(s, -0.05, -0.05, 13.43, 7.6, DARK_BG)
    tf = textbox(s, 0.72, 1.62, 10.4, 0.35)
    run(para(tf, first=True), "ANTI-JAM MILESTONE  ·  P12b + PHASE O/O2  ·  2026-09-12",
        12.5, SANS, True, DARK_ACCENT)
    tf = textbox(s, 0.72, 2.19, 11.2, 2.3)
    for i, line in enumerate(["Mode C anti-jam:", "what works, what doesn't,",
                              "and what ships"]):
        run(para(tf, first=(i == 0), space_after=0), line, 42, SERIF, True, DARK_INK)
    tf = textbox(s, 0.72, 4.42, 10.6, 0.95)
    run(para(tf, first=True),
        "Three campaigns on the covariance-domain anti-jam stack: the drifting-jammer "
        "failure diagnosed and fixed, the on/off predictor found never to be running and "
        "repaired, a recommendation reversed by measurement — and a final amplitude sweep "
        "that names what is still on the table.",
        15.5, SANS, False, DARK_INK2)
    stats = [("6,246", "closed-loop runs"), ("7 × 5", "arrays × target positions"),
             ("0", "failed case-seeds"), ("59 / 59", "gates passing")]
    for i, (v, l) in enumerate(stats):
        statbox(s, 0.72 + i * 2.86, 5.72, 2.6, 1.15, v, l,
                DARK_INK, DARK_PANEL, DARK_INK2)


def slide_conclusions(prs):
    s = blank(prs)
    chrome(s, 2, "The bottom line", "Three conclusions")
    stats = [("2 → 11", "of 15 drifting-jammer cells passing, after the fix", INK),
             ("73.5 → 80.6", "on/off score once the predictor was made to fire", INK),
             ("8.1 dB", "of achievable gain still unclaimed at strong signal", FAIL)]
    for i, (v, l, c) in enumerate(stats):
        statbox(s, 0.72 + i * 4.03, 2.15, 3.75, 1.35, v, l, c, value_size=29)
    bullets(s, 0.72, 3.85, 11.4, 2.9, [
        ("Both algorithm failures were arithmetic, not tuning.",
         "The drift gap is angle lag of exactly one covariance horizon. The on/off predictor "
         "was never firing at all — two independent blockers. Both are fixed and both fixes "
         "ship."),
        ("Measurement reversed one of our own recommendations.",
         "Phase O advised against enabling the on/off repair on the strength of a regression "
         "that turned out to be a run-length artifact. Re-run at 8 cycles it improves or ties "
         "on every array."),
        ("What limits delivery is steering accuracy, not nulling.",
         "Calibration error past ~1° RMS collapses the stack, and at strong signal up to "
         "8.1 dB of achievable gain goes unclaimed. That is where the next work belongs."),
    ])


def slide_scope(prs):
    s = blank(prs)
    chrome(s, 3, "Scope", "What Mode C is, and what we are testing")
    prose(s, 0.72, 2.1, 5.7, 3.6, [
        "Mode C is the observation regime where the receiver provides per-element complex "
        "snapshots. The algorithm sees the raw array output and nothing else.",
        "It is the richer of the two regimes in this milestone — Mode S provides only a "
        "scalar SINR reading — and it is the one a real digital-beamforming receiver offers.",
        [("Scope is deliberately narrow: ", False), ("one jammer", True),
         (", unknown direction, either drifting slowly or switching on and off, in 2-D "
          "(θ, φ), against real CST element patterns.", False)],
    ])
    rect(s, 6.9, 2.1, 5.7, 1.15, PANEL)
    tf = textbox(s, 7.15, 2.42, 5.2, 0.6)
    run(para(tf, first=True),
        "SINR(w) = σs²|wᴴe_s|² / (σj²|wᴴe_j|² + σn²‖w‖²)", 14, MONO, False, INK)
    prose(s, 6.9, 3.5, 5.7, 2.2, [
        "Each element is multiplied by a complex weight and summed. The weights decide how "
        "much of the wanted signal survives and how much of the jammer does. Everything in "
        "this deck is a statement about how well a given rule chooses them.",
    ])
    footnote(s, "Algorithms may consume only the snapshots. The true jammer angle is never "
                "passed to an algorithm — only to the oracle, which exists to bound them.")


def slide_algorithms(prs):
    s = blank(prs)
    chrome(s, 4, "The algorithms", "Four algorithms, chosen to bracket the problem")
    gf = table(s, 0.72, 2.15, 11.8, [
        ["", "knows the jammer angle?", "what it is for"],
        [("oracle", INK, False), "yes, exactly", "the ceiling — not deliverable"],
        [("lcmv", INK, False), "no", "the reactive baseline"],
        [("predict", INK, False), "no — estimates it", "the candidate"],
        [("predict + repair", INK, False), "no — estimates it", "what now ships"],
    ], [2.6, 3.4, 5.8], row_h=0.44, left_cols=(1, 2))
    for ri in range(1, 5):
        for ci in (1, 2):
            p = gf.table.cell(ri, ci).text_frame.paragraphs[0]
            p.alignment = PP_ALIGN.LEFT
            for r in p.runs:
                r.font.bold = False
                r.font.color.rgb = INK2
    prose(s, 0.72, 4.6, 11.4, 2.0, [
        "They are not four attempts at the same thing. The oracle sets the ceiling and cannot "
        "be fielded, because it is handed the truth the others must infer. lcmv is the honest "
        "reactive baseline: it nulls whatever is loud, with no idea where the jammer is. "
        "predict adds a direction estimate and a model of where the jammer is going. The "
        "distance between the last two and the oracle is the whole subject of this deck.",
    ])


def slide_how(prs):
    s = blank(prs)
    chrome(s, 5, "How it works", "The reactive nuller needs no idea where the jammer is")
    prose(s, 0.72, 2.1, 11.4, 0.6, [
        "Minimise total output power, subject to holding unit gain on the target direction:",
    ])
    rect(s, 0.72, 2.72, 11.8, 0.95, PANEL)
    tf = textbox(s, 1.0, 3.03, 11.2, 0.5)
    run(para(tf, first=True), "w = R⁻¹e(θs) / (e(θs)ᴴR⁻¹e(θs))", 18, MONO, False, INK)
    bullets(s, 0.72, 4.0, 11.4, 2.3, [
        ("Anything loud that is not the target gets nulled automatically.",
         "The jammer's direction never appears in the formula — it appears in R, the measured "
         "covariance, which is all the receiver needs."),
        ("R is estimated with exponential forgetting,",
         "λ = 0.90, so the estimate looks back about 10 snapshots. That horizon is the source "
         "of the drift failure on slide 8, and the reason a null survives a jammer switching "
         "off on slide 14."),
    ])
    footnote(s, "This is MVDR/LCMV in its MPDR form — the wanted signal is inside R. "
                "Convenient in practice, and the cause of the last finding in this deck.")


def slide_metric(prs):
    s = blank(prs)
    chrome(s, 6, "Method", "One number per cell, normalised against the achievable")
    rect(s, 0.72, 2.15, 11.8, 0.95, PANEL)
    tf = textbox(s, 1.0, 2.46, 11.2, 0.5)
    run(para(tf, first=True),
        "track_score = 100 · mean( (SINR_oracle − SINR) ≤ 3 dB )", 17, MONO, False, INK)
    bullets(s, 0.72, 3.4, 11.4, 3.0, [
        ("The fraction of the run spent within 3 dB of the perfect-knowledge beamformer.",
         "Because it is normalised per cell, a 6-element patch array and a 20-element dipole "
         "array can be scored on one axis — each against its own potential, not against each "
         "other."),
        ("An array that can do nothing scores 100 for doing nothing,",
         "which is correct: the 1-element array's quiescent beam is its optimum, and the "
         "metric says so rather than marking it failed."),
        ("90 is the pass mark.",
         "A judgement, not a derivation — stated here so the pass counts elsewhere in the "
         "deck can be re-read against a different one."),
    ])
    footnote(s, "Availability, the previous headline metric, saturates on easy cells: it moved "
                "+0.05 pp for a change worth 35 points here. It cannot rank algorithms.")


def slide_drift_result(prs):
    s = blank(prs)
    chrome(s, 7, "Result 1", "Static is fine. Drifting is not.")
    table(s, 0.72, 2.15, 8.6, [
        ["Scenario", "lcmv", "passing", "predict", "passing"],
        ["STATIC", "98.3", ("15 / 15", PASS, True), "98.3", ("15 / 15", PASS, True)],
        ["DRIFT", "55.2", ("2 / 15", FAIL, True), "55.0", ("2 / 15", FAIL, True)],
        ["ON/OFF", "78.1", ("8 / 15", FAIL, True), "78.1", ("8 / 15", FAIL, True)],
    ], [2.6, 1.5, 1.5, 1.5, 1.5], row_h=0.44)
    prose(s, 0.72, 4.3, 11.4, 2.2, [
        [("Thirteen of fifteen drift cells fail, on all three arrays. And note the second "
          "column pair: ", False),
         ("predict scores the same as lcmv, to a tenth of a point", True),
         (" — the algorithm that was supposed to handle a moving jammer was contributing "
          "nothing. That is the thread this campaign pulled.", False)],
        [("The ON/OFF row hides the same problem in a different place, and it took a "
          "second campaign to see it.", True)],
    ])


def slide_diagnosis(prs):
    s = blank(prs)
    chrome(s, 8, "Result 2 — diagnosis",
           "Aiming one covariance horizon into the past")
    picture(s, "lag_sweep.png", 0.72, 2.35, 6.3, 3.48)
    caption(s, 0.72, 5.95, 6.3, 0.6,
            "Steering a hard null from scenario truth, delayed by τ steps. At τ = 10 it "
            "reproduces the measured failure exactly.")
    bullets(s, 7.35, 2.35, 5.2, 3.6, [
        ("The forgetting factor sets the lag.",
         "λ = 0.90 means R is an average over the last ~10 steps, so the direction it "
         "encodes is ~10 steps old."),
        ("At 10°/s that is 30° of angle error,",
         "which is several beamwidths on these arrays — the null lands where the jammer "
         "was, not where it is."),
        ("The test was constructed to be decisive.",
         "Feeding the null a deliberately delayed truth reproduces the failure curve; no "
         "other explanation was needed."),
    ], size=14.5)


def slide_lag_budget(prs):
    s = blank(prs)
    chrome(s, 9, "Result 2 — the lag budget", "The lead is derived, not tuned")
    prose(s, 0.72, 2.1, 11.4, 0.8, [
        "MUSIC is computed from the same exponentially-forgotten covariance, so it inherits "
        "that horizon. The rest is the delay between deciding a weight and applying it.",
    ])
    rect(s, 0.72, 3.0, 11.8, 0.95, PANEL)
    tf = textbox(s, 1.0, 3.31, 11.2, 0.5)
    run(para(tf, first=True), "lead = 10 (covariance) + 4 (application) = 14 steps",
        18, MONO, False, INK)
    bullets(s, 0.72, 4.3, 11.4, 2.0, [
        ("No free parameter was fitted.",
         "Both terms come from configuration that was already fixed: λ gives the first, the "
         "pipeline depth gives the second. A tuned lead would have been a different kind of "
         "result — and a much weaker one."),
        ("It is capped and gated, not applied blindly.",
         "The lead is scaled by the observed toggle period where one exists, and the "
         "predictor is gated on estimated angular speed so a static jammer never sees it."),
    ])


def slide_cv_fix(prs):
    s = blank(prs)
    chrome(s, 10, "Result 3 — the fix", "A constant-velocity tracker on the jammer's angle")
    picture(s, "cv_tracker.png", 0.72, 2.05, 7.39, 3.45)
    caption(s, 0.72, 5.62, 7.39, 0.6,
            "A Kalman filter on [θ, θ̇, φ, φ̇], fed by the MUSIC angle the algorithm already "
            "computes, predicting 14 steps ahead.")
    bullets(s, 8.45, 2.05, 4.1, 3.6, [
        ("Nothing new is measured.", "The filter consumes an estimate the stack already "
         "produced; no extra observation is required of the receiver."),
        ("Mirror folding is handled explicitly.", "On symmetric arrays θ and 180°−θ are "
         "indistinguishable; the filter detects the degeneracy and folds rather than "
         "chasing a phantom."),
        ("It self-disables.", "Below a minimum estimated speed the prediction is declared "
         "invalid and the reactive solution stands."),
    ], size=13.5, gap=8)


def slide_cv_campaign(prs):
    s = blank(prs)
    chrome(s, 11, "Result 3 — campaign", "+35 points on drift, and nothing else moves")
    table(s, 0.72, 2.15, 8.6, [
        ["Scenario", "lcmv", "predict", "predict + CV", "Δ"],
        ["STATIC", "98.3", "98.3", "98.3", ("0.0", INK3, False)],
        ["DRIFT", "55.2", "55.0", ("90.1", INK, True), ("+35.1", PASS, True)],
        ["ON/OFF", "78.1", "78.1", "78.1", ("0.0", INK3, False)],
    ], [2.6, 1.5, 1.5, 1.5, 1.5], row_h=0.44)
    prose(s, 0.72, 4.3, 11.4, 2.2, [
        [("The zeros are not luck.", True),
         (" The predictor is gated on estimated angular speed, so on a static or a switching "
          "jammer it never engages and the weights are bit-for-bit the reactive ones. A fix "
          "that improves its target case and provably does not touch the others is the "
          "cheapest kind to ship.", False)],
        [("Drifting-jammer cells at the pass mark: 2 → 11 of 15.", True)],
    ])


def slide_calibration(prs):
    s = blank(prs)
    chrome(s, 12, "Result 4 — the real constraint",
           "Calibration, not drift, is what limits fielding")
    picture(s, "calibration_sweep.png", 0.72, 2.3, 6.3, 3.48)
    caption(s, 0.72, 5.9, 6.3, 0.7,
            "Per-element calibration error, 12 random draws per level, median with p10–p90 "
            "band. Until this phase the simulator computed signal power from the same "
            "steering vector it handed the beamformer, so this was unmeasurable.")
    bullets(s, 7.35, 2.3, 5.2, 3.6, [
        ("Usable below ~1° RMS per-element phase error;", "by 3° everything is at the floor."),
        ("The spread is enormous.", "At 2° the p10–p90 range is 1.1–76.8. One array build "
         "could be fine and the next unusable."),
        ("Unchanged by the predictor,", "which keeps its 2–4× advantage throughout. This is "
         "an exposure of the whole approach, not of one algorithm."),
        ("Slide 17 measures the same weakness from the other direction",
         "— and finds it larger than this chart alone suggests."),
    ], size=14, gap=9)


def slide_onoff_never_ran(prs):
    s = blank(prs)
    chrome(s, 13, "Phase O — diagnosis",
           "The algorithm built for on/off jammers was never running")
    stats = [("≤ 1.3%", "of steps where the anticipatory branch fired, before the fix", FAIL),
             ("0.46 → 0.13", "presence error against the true duty cycle, after", INK),
             ("2", "independent blockers, both arithmetic", INK)]
    for i, (v, l, c) in enumerate(stats):
        statbox(s, 0.72 + i * 4.03, 2.15, 3.75, 1.35, v, l, c, value_size=29)
    bullets(s, 0.72, 3.85, 11.4, 2.9, [
        ("The presence test read a covariance that could not see the toggle.",
         "Detection ran on the λ = 0.90 covariance, whose 10-step memory smears an on/off "
         "edge into nothing. A second, fast-forgetting covariance is now kept purely for "
         "presence."),
        ("The anticipatory lead was shorter than the toggle it was anticipating.",
         "The fixed 14-step lead is now scaled by the observed period and capped, so it "
         "reaches across an OFF window instead of falling inside it."),
        ("Neither was a tuning problem, and neither was visible in the score.",
         "The ON/OFF row had sat at 78.1 for three algorithms across two campaigns — "
         "identical numbers that should have been read as a signal much earlier."),
    ])


def slide_every_array(prs):
    s = blank(prs)
    chrome(s, 14, "Phase O — coverage", "Every array now runs, and the impossible ones say so")
    table(s, 0.72, 2.15, 11.8, [
        ["Array", "elements", "before", "after", "note"],
        ["spacing0.6", "16", ("ran", INK3, False), "80.1", "hardest geometry in the suite"],
        ["ManyDipoles", "20", ("ran", INK3, False), "81.0", "mirror-degenerate"],
        ["Monopoles", "14", ("ran", INK3, False), "82.6", ""],
        ["patchs_with_monopoles", "6", ("ran", INK3, False), "80.4", "very broad beam"],
        ["Dipole", "1", ("ERROR", FAIL, True), ("100.0", PASS, True),
         "no degrees of freedom — doing nothing is optimal"],
        ["patch_back2back", "2", ("ERROR", FAIL, True), ("50.6", PASS, True),
         "MUSIC infeasible; degrades to reactive"],
    ], [3.3, 1.3, 1.3, 1.3, 4.6], row_h=0.40, size=13, left_cols=(4,))
    for gf in [sh for sh in s.shapes if sh.has_table]:
        for ri in range(1, 7):
            p = gf.table.cell(ri, 4).text_frame.paragraphs[0]
            p.alignment = PP_ALIGN.LEFT
            for r in p.runs:
                r.font.bold = False
                r.font.size = Pt(12)
                r.font.color.rgb = INK3
    prose(s, 0.72, 5.4, 11.4, 1.2, [
        [("Two arrays used to crash the whole stack.", True),
         (" An opt-in loading feature threw instead of degrading, and MUSIC's hardcoded "
          "source count threw on any array with too few elements. Both now warn and fall "
          "back. A geometry preflight additionally refuses 11 of 35 (array, target) pairs "
          "as not real tests, each with a stated reason — rather than scoring them as "
          "failures.", False)],
    ])


def slide_reversal(prs):
    s = blank(prs)
    chrome(s, 15, "Phase O2 — the reversal",
           "The earlier recommendation was wrong, and measurement said so")
    table(s, 0.72, 2.15, 7.4, [
        ["", "mean score", "cells passing"],
        ["no repair", "73.5", "24 / 90"],
        ["+ graded release", "76.1", "29 / 90"],
        [("+ binary repair — ships", INK, True), ("80.6", INK, True), ("32 / 90", PASS, True)],
        ["best-of-both (oracle)", "81.5", "32 / 90"],
    ], [3.4, 2.0, 2.0], row_h=0.42)
    bullets(s, 8.35, 2.15, 4.2, 4.0, [
        ("Phase O measured a regression and advised against enabling.",
         "That was a run-length artifact: period detection needs 3 cycles and those runs "
         "were 4 long."),
        ("Re-run at 8 cycles it improves or ties on every array.",
         "Now enabled by default."),
        ("The graded release was built, measured and rejected.",
         "It rescues the binary tail but loses on 45 of 90 cells."),
    ], size=13.5, gap=8)
    footnote(s, "Best-of-both is an oracle no real switching rule could beat, so the entire "
                "family of release policies is bounded at +0.9 mean and zero extra passing "
                "cells. That bound is the useful output of having built the graded one.")


def slide_videos(prs):
    s = blank(prs)
    chrome(s, 16, "Evidence", "What the repair looks like, cell by cell")
    table(s, 0.72, 2.15, 11.8, [
        ["Video", "array / geometry", "lcmv", "predict", "+ repair"],
        ["A_easy", "patchs_with_monopoles, T = 25 s", "96.8", "96.8",
         ("97.9", PASS, True)],
        ["B_repair_pays", "spacing0.6 (45,150), T = 10 s", "26.2", "30.4",
         ("58.3", INK, True)],
        ["C_hard", "spacing0.6_disturbed3, T = 4 s, close jammer", "35.4", "35.4",
         ("62.7", INK, True)],
        ["D_mirror", "ManyDipoles — mirror-degenerate", "93.1", "93.3",
         ("96.8", PASS, True)],
        ["E_infeasible", "patch_back2back, 2 elements", "98.0", "98.0", "98.0"],
    ], [2.5, 5.1, 1.4, 1.4, 1.4], row_h=0.40, size=13, left_cols=(1,))
    for gf in [sh for sh in s.shapes if sh.has_table]:
        for ri in range(1, 6):
            p = gf.table.cell(ri, 1).text_frame.paragraphs[0]
            p.alignment = PP_ALIGN.LEFT
            for r in p.runs:
                r.font.bold = False
                r.font.size = Pt(12)
                r.font.color.rgb = INK3
    prose(s, 0.72, 4.85, 11.4, 1.8, [
        [("Each video shows one radiation pattern per algorithm on a ", False),
         ("shared", True),
         (" colour scale, above a common output-SINR trace with the jammer-ON windows "
          "shaded. The shared scale is the point: per-panel autoscaling would make a "
          "beamformer that has thrown away 10 dB of gain look identical to one that has "
          "not.", False)],
        [("E_infeasible is a null result on purpose", True),
         (" — MUSIC cannot run on two elements, and all three columns agreeing is the "
          "correct outcome.", False)],
    ])


def slide_amplitude(prs):
    s = blank(prs)
    chrome(s, 17, "Amplitude sweep — the last measurement",
           "The gap that is left is a strong-signal gap")
    table(s, 0.72, 2.15, 7.2, [
        ["Array", "oracle gains", "we capture", "lost"],
        ["spacing0.6", "20.0", "11.9", ("8.1", FAIL, True)],
        ["patchs_with_monopoles", "20.0", "17.3", "2.7"],
        ["ManyDipoles", "20.0", "19.7", ("0.3", PASS, True)],
        ["Monopoles", "20.0", "19.8", ("0.2", PASS, True)],
    ], [3.0, 1.4, 1.4, 1.4], row_h=0.42)
    caption(s, 0.72, 4.4, 7.2, 0.5,
            "dB gained when the wanted signal goes from 0 to 20 dB above the noise floor.")
    bullets(s, 8.15, 2.15, 4.4, 4.2, [
        ("Raising the jammer 20 dB is nearly free.",
         "At most 2.1 dB of the achievable value. The nulling works."),
        ("Raising the wanted signal 20 dB is not.",
         "MPDR self-cancellation: the wanted signal sits inside R, so under steering "
         "mismatch the solver spends degrees of freedom cancelling it."),
        ("The oracle is immune,", "because it is handed the true steering vector."),
    ], size=13.5, gap=8)
    prose(s, 0.72, 5.05, 7.2, 1.6, [
        [("Two honest caveats.", True),
         (" The worst case is the highest-directivity array, which fits the mismatch story "
          "— a narrow beam pays more for the same angular error — but four arrays at one "
          "geometry each is an observation, not a controlled sweep. And it means every "
          "other number in this deck is a weak-signal number.", False)],
    ], size=12.5)


def slide_status(prs):
    s = blank(prs)
    chrome(s, 18, "Status", "Where this leaves the milestone")
    tf = textbox(s, 0.72, 2.1, 5.7, 0.3)
    run(para(tf, first=True), "CLOSED THIS PHASE", 11, SANS, True, ACCENT)
    bullets(s, 0.72, 2.5, 5.7, 3.8, [
        ("Drifting jammer —", "diagnosed and fixed, +35 points"),
        ("On/off predictor —", "repaired, and enabled by default"),
        ("Release policy —", "closed, bounded at +0.9"),
        ("Every array runs —", "and the impossible ones report why"),
        ("Calibration mismatch —", "now measurable"),
        ("Amplitude sweep —", "run; it is what surfaced the finding on slide 17"),
    ], size=13.5, gap=7)
    tf = textbox(s, 6.9, 2.1, 5.7, 0.3)
    run(para(tf, first=True), "STILL OPEN", 11, SANS, True, FAIL)
    bullets(s, 6.9, 2.5, 5.7, 3.8, [
        ("Steering mismatch at strong signal —", "up to 8.1 dB unclaimed; the largest "
         "named number left"),
        ("Why the loss ranges 8.1 → 0.2 dB across arrays —", "needs a mismatch × σs sweep"),
        ("guard_deg is still a global 5° in config —", "derive it per array"),
        ("Not qualified above ~10°/s drift —", "the predictor extends the envelope, it "
         "does not remove the limit"),
        ("Multiple jammers; polarization mismatch; real-time compute", ""),
    ], size=13.5, gap=7)
    footnote(s, "The CV drift predictor ships disabled; the on/off repair ships enabled after "
                "the O2 reversal. 59 of 59 anti-jam gates pass.  Next: steering-mismatch "
                "robustness — worth up to 8.1 dB, against the +0.9 dB a perfect release "
                "policy was worth.", top=6.5)


def main():
    prs = Presentation()
    prs.slide_width = Emu(12191695)
    prs.slide_height = Emu(6858000)
    for fn in (slide_title, slide_conclusions, slide_scope, slide_algorithms,
               slide_how, slide_metric, slide_drift_result, slide_diagnosis,
               slide_lag_budget, slide_cv_fix, slide_cv_campaign, slide_calibration,
               slide_onoff_never_ran, slide_every_array, slide_reversal,
               slide_videos, slide_amplitude, slide_status):
        fn(prs)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    prs.save(OUT)
    print("wrote %s (%d slides)" % (OUT, len(prs.slides)))


if __name__ == "__main__":
    main()
