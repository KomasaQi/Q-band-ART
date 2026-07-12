from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, JpegImagePlugin  # noqa: F401


HERE = Path(__file__).resolve().parent
BACKGROUND = HERE / "qband_scheduler_bg.png"
OUT_PNG = HERE / "qband_scheduler.png"
OUT_PDF = HERE / "qband_scheduler.pdf"


def load_font(name: str, size: int) -> ImageFont.FreeTypeFont:
    candidates = [
        Path("C:/Windows/Fonts") / name,
        Path("C:/Windows/Fonts/arial.ttf"),
        Path("C:/Windows/Fonts/calibri.ttf"),
    ]
    for path in candidates:
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default(size=size)


def hex_rgba(color: str, alpha: int) -> tuple[int, int, int, int]:
    color = color.lstrip("#")
    return tuple(int(color[i : i + 2], 16) for i in (0, 2, 4)) + (alpha,)


def arrow(draw: ImageDraw.ImageDraw, start, end, fill, width: int = 5):
    draw.line([start, end], fill=fill, width=width)
    x1, y1 = start
    x2, y2 = end
    dx, dy = x2 - x1, y2 - y1
    length = max((dx * dx + dy * dy) ** 0.5, 1.0)
    ux, uy = dx / length, dy / length
    px, py = -uy, ux
    head = width * 4.2
    p1 = (x2, y2)
    p2 = (x2 - head * ux + 0.55 * head * px, y2 - head * uy + 0.55 * head * py)
    p3 = (x2 - head * ux - 0.55 * head * px, y2 - head * uy - 0.55 * head * py)
    draw.polygon([p1, p2, p3], fill=fill)


def rounded_label(
    draw: ImageDraw.ImageDraw,
    box,
    title: str,
    lines: list[str],
    accent: str,
    title_font: ImageFont.FreeTypeFont,
    body_font: ImageFont.FreeTypeFont,
):
    x, y, w, h = box
    radius = int(min(w, h) * 0.12)
    draw.rounded_rectangle(
        [x, y, x + w, y + h],
        radius=radius,
        fill=(255, 255, 255, 232),
        outline=hex_rgba(accent, 230),
        width=max(2, int(w * 0.01)),
    )
    draw.rounded_rectangle(
        [x, y, x + int(w * 0.04), y + h],
        radius=radius,
        fill=hex_rgba(accent, 235),
    )
    pad_x = int(w * 0.085)
    pad_y = int(h * 0.16)
    draw.text((x + pad_x, y + pad_y), title, font=title_font, fill=hex_rgba(accent, 255))
    body_y = y + pad_y + int(h * 0.34)
    for line in lines:
        draw.text((x + pad_x, body_y), line, font=body_font, fill=(34, 43, 53, 255))
        body_y += int(body_font.size * 1.25)


def main():
    base = Image.open(BACKGROUND).convert("RGBA")
    w, h = base.size
    overlay = Image.new("RGBA", base.size, (255, 255, 255, 0))
    draw = ImageDraw.Draw(overlay)

    title_font = load_font("arialbd.ttf", max(25, int(w / 62)))
    body_font = load_font("arial.ttf", max(19, int(w / 86)))
    small_font = load_font("arial.ttf", max(17, int(w / 98)))

    # Slight veil behind annotations keeps the generated details visible but readable.
    draw.rectangle([0, 0, w, h], fill=(255, 255, 255, 30))

    labels = [
        (
            (0.395, 0.050, 0.235, 0.128),
            "Offline Q-band library",
            ["bands B_i | Q_i | Ebar_i", "target weights w_i^tar"],
            "#1b5e91",
        ),
        (
            (0.070, 0.222, 0.210, 0.120),
            "Preview path risk",
            ["kappa_ref -> a_y,ref", "look-ahead excitation"],
            "#2b7fc5",
        ),
        (
            (0.293, 0.392, 0.190, 0.118),
            "Predicted response",
            ["A_y(U) from model", "ESO residual d_hat"],
            "#143e6e",
        ),
        (
            (0.516, 0.392, 0.198, 0.128),
            "Band-energy risk",
            ["R_i^pre = [E_i/Ebar_i - 1]+", "danger band -> penalty"],
            "#c85f16",
        ),
        (
            (0.722, 0.315, 0.200, 0.146),
            "Dynamic scheduling",
            ["lambda_i = preview + feedback", "w_i(t) rate-smoothed"],
            "#087d7b",
        ),
        (
            (0.845, 0.514, 0.130, 0.105),
            "NMPC steering",
            ["delta_f command", "Q-band penalty"],
            "#087d7b",
        ),
        (
            (0.420, 0.735, 0.350, 0.128),
            "Feedback monitor and memory",
            ["band gain G_j^fb", "memory mu_j(t) decays"],
            "#0a7775",
        ),
    ]

    for rx, title, lines, accent in labels:
        x, y, bw, bh = rx
        rounded_label(
            draw,
            (int(x * w), int(y * h), int(bw * w), int(bh * h)),
            title,
            lines,
            accent,
            title_font,
            body_font,
        )

    dark_blue = hex_rgba("#113b6c", 235)
    teal = hex_rgba("#087d7b", 235)
    orange = hex_rgba("#c85f16", 235)

    # Callout arrows for the logic that is not fully explicit in the generated image.
    arrow(draw, (int(0.510 * w), int(0.180 * h)), (int(0.605 * w), int(0.365 * h)), dark_blue, max(4, int(w / 380)))
    arrow(draw, (int(0.267 * w), int(0.302 * h)), (int(0.312 * w), int(0.428 * h)), dark_blue, max(4, int(w / 400)))
    arrow(draw, (int(0.708 * w), int(0.455 * h)), (int(0.752 * w), int(0.405 * h)), orange, max(4, int(w / 390)))
    arrow(draw, (int(0.900 * w), int(0.610 * h)), (int(0.792 * w), int(0.735 * h)), teal, max(4, int(w / 390)))
    arrow(draw, (int(0.758 * w), int(0.735 * h)), (int(0.180 * w), int(0.735 * h)), teal, max(4, int(w / 390)))

    # Small explanatory strip: concise enough to remain legible in a two-column figure.
    strip = [int(0.295 * w), int(0.900 * h), int(0.865 * w), int(0.962 * h)]
    draw.rounded_rectangle(strip, radius=int(0.014 * w), fill=(255, 255, 255, 232), outline=teal, width=max(2, int(w / 600)))
    strip_text = "Only preview-excited or feedback-amplified bands are weighted up; release is slowed by residual-memory decay."
    draw.text((strip[0] + int(0.018 * w), strip[1] + int(0.018 * h)), strip_text, font=small_font, fill=(28, 45, 55, 255))

    composed = Image.alpha_composite(base, overlay).convert("RGB")
    composed.save(OUT_PNG, optimize=True, quality=95)
    composed.save(OUT_PDF, resolution=300.0)
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
