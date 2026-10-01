"""Regenerate the photographs of this example.

They are generated rather than photographed, and that is the whole point: the shear is put
in by this script, so after running the case you can check the answer instead of trusting
it. Nothing here imports PIV-NP -- it only writes PNG files -- so it is readable on its own
and easy to change.

    python examples/piv-from-images/make_example.py

Change ``SLIP_PER_STEP`` to shear the block harder, or ``STEPS`` to make the test longer,
and remember to change ``total_steps`` in the ``.PAR`` to match.
"""

from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent

WIDTH, HEIGHT = 256, 192
#: the block of soil in the first photograph: left, right, top, bottom, in pixels
LEFT, RIGHT, TOP, BOTTOM = 32, 224, 32, 160
STEPS = 10
#: how many pixels further the top of the block slides than the base, per step. Kept small
#: on purpose: a displacement above about a quarter of the interrogation window is more than
#: a correlation can follow.
SLIP_PER_STEP = 2.0

#: what the ``.PIV`` declares, repeated here only so this script can print the answer in mm
SCALE = 0.0005


def speckle(seed: int = 20251001) -> np.ndarray:
    """Blurred noise: the texture of sand, and something a correlation can lock onto.

    The blur matters. White noise has a correlation peak one pixel wide, and a three-point
    sub-pixel fit has nothing to fit on it -- a PIV run on white noise returns whole pixels
    and looks broken. Real photographs of soil are never white noise.
    """
    rng = np.random.default_rng(seed)
    # built at twice the size and halved, which blurs and avoids a visible pixel grid
    rough = rng.normal(0.0, 1.0, (HEIGHT * 2, WIDTH * 2))
    weights = np.array([1.0, 4.0, 7.0, 4.0, 1.0])
    weights /= weights.sum()
    out = rough
    for axis in (0, 1):
        length = out.shape[axis]
        edge = np.clip(np.arange(-2, length + 2), 0, length - 1)
        padded = np.take(out, edge, axis=axis)
        blurred = np.zeros_like(out)
        for offset, weight in enumerate(weights):
            blurred += weight * np.take(padded, np.arange(offset, offset + length), axis=axis)
        out = blurred
    out = out[::2, ::2]
    out = (out - out.min()) / (out.max() - out.min())
    return 40.0 + 200.0 * out          # 40 to 240, well clear of the black background


def shear_after(step: int) -> float:
    """Accumulated shear after ``step`` steps: the engineering strain du/dy."""
    return step * SLIP_PER_STEP / (BOTTOM - TOP)


def render(reference: np.ndarray, step: int) -> np.ndarray:
    """The photograph at a step: the block sheared, on a black background.

    Simple shear, so the answer is one number everywhere. A point at height ``y`` has slid
    sideways by ``shear x (BOTTOM - y)``: the base of the block stays where it is and the top
    moves most, as soil on a slope does. Each photograph is built by asking, for every one of
    its pixels, where in the undeformed block that pixel came from -- which is the way round
    that leaves no gaps.
    """
    shear = shear_after(step)
    y, x = np.mgrid[0:HEIGHT, 0:WIDTH].astype(np.float64)
    came_from_x = x - shear * (BOTTOM - y)
    came_from_y = y

    inside = ((came_from_x >= LEFT) & (came_from_x <= RIGHT - 1)
              & (came_from_y >= TOP) & (came_from_y <= BOTTOM - 1))
    x0 = np.clip(np.floor(came_from_x).astype(int), 0, WIDTH - 2)
    y0 = np.clip(np.floor(came_from_y).astype(int), 0, HEIGHT - 2)
    fx, fy = came_from_x - x0, came_from_y - y0
    value = ((1 - fy) * ((1 - fx) * reference[y0, x0] + fx * reference[y0, x0 + 1])
             + fy * ((1 - fx) * reference[y0 + 1, x0] + fx * reference[y0 + 1, x0 + 1]))
    return np.where(inside, value, 0.0)


def main() -> None:
    folder = HERE / "images"
    folder.mkdir(parents=True, exist_ok=True)
    reference = speckle()
    for number in range(1, STEPS + 2):          # one photograph more than there are steps
        frame = render(reference, number - 1)
        Image.fromarray(np.round(frame).astype(np.uint8), mode="L").save(
            folder / f"shear_{number:03d}.png")

    total = shear_after(STEPS)
    print(f"{STEPS + 1} photographs of {WIDTH}x{HEIGHT} written to {folder}")
    print(f"shear per step                 : {shear_after(1):.6f}")
    print(f"total shear after {STEPS} steps      : {total:.6f}")
    print(f"the top of the block slides    : {SLIP_PER_STEP * STEPS:.1f} px "
          f"= {SLIP_PER_STEP * STEPS * SCALE * 1000:.2f} mm")
    print(f"equivalent strain it implies   : {total / np.sqrt(3):.6f}")


if __name__ == "__main__":
    main()
