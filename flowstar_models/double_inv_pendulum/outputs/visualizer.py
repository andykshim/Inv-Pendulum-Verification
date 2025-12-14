import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.widgets import Slider
import re

def parse_plt_file(path):
    """
    Parse a Flow* GNUPLOT .plt file containing multiple octagons
    separated by blank lines. Returns:

        frames: list of list-of-(x,y) tuples
        xlabel: string or None
        ylabel: string or None
    """

    frames = []
    current = []

    xlabel = None
    ylabel = None

    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()

            # Extract xlabel/ylabel if present
            if line.startswith("set xlabel"):
                # set xlabel "x"
                try:
                    xlabel = line.split('"')[1]
                except:
                    pass
                continue

            if line.startswith("set ylabel"):
                try:
                    ylabel = line.split('"')[1]
                except:
                    pass
                continue

            # Skip gnuplot commands
            if not line or line.startswith("set ") or line.startswith("plot"):
                # Blank line → finish current frame
                if line == "" and current:
                    frames.append(current)
                    current = []
                continue

            # Parse numeric data lines: two floating point numbers
            try:
                x, y = map(float, line.split())
                current.append((x, y))
            except:
                # Skip anything else
                continue

    # Last frame if missing trailing blank line
    if current:
        frames.append(current)

    return frames, xlabel, ylabel

class SimpleVisualizer:
    def __init__(self, plt_file):
        self.polygons, self.xlabel, self.ylabel = parse_plt_file(plt_file)
        self.num_frames = len(self.polygons)
        self.cur = 0

        all_pts = np.vstack(self.polygons)
        self.min_x, self.max_x = all_pts[:,0].min(), all_pts[:,0].max()
        self.min_y, self.max_y = all_pts[:,1].min(), all_pts[:,1].max()
        pad_x = 0.1*(self.max_x-self.min_x)
        pad_y = 0.1*(self.max_y-self.min_y)
        self.min_x -= pad_x; self.max_x += pad_x
        self.min_y -= pad_y; self.max_y += pad_y

        # create figure
        self.fig, self.ax = plt.subplots(figsize=(8,6))
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)

        # slider (optional but helpful)
        slider_ax = self.fig.add_axes([0.15, 0.02, 0.7, 0.03])
        self.slider = Slider(slider_ax, "Frame", 0, self.num_frames - 1,
                             valinit=0, valstep=1)
        self.slider.on_changed(self.on_slider)

        self.draw()

    def draw(self):
        self.ax.clear()
        poly = self.polygons[self.cur]
        self.ax.add_patch(Polygon(poly, closed=True, facecolor='lightblue', edgecolor='blue'))
        self.ax.set_xlim(self.min_x, self.max_x)
        self.ax.set_ylim(self.min_y, self.max_y)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_title(f"Frame {self.cur+1}/{self.num_frames}")
        self.fig.canvas.draw_idle()

    def on_key(self, event):
        if event.key == 'right':
            self.cur = min(self.cur + 1, self.num_frames - 1)
            self.draw()
        elif event.key == 'left':
            self.cur = max(self.cur - 1, 0)
            self.draw()

    def on_slider(self, val):
        self.cur = int(val)
        self.draw()

    def show(self):
        plt.show()


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("plt_file")
    args = parser.parse_args()
    SimpleVisualizer(args.plt_file).show()

if __name__ == "__main__":
    main()
