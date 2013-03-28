from matplotlib.path import Path
from matplotlib.patches import BoxStyle
import matplotlib.pyplot as plt

# we may derive from matplotlib.patches.BoxStyle._Base class.
# You need to overide transmute method in this case.

class MyStyle(BoxStyle._Base):
    def __init__(self, pad=0.3):
        self.pad = pad
        super(MyStyle, self).__init__()

    def transmute(self, x0, y0, width, height, mutation_size):
        """
        Given the location and size of the box, return the path of
        the box around it.

         - *x0*, *y0*, *width*, *height* : location and size of the box
         - *mutation_size* : a reference scale for the mutation.

        Often, the *mutation_size* is the font size of the text.
        You don't need to worry about the rotation as it is
        automatically taken care of.
        """

        # padding
        pad = mutation_size * self.pad

        # width and height with padding added.
        width, height = width + 2.*pad, \
                        height + 2.*pad,

        # boundary of the padded box
        x0, y0 = x0-pad, y0-pad,
        x1, y1 = x0+width, y0 + height

        cp = [(x0, y0),
              (x1, y0), (x1, y1), (x0, y1),
               (x0, y0),
              (x0, y0)]

        com = [Path.MOVETO,
               Path.LINETO, Path.LINETO, Path.LINETO,
               Path.LINETO,
               Path.CLOSEPOLY]

        path = Path(cp, com)

        return path

def text(x,y,rot,t,size=8):
    # register the custom style
    BoxStyle._style_list["angled"] = MyStyle

    plt.gca().text(x,y, t, size=size, va="center", ha="center", rotation=rot,
            bbox=dict(boxstyle="angled,pad=0.0",facecolor='white',edgecolor='white'))

    del BoxStyle._style_list["angled"]

if __name__=='__main__':
    plt.figure(1, figsize=(3,3))
    plt.subplot(111)
    plt.plot([0,1],[0,1])
    text(0.5,0.5,45,'halloj',size=14)
    plt.show()
