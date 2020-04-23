import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# To save the animation, use e.g.
#
# ani.save("movie.mp4")
#
# or
#
# from matplotlib.animation import FFMpegWriter
# writer = FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("movie.mp4", writer=writer)
class Grapher:
    def __init__(self, days: int, y):
        self.days = days
        self.y = y

    def show(self):
        x = np.arange(0, self.days, 1)

        fig, ax = plt.subplots()
        line, = ax.plot(x, self.y, color='k')

        def update(num, x, y, line):
            line.set_data(x[:num], y[:num])
            line.axes.axis([0, 10, 0, 10])
            return line,

        ani = animation.FuncAnimation(fig, update, len(x), fargs=[x, y, line],
                                      interval=50, blit=True)
        # ani.save('test.gif')
        plt.show()

# Test:
# days = 10
# y = np.arange(0, days, 1)
# g = Grapher(days, y)
# print(y[:5])
# print(np.pad(y[:5], (0, days-5), mode='constant', constant_values=0))
# g.show()