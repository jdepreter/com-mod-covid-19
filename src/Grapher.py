import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

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
    def __init__(self, days: int, y: list):
        """
        :param days: days of the simulation
        :param y: list containing list of y values
        y = [ ndarray-susceptible, ndarray-infected ... ]
        """
        self.days = days
        self.y = y

    def show(self):
        """
        Plot the graph
        :return: nothing
        """
        x = np.arange(0, self.days, 1)

        fig, ax = plt.subplots()

        for yvalues in self.y:
            line, = ax.plot(x, yvalues, color='k')

            def update(num, x, y, line):
                line.set_data(x[:num], y[:num])
                line.axes.axis([0, 10, 0, 20])
                return line,

            ani = animation.FuncAnimation(fig, update, len(x), fargs=[x, yvalues, line],
                                          interval=200, blit=True)
            # ani.save('test.gif')
        mywriter = animation.FFMpegFileWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        ani.save('myAnimation.mp4', writer=mywriter)

        plt.show()

# Test:
days = 10
y = [np.arange(0, days, 1), np.arange(0, days*2, 2) ]
g = Grapher(days, y)
# print(y[:5])
# print(np.pad(y[:5], (0, days-5), mode='constant', constant_values=0))
g.show()