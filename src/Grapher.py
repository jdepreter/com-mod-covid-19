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

    def animate(self):
        """
        Plot the graph
        :return: nothing
        """
        fig = plt.figure()
        ax = plt.axes(xlim=(0, self.days), ylim=(0, 10000))

        lines = [plt.plot([], [])[0] for _ in range(len(self.y))]  # lines to animate

        def init():
            # init lines
            for line in lines:
                line.set_data([], [])

            return lines  # return everything that must be updated

        def animate(i):
            # animate lines

            for j, line in enumerate(lines):
                line.set_data(range(i), self.y[j][:i])
                print(i)

            return lines  # return everything that must be updated

        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=len(self.y[0]), interval=100, blit=True)

        # mywriter = animation.FFMpegFileWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        anim.save('animation.gif', writer='imagemagick', fps=24)
        # plt.show()


# Test:
# days = 10
# y = [np.arange(0, days, 1), np.arange(0, days * 2, 2)]
# g = Grapher(days, y)
# # print(y[:5])
# # print(np.pad(y[:5], (0, days-5), mode='constant', constant_values=0))
# g.animate()
