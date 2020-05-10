import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from os.path import dirname, abspath
img_folder = os.path.join(dirname(dirname(abspath(__file__))), 'img')
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


def animate(y: list, labels: list, save=False, display=True, y_label='', x_label='Days since first case', name='test'):
    """
    Plot the graph
    :return: nothing
    """
    fig = plt.figure()
    ax = plt.axes()

    lines = [plt.plot([], [], label=labels[_])[0] for _ in range(len(y))]  # lines to animate

    def init():
        # init lines
        for line in lines:
            line.set_data([], [])

        return lines  # return everything that must be updated

    def animate(i):
        # animate lines

        for j, line in enumerate(lines):
            line.set_data(range(i), y[j][:i])

        ax.relim()
        ax.autoscale_view()
        return lines, ax  # return everything that must be updated

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(y[0]), interval=100, blit=False)

    # mywriter = animation.FFMpegFileWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    plt.legend(loc="upper left")
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if save:
        anim.save('%s/%s.gif' % (img_folder, name), writer='imagemagick', fps=10)

    if display:
        plt.show()


def plot(y: list, labels: list, name='temp', y_label=''):
    plt.clf()
    for i, line in enumerate(y):
        plt.plot(line, label=labels[i])

    plt.xlabel('Days since first case')
    plt.ylabel(y_label)
    plt.title(name, pad=15)
    plt.legend()
    plt.savefig(img_folder + '/' + name)

# Test:
# days = 10
# y = [np.arange(0, days, 1), np.arange(0, days * 2, 2)]
# g = Grapher(days, y)
# # print(y[:5])
# # print(np.pad(y[:5], (0, days-5), mode='constant', constant_values=0))
# g.animate("test")
