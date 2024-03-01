from matplotlib.pylab import gca, figure, plot, subplot, title, xlabel, ylabel, xlim,show

def draw_plot(data,plot_title):
    plot(range(len(data)),data,alpha=0.8,color='red')
    title(plot_title)
    xlabel("Samples")
    ylabel("Signal")
    xlim((0,len(data)-1))

def draw_segments(segments):
    ax = gca()
    for segment in segments:
        line = Line2D((segment[0],segment[2]),(segment[1],segment[3]))
        ax.add_line(line)

def draw_xvline_by_segments(segments):
    for seg in segments:
        x0 = seg[0]
        plt.axvline(x=x0)
