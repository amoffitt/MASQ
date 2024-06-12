import matplotlib
import matplotlib.pyplot as plt

WIDTH = 0.25


def plot_rollup_results(
    sample: str, output_file: str,
    num_obs_tags: list[int],
    num_uniq_tags: list[int],
    num_collapsed_tags: list[int],
) -> None:
    # Graph for each region, total tags, unique tags, collapsed tags
    regcount = len(num_obs_tags)
    assert regcount == len(num_uniq_tags) == len(num_collapsed_tags)

    fig = plt.figure(figsize=(50,10))
    fig.suptitle(sample+"\n"+
                "Results of Tag Rollup - 1 error allowed", fontsize=16)
    x = list(range(regcount))
    ax = plt.subplot(111)
    ax.bar(
        x, num_obs_tags, WIDTH,
        alpha=0.5, color='purple', label='Total')
    ax.bar(
        [p + WIDTH for p in x], num_uniq_tags, WIDTH,
        alpha=0.5, color='green', label='Unique')
    ax.bar(
        [p + WIDTH*2 for p in x], num_collapsed_tags, WIDTH,
        alpha=0.5, color='blue', label='Collapsed')
    ax.legend(['Total','Unique','Collapsed'], loc='upper left')
    ax.ticklabel_format(style='plain')
    ax.get_yaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    plt.xticks([p + WIDTH for p in x],range(regcount))
    plt.xlim([min(x)-WIDTH,max(x)+WIDTH*4])
    plt.xlabel("Region",fontsize=14)
    plt.ylabel("Tag Counts",fontsize=14)
    plt.savefig(output_file, dpi=200, facecolor='w', edgecolor='w',
                papertype=None, format=None,
                transparent=False)
    plt.close()
