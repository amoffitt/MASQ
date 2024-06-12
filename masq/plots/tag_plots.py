
import numpy as np
import matplotlib.pyplot as plt


def plot_number_of_reads_per_tag(
    numreads: np.ndarray,
    numtags: np.ndarray,
    totalreads: int,
    totaltags: int,
    sample: str = "",
    filename: str = "number_reads_per_tag.png",
    logscale: bool = True,
) -> None:

    fig = plt.figure(figsize=(50,10))
    fig.suptitle(
        f"{sample}\n"
        f"Number of Reads Per Tag\n"
        f"Total Reads: {totalreads}\n"
        f"Total Tags: {totaltags}", fontsize=16)

    plt.bar(numreads,numtags,0.8,log=logscale)

    plt.xlabel("Reads Per Tag",fontsize=14)
    plt.ylabel("Number of Tags",fontsize=14)
    plt.savefig(filename, dpi=200, facecolor='w', edgecolor='w',
                papertype=None, format=None,
                transparent=False)
    plt.close()


def plot_at_least_x_reads_per_tag(
    numreads: np.ndarray,
    numtags: np.ndarray,
    totalreads: int,
    totaltags: int,
    sample: str = "",
    filename: str = "atleastX_reads_per_tag.png",
) -> None:
    fig = plt.figure(figsize=(50,10))
    fig.suptitle(
        f"{sample}\n"
        f"Number of Tags with at Least X Reads\n"
        f"Total Reads: {totalreads}\n"
        f"Total Tags: {totaltags}", fontsize=16)
    x=[]
    y=[]
    for xbin in range(1,100):
        x.append(xbin)
        y.append(sum(numtags[numreads>=xbin]))
    plt.bar(x,y,0.8,color='green')
    plt.xlabel("Reads Per Tag",fontsize=14)
    plt.ylabel("Number of Tags with at least X Reads Per Tag",fontsize=14)
    plt.savefig(
        filename, dpi=200, facecolor='w', edgecolor='w',
        papertype=None, format=None,
        transparent=False)
    plt.close()
