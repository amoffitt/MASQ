# QC Metric Plots
# Inputs - all final reports
# Output - plots of each sample, by locus, for a set of metrics
# Some automatic QC rules are applied, and QC failed loci are annotated

# Libraries
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
library(scales)

# Collect list of report files to plot 
args = commandArgs(trailingOnly=TRUE)
outplot = args[1]
outtxt = args[2]
badlocifile = args[3]
protocol = args[4]
infiles = args[-c(1, 2, 3, 4)]

# Number of files
N = length(args) - 4

# Load files and get sample name
load_table = function(tablefile) {
    D=read.table(tablefile, sep="\t", header=T, stringsAsFactors=F)
    D$gt2templates=apply(D[, c("A2", "C2", "G2", "T2")], 1, sum)
    return(D)
}
get_samplename = function(tablefile) {
    Dname=sub('\\.final_report.txt$', '', basename(tablefile))
    return(Dname)
}

# Function to get plots for 1 file
qc_plots_1file = function(D, Dname, qc_fail) {
    qcindex = which(D$loc %in% qc_fail)
    # Number of Tags
    p1 = ggplot(data=D) +
    geom_bar(aes(x=factor(loc), y=Number.of.Tags + 1), stat="identity", fill="lightpink") +
    xlab("Locus") + 
    ylab("Total Templates") + 
    scale_y_continuous(trans="log10", breaks=10^(0:7), labels=math_format(expr=10^.x)(0:7)) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.title = element_text(size=10)) +
    ggtitle(paste0(Dname, "\n", sprintf("%0.2f", mean(D$Number.of.Tags)))) +
    #       "\t",
    #       sprintf("%0.2f",mean(D$Yield:.Rolled-Up.Tags))))+ 
    geom_hline(yintercept=0.10 * median(D$Number.of.Tags), lty=3, col='black')
    # Primer Alignment Rate
    p4 = ggplot(data=D) + 
    geom_bar(aes(x=factor(loc), y=Percent.of.Assigned.that.Pass), stat="identity", fill="orange") +
    xlab("Locus") + 
    ylab("Primer Alignment Rate") + 
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.title = element_text(size=10)) +
    ggtitle(paste0(Dname, "\n", sprintf("%0.2f", mean(D$Percent.of.Assigned.that.Pass)))) +
    ylim(c(0, 1.05)) +
    geom_hline(yintercept=0.75, lty=3, col='red')
    # Uncorrected Variant AF
    p6 = ggplot(data=D) + 
    geom_bar(aes(x=factor(loc), y=VarAF), stat="identity", fill="darkorchid3") + 
    xlab("Locus") + 
    ylab("Uncorrected Variant AF") + 
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=6),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        plot.title = element_text(size=10)) +
    ggtitle(paste0(Dname, "\n", sprintf("%0.2f", mean(D$VarAF))))

    if (protocol != "standard PCR") {
        p2 = ggplot(data=D)+
        geom_bar(aes(x=factor(loc), y=gt2templates), stat="identity", fill="gold2") +
        xlab("Locus") + 
        ylab("Aligned >=2 RPT Templates") + 
        scale_y_continuous(trans = "log10", breaks=10^(0:7), labels=math_format(expr=10^.x)(0:7)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=8),
            axis.text.y = element_text(size=8),
            axis.title.x = element_text(size=10),
            axis.title.y = element_text(size=10),
            plot.title = element_text(size=10)) +
        ggtitle(paste0(Dname, "\n", sprintf("%0.2f", mean(D$gt2templates)))) +
        geom_hline(yintercept=0.10*median(D$gt2templates), lty=3, col='red')


        # Reads Per Template
        p3 = ggplot(data=D)+ 
        geom_bar(aes(x=factor(loc), y=Avg.Reads.Per.Rolled.Up.Tag), stat="identity", fill="green") +
        xlab("Locus") + 
        ylab("Average Reads Per Tags") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=8),
            axis.text.y = element_text(size=8),
            axis.title.x = element_text(size=10),
            axis.title.y = element_text(size=10),
            plot.title = element_text(size=10)) +
        ggtitle(paste0(Dname, "\n", sprintf("%0.2f", mean(D$Avg.Reads.Per.Rolled.Up.Tag)))) +
        geom_hline(yintercept=5, lty=3, col='black')

        # Sequence Alignment Rate
        p5 = ggplot(data=D)+ 
        geom_bar(aes(x=factor(loc), y=Fraction.Aligned.Reads), stat="identity", fill="lightblue") +
        xlab("Locus") + 
        ylab("Sequence Alignment Rate") + 
        theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=6),
            axis.text.y = element_text(size=8),
            axis.title.x = element_text(size=10),
            axis.title.y = element_text(size=10),
            plot.title = element_text(size=10)) +
        ggtitle(paste0(Dname, "\n", sprintf("%0.2f", mean(D$Fraction.Aligned.Reads)))) +
        ylim(c(0, 1.05))+
        geom_hline(yintercept=0.75, lty=3, col='red')
    
        # IF THERE ARE QC FAIL LOCI, ANNOTATE THEM
        if (length(qcindex)>0) {
            p1 = p1 + annotate('point', x=qcindex, y=1.05*max(D$Number.of.Tags), size=5, color='red', shape='*')
            p2 = p2 + annotate("point", x=qcindex, y = 1.05*max(D$gt2templates), size=5, color='red', shape='*')
            p3 = p3 + annotate("point", x=qcindex, y = 1.05*max(D$Avg.Reads.Per.Rolled.Up.Tag), size=5, color='red', shape='*')
            p4 = p4 + annotate("point", x=qcindex, y = 1.02, size=5, color='red', shape='*')
            p5 = p5 + annotate("point", x=qcindex, y = 1.02, size=5, color='red', shape='*')
            p6 = p6 + annotate("point", x=qcindex, y = 1.05*max(D$VarAF), size=5, color='red', shape='*')
    
        }
        onesamp_plots = list(p1, p2, p3, p4, p5, p6)
        return(onesamp_plots)
    } else {
        if (length(qcindex) > 0) {
            p1 = p1 + annotate('point', x=qcindex, y=1.05 * max(D$Number.of.Tags), size=5, color='red', shape='*')
            p4 = p4 + annotate("point", x=qcindex, y=1.02, size=5, color='red', shape='*')
            p6 = p6 + annotate("point", x=qcindex, y=1.05 * max(D$VarAF), size=5, color='red', shape='*')
        }
        onesamp_plots = list(p1, p4, p6)
        return(onesamp_plots)
    }
}
  
# QC filters...
qc_filters = function(alltables) {
    locusnames = alltables[[1]]$loc
  
    get_primer_rate = function(D) {
        return(D$Percent.of.Assigned.that.Pass)
    }
    combinedcols = do.call(cbind, lapply(alltables, get_primer_rate))
    max_primer_rate = apply(combinedcols, 1, max)
  
    get_align_rate = function(D) {
        return(D$Fraction.Aligned.Reads)
    }
    combinedcols = do.call(cbind, lapply(alltables, get_align_rate))
    max_align_rate = apply(combinedcols, 1, max)
    min_align_rate = apply(combinedcols, 1, min)
  
    get_total_templates = function(D) {
        return(D$gt2templates)
    }
    combinedcols = do.call(cbind, lapply(alltables, get_total_templates))
    medcounts = apply(combinedcols, 2, median)
    medcountmat = matrix(medcounts, nrow=nrow(combinedcols), ncol=length(medcounts), byrow=TRUE)
    n_okcount = apply(combinedcols > (medcountmat * 0.1), 1, sum)
    
    # return(locusnames[max_primer_rate<0.75 | max_align_rate<0.75 ])
    return(locusnames[max_primer_rate < 0.75 | min_align_rate < 0.5 | max_align_rate < 0.75 | n_okcount == 0])
}

# Load all files
alltables = lapply(infiles, load_table)
allnames = unlist(lapply(infiles, get_samplename))

# QC fail loci
qc_fail = qc_filters(alltables)

# Add bad loci
bl = read.table(badlocifile, header=F, col.names=c("Locus"))
if (nrow(bl) > 0) {
    qc_fail = unique(c(qc_fail, bl$Locus))
}

# Write to file
write.table(qc_fail, file=outtxt, quote=F, sep="\t", col.names=F, row.names=F)

# Plot all samples, with annotated QC fails
allplots = mapply(qc_plots_1file, alltables, allnames, MoreArgs=list(qc_fail=qc_fail))

# Plot grid
# Rows - samples
# Columns - metrics
grDevices::pdf(NULL)
if (protocol == "standard PCR") {
    p = plot_grid(plotlist=allplots, nrow=N, ncol=3)
    # generates Rplot.pdf, if null device is not opened (in some version of ggplot / R / cowplot)
    save_plot(outplot, p, base_height=3*N, base_width=18)
} else {
    p=plot_grid(plotlist=allplots, nrow=N, ncol=6) 
    # generates Rplot.pdf, if null device is not opened (in some version of ggplot / R / cowplot)
    save_plot(outplot, p, base_height = 3*N, base_width = 24)
}
grDevices::dev.off()
