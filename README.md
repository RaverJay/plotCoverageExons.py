# plotCoverageExons.py
Plot the coverage at some specific region from multiple bam files:

    ./plotCoverageExons.py [-p <geneposition>] [-o <outfile>.pdf] [-options] <bamfile1> [<bamfile2> ...]

This grew over time, excuse the bad code

- Uses pysam to get pile up reads covering the selected regions
- Plots the coverage using matplotlib
- Allows grouping of samples and other customization

### Example

TODO

### Full usage

    ./plotCoverageExons.py [-p <geneposition>] [-options] <bamfile1> [<bamfile2> ...]

    -p <geneposition>          specify genomic position with: <chromosome>:<start>-<end>
    -z <sizefactor1>,<sf2>,... specify normalization size factors for each bamfile in order
    -s <+|->                   specify strand ( + or - ) to be included in plot title
    -e <start>-<end>,<s-e>,... select feature regions to be plotted exclusively
    -l <label1>,<l2>,...       specify labels for the regions
    -c <condition1>,<c2>,...   specify patterns to group bamfiles by filename
    -n <name1>,<n2>,...        specify names for the condition groups
    -u <uname>                 include unmatched bamfiles under name <uname>
    -a                         average samples of the same condition, default off (plot all samples)
    -L                         set logarithmic scale for the y-axis
    -S                         add smoothed curve (rolling mean of 49-mers, twice applied)
    -f                         enable samtools-like filtering of reads (default no filtering)
    -o <file>                  specify output filename and thus output format
    -P <file>                  only write coverage data to pickle file
