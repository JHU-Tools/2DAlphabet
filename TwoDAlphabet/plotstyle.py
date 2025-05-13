'''Container for matplotlib + mplhep format styles'''
# Options for the figure object (analagous to TCanvas)
fig_style = {
    'figsize': (10, 10),
    'dpi': 100
}
ratio_fig_style = {
    'figsize': (10, 10),
    'gridspec_kw': {'height_ratios': (3, 1)},
    'dpi': 100
}

# Options for plotting
stack_style = {
    'edgecolor': (0, 0, 0, 0.5),
}
hatch_style = {
    'facecolor': 'none',
    'edgecolor': (0, 0, 0, 0.5),
    'linewidth': 0,
    'hatch': '///',
}
errorbar_style = {
    'linestyle': 'none',
    'marker': '.',      # display a dot for the datapoint
    'elinewidth': 2,    # width of the errorbar line
    'markersize': 20,   # size of the error marker
    'capsize': 0,       # size of the caps on the errorbar (0: no cap fr)
    'color': 'k',       # black 
}
shaded_style = {
    'facecolor': (0,0,0,0.3),
    'linewidth': 0
}

# Conversion from matplotlib color naming to ROOT TColor naming
mpl_to_root_colors = {
    'white': 0,
    'black': 1,
    'gray': 920, 
    'red': 632, 
    'green': 416, 
    'blue': 600,
    'yellow': 400,
    'magenta': 616,
    'cyan': 432,
    'orange': 800, 
    'spring': 820,
    'teal': 840, 
    'azure': 860,
    'violet': 880,
    'pink': 900
}
# And a helper function to return the matplotlib color name from ROOT TColor
def root_to_matplotlib_color(TColor):
    return list(mpl_to_root_colors.keys())[list(mpl_to_root_colors.values()).index(TColor)]
