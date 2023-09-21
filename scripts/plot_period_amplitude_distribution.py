import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

"""
    This script plots the joint distribution of period and amplitude for a given dataset.

    The dataset is assumed to be a CSV file with at least the following columns:
        - Frequency
        - Amplitude

    The script will save the plot as a PNG file in the given output directory.

Author: Franco Tavella
Date: 09/19/2023
"""

# Read data_file and output_dir as script arguments using sys
if len(sys.argv) != 3:
    raise Exception("Usage: python plot_period_amplitude_distribution.py <data_file> <output_file_with_path>")
else:
    data_file = sys.argv[1]
    output_file_with_path = sys.argv[2]

# Read in data
df = pd.read_csv(data_file)
df['Frequency'] = np.log10(df['Frequency'])
df['Average Amplitude'] = df['Amplitude']

# Plot
sns.jointplot(data=df, x="Frequency", y="Average Amplitude", kind="hex", xlim=(-2, 2), ylim=(-0.01, 1.0),
              gridsize=50, space=0,)

# Set x-axis tick labels as powers of 10
plt.xticks([-2, -1, 0, 1, 2], [r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$', r'$10^{2}$'])

plt.savefig(output_file_with_path, dpi=300, bbox_inches='tight')