from __future__ import print_function
import json
import numpy as np
import matplotlib.pyplot as plt

# Load previously computed results from the JSON file
filename = 'grb_data.json'
print('Reading', filename)
data = json.load(open(filename))
grb_info = data['grb_info']
sigma = np.array(data['sigma'])
sigma = np.nan_to_num(sigma)
sigma = np.clip(sigma, 0., 1.e8)
t_start = data['t_start']
t_stop = data['t_stop']

# Plot the significance for each time interval
plt.figure()
x, y = np.meshgrid(t_start, t_stop)
plt.pcolormesh(x, y, sigma, cmap='hot')
plt.ylim(t_start[-1], t_start[0])

plt.xlabel('Observation stop time (s)')
plt.ylabel('Observation start time (s)')
plt.xlim(20, 172800)
plt.xscale('log')
plt.yscale('log')
plt.colorbar(label='Detection significance')
plt.title('{} (z={})'.format(grb_info['name'], grb_info['redshift']))
plt.tight_layout()

filename = 'grb_twindow.png'
print('Writing', filename)
plt.savefig(filename, format='png', dpi=300)

plt.show()
