import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Load the Fig.4g image
img = plt.imread('data/fig04-g.png')

fig, ax = plt.subplots(figsize=(8, 6))
ax.imshow(img)
ax.axis('off')
ax.set_title('Click on each cholera data point (2000–2016)\nClose the window when done')

# List to store clicked pixel coordinates
clicked_pixels = []

def onclick(event):
    if event.xdata is not None and event.ydata is not None:
        clicked_pixels.append((event.xdata, event.ydata))
        ax.plot(event.xdata, event.ydata, 'ro')
        fig.canvas.draw()

# Connect the click event
cid = fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()

# For demonstration, assume you clicked points in order from 2000 to 2016
# You will need to manually input the pixel->data transformation parameters:
# Enter the pixel coordinates and actual data values for two reference points:
# Example reference: bottom-left corner of data area corresponds to (year=2000, value=-6000)
#                    top-right  corner corresponds to (year=2016, value=8000)
# Replace these with your measured values:
px_min, py_min, year_min, val_min = 107, 306, 2000, -6000
px_max, py_max, year_max, val_max = 336, 33,  2017,  8000

# Transform pixel coords to data coords
years = np.linspace(year_min, year_max, len(clicked_pixels))
detrended = []
for (px, py), year in zip(clicked_pixels, years):
    # Linear mapping from pixel to data
    year_val = ((px - px_min) / (px_max - px_min)) * (year_max - year_min) + year_min
    val = ((py - py_min) / (py_max - py_min)) * (val_max - val_min) + val_min
    detrended.append(val)

# Build DataFrame
df = pd.DataFrame({
    'year': np.round(years).astype(int),
    'detrended_cholera_cases': detrended
})

# Show the extracted table
# import ace_tools as tools; tools.display_dataframe_to_user(name="Tanzania Cholera (Detrended) Time Series", dataframe=df)
print(df)

# Save to CSV for further analysis
df.to_csv('data/tanzania_cholera_detrended.csv', index=False)