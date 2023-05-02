"""
Stitch together all the pictures and make an animation
"""

import os
import imageio

png_dir = os.path.join(".", "images", "animation")

png_files = sorted([os.path.join(png_dir, f) for f in os.listdir(png_dir) if f.endswith('.png')])

# Create a writer object to write the video
writer = imageio.get_writer('output.mp4', fps=60)

# Loop over the PNG files and add them to the video
print(png_files)
for png_file in png_files:
    img = imageio.imread(png_file)
    writer.append_data(img)

# Close the writer to finalize the video
writer.close()
