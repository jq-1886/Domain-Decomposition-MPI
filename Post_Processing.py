import os
import glob

import numpy as np
import matplotlib.pyplot as plt

from PIL import Image
from natsort import natsorted

# find the diretory with the most recent mtime
latest_subdirectory = max(glob.glob(os.path.join('./output/', '*/')),
                          key=os.path.getmtime)

# creating a list of all folders within latest_subdirectory
all_iteration_directories = os.listdir(latest_subdirectory)

# creating a the path to the images folder
image_path = latest_subdirectory + 'images'
# creating the image directory
os.mkdir(image_path)

# creating a loop to move over each folder in latest_subdirectory
# and merge .dat files into one
for iteration_directory in all_iteration_directories:
    it_path = latest_subdirectory + iteration_directory
    # using natsorted to sort it_file_list naturally
    it_file_list = natsorted(os.listdir(it_path))
    merged_file = open(it_path + '/' + iteration_directory +
                       '_merged.dat', 'w')
    for it_file in it_file_list:
        # print all entries that are files
        it_file_path = it_path + '/' + it_file
        with open(it_file_path) as read_file:
            merged_file.write(read_file.read())
        read_file.close()
    merged_file.close()

    # converting directory name into a string
    date_time_string = latest_subdirectory[9:11] + '/' + latest_subdirectory[
        11:13] + '/' + latest_subdirectory[
        13:15] + ' at ' + latest_subdirectory[
        16:18] + ':' + latest_subdirectory[
        18:20] + ':' + latest_subdirectory[20:22]

    # setting up the plot
    plot_data = np.loadtxt(it_path + '/' + iteration_directory + '_merged.dat')
    plot_rows, plot_cols = plot_data.shape
    x = np.linspace(0, 1, plot_rows)
    y = np.linspace(0, 1, plot_cols)

    # setting up the meshgrid and plot data
    X, Y = np.meshgrid(x, y)
    Z = plot_data

    # populating fig and ax
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')
    cont = ax.contourf(X, Y, Z, cmap='Blues')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.text(0, -0.1, 'Iteration ' + iteration_directory + ' performed on a ' +
            str(plot_rows) + 'x' + str(plot_cols) + ' grid on ' +
            date_time_string, fontsize=12, verticalalignment='bottom',
            horizontalalignment='left', color='black')
    # creating the colorbar key
    fig.colorbar(cont)

    # populating tiles and saving file
    ax.set_title('Iteration ' + iteration_directory, fontsize=18)
    fig.savefig(image_path + '/' + 'Iteration_' + iteration_directory +
                '_plot.png')
    plt.close()

# creating a list to store the .png images before gif conversion
gif_frames = []
# populating gif_image_list with natural sorting
gif_image_list = natsorted(glob.glob(image_path + '/' + '*.png'))

# creating a string to name the gif
gif_date_string = latest_subdirectory[9:22]

# opening the images and appending to gif_frames list
for image in gif_image_list:
    new_gif_frame = Image.open(image)
    gif_frames.append(new_gif_frame)

# saving gif
gif_frames[0].save(image_path + '/' + 'Simulation_Animation_' +
                   gif_date_string + '.gif', format='GIF',
                   append_images=gif_frames[1:], save_all=True, duration=300,
                   loop=0)
