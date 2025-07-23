import math
from PIL import Image, ImageDraw, ImageFont
import os
import argparse # Import the argparse module
import numpy as np
import sys
import matplotlib.pyplot as plt

def dag_fill_screen_params(params):
    """
    Given a dict with some of the following keys:
      width_px, width_cm, width_deg,
      height_px, height_cm, height_deg,
      pix_per_cm, pix_per_deg,
      screen_distance_cm,
      aspect_ratio (w/h),
    fills in as many missing ones as possible and returns a new dict.

    All angles are in degrees.
    """
    # copy input to avoid side-effects
    p = params.copy()

    # helper to compute width_cm from width_deg
    def deg_to_cm(deg, dist):
        return 2 * dist * np.tan(np.deg2rad(deg) / 2)

    # cm to deg
    def cm_to_deg(cm, dist):
        return np.rad2deg(2 * np.arctan(cm / (2 * dist)))
    full_set = set([
        'width_px', 'height_px', 'width_deg', 'height_deg', 'width_px', 'height_px', 
        'aspect_ratio', 'pix_per_cm', 'pix_per_deg', 'screen_distance_cm'
    ])
    unit_list = ['cm', 'px', 'deg']
    ratios = set(['aspect_ratio', 'pix_per_cm', 'pix_per_deg', 'deg_per_cm'])
    # iterative inference
    loop_count = 0
    while not set(p.keys()).issubset(ratios):

        # aspect_ratio
        if 'aspect_ratio' not in p:
            for unit in unit_list:
                if (f'width_{unit}' in p) and (f'height_{unit}' in p):
                    p['aspect_ratio'] = p[f'width_{unit}']/p[f'height_{unit}']
                    break
        # px/cm conversion
        if 'pix_per_cm' not in p:
            for merid in ['width', 'height']:
                if (f'{merid}_px' in p) and (f'{merid}_cm' in p):
                    p['pix_per_cm'] = p[f'{merid}_px'] / p[f'{merid}_cm']
        # pix/deg conversion
        if 'pix_per_deg' not in p:
            for merid in ['width', 'height']:
                if (f'{merid}_px' in p) and (f'{merid}_deg' in p):
                    p['pix_per_deg'] = p[f'{merid}_px'] / p[f'{merid}_deg']
                    break

        # pix/deg conversion
        if 'deg_per_cm' not in p:
            for merid in ['width', 'height']:
                if (f'{merid}_deg' in p) and (f'{merid}_cm' in p):
                    p['deg_per_cm'] = p[f'{merid}_deg'] / p[f'{merid}_cm']
                    break
        loop_count +=1
        if loop_count>20:
            break
    
    loop_count = 0
    while set(p.keys()) != full_set:
        # have width have aspect ratio missing height
        if 'aspect_ratio' in p:
            for unit in unit_list:
                if (f'width_{unit}' in p) and (f'height_{unit}' not in p):
                    p[f'height_{unit}'] = p[f'width_{unit}'] / p['aspect_ratio']

        # have height have aspect ratio missing width
        if 'aspect_ratio' in p:
            for unit in unit_list:
                if (f'width_{unit}' not in p) and (f'height_{unit}' in p):
                    p[f'width_{unit}'] = p[f'height_{unit}'] * p['aspect_ratio']


        # height_cm from width_cm & aspect
        if 'height_cm' not in p and 'width_cm' in p and 'aspect_ratio' in p:
            p['height_cm'] = p['width_cm'] / p['aspect_ratio']
        if 'width_cm' not in p and 'height_cm' in p and 'aspect_ratio' in p:
            p['width_cm'] = p['height_cm'] * p['aspect_ratio']

            
        if 'width_px' not in p and 'pix_per_cm' in p and 'width_cm' in p:
            p['width_px'] = int(round(p['pix_per_cm'] * p['width_cm']))
        if 'height_px' not in p and 'pix_per_cm' in p and 'height_cm' in p:
            p['height_px'] = int(round(p['pix_per_cm'] * p['height_cm']))

        # px/deg conversion

        if 'width_px' not in p and 'pix_per_deg' in p and 'width_deg' in p:
            p['width_px'] = int(round(p['pix_per_deg'] * p['width_deg']))

        # deg<->cm using distance
        if 'distance_cm' in p:
            d = p['distance_cm']
            # width_deg from width_cm
            if 'width_deg' not in p and 'width_cm' in p:
                p['width_deg'] = cm_to_deg(p['width_cm'], d)
            # width_cm from width_deg
            if 'width_cm' not in p and 'width_deg' in p:
                p['width_cm'] = deg_to_cm(p['width_deg'], d)
            # height_deg
            if 'height_deg' not in p and 'height_cm' in p:
                p['height_deg'] = cm_to_deg(p['height_cm'], d)
            if 'height_cm' not in p and 'height_deg' in p:
                p['height_cm'] = deg_to_cm(p['height_deg'], d)

        # pix_per_deg from pix_per_cm
        if 'pix_per_deg' not in p and 'pix_per_cm' in p and 'distance_cm' in p:
            # derive via cm/deg
            if 'width_deg' in p:
                # prefer width
                p['pix_per_deg'] = p['pix_per_cm'] * (p['width_cm'] / p['width_deg'])
            else:
                # compute average cm_per_deg
                cm_per_deg = deg_to_cm(1, p['distance_cm'])
                p['pix_per_deg'] = p['pix_per_cm'] * cm_per_deg

        # 
        loop_count +=1
        if loop_count > 10:
            print(p)
            return p
            # bloop
    return p

import matplotlib.pyplot as plt
import os

def create_plot_png(aspect_ratio: float, horizontal_deg: float, vertical_deg: float, output_filename: str = 'output_plot.png'):
    """
    Generates a PNG image using Matplotlib with specified aspect ratio,
    horizontal and vertical degrees, and lines through the center.

    Args:
        aspect_ratio (float): The desired aspect ratio of the figure (width / height).
        horizontal_deg (float): The total span of the x-axis in degrees.
        vertical_deg (float): The total span of the y-axis in degrees.
        output_filename (str): The name of the output PNG file.
    """
    # Define a base width for the figure in inches.
    # The height will be calculated based on the aspect ratio.
    base_width_inches = 10
    figure_height_inches = base_width_inches / aspect_ratio

    # Create a figure and an axes object with the specified figure size.
    # The figsize argument takes (width, height) in inches.
    fig, ax = plt.subplots(figsize=(base_width_inches, figure_height_inches))
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    
    ax.yaxis.set_ticks_position('left')

    # Calculate the x and y limits based on the total horizontal and vertical degrees.
    # The plot will be centered around (0,0).
    x_min = -horizontal_deg / 2
    x_max = horizontal_deg / 2
    y_min = -vertical_deg / 2
    y_max = vertical_deg / 2

    # Set the limits for the x and y axes.
    ax.set_xticks(np.arange(int(x_min), int(x_max)+1))    
    ax.set_yticks(np.arange(int(y_min), int(y_max)+1))    
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Draw a horizontal line through the middle (y=0).
    # ax.axhline(0, color='red', linestyle='--', linewidth=1.5, label='Horizontal Center Line')

    # Draw a vertical line through the middle (x=0).
    # ax.axvline(0, color='blue', linestyle='--', linewidth=1.5, label='Vertical Center Line')

    # Add a title and labels for clarity.
    # ax.set_title(f'Plot with Aspect Ratio {aspect_ratio:.2f} and Extents')
    # ax.set_xlabel(f'Horizontal Degrees (Span: {horizontal_deg:.1f})')
    # ax.set_ylabel(f'Vertical Degrees (Span: {vertical_deg:.1f})')

    # Add a grid for better visualization.
    ax.grid(True, linestyle=':', alpha=0.7)

    # Add a legend to explain the lines.
    # ax.legend()

    # Ensure tight layout to prevent labels/titles from overlapping.

    # Save the figure as a PNG file.
    plt.savefig(output_filename, dpi=300, bbox_inches='tight', pad_inches=0) # dpi (dots per inch) controls image resolution
    print(f"Plot saved successfully as '{output_filename}'")

    # Close the plot to free up memory.
    plt.close(fig)


if __name__ == "__main__":
    argv = sys.argv[1:]
    p_kwargs = {}
    for i,arg in enumerate(argv):        
        if '--' in arg:
            this_param = arg.replace('--', '')
            p_kwargs[this_param] = float(argv[i+1])
    params = dag_fill_screen_params(p_kwargs)
    print(params)

    # create_plot_png(
    #     aspect_ratio=params['aspect_ratio'],
    #     horizontal_deg=params['width_deg'],
    #     vertical_deg=params['height_deg'],
    # )


# {'height_px': 1080.0, 'width_px': 1920.0, 'distance_cm': 220.0, 'width_cm': 69.0, 'aspect_ratio': 1.7777777777777777, 'pix_per_cm': 27.82608695652174, 'height_cm': 38.8125, 'width_deg': 17.82486993904548, 'height_deg': 10.082051873839621, 'pix_per_deg': 107.7146709381721}