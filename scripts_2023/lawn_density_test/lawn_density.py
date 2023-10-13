import os
import cv2
import numpy as np
import matplotlib.pyplot as plt
import re
import csv
import argparse


# Function Definitions
def load_images(directory):
    """Load and sort images by natural order."""
    def atoi(text):
        return int(text) if text.isdigit() else text
    
    def natural_keys(text):
        return [atoi(c) for c in re.split(r'(\d+)', text)]

    return [cv2.imread(os.path.join(directory, filename)) for filename in sorted(os.listdir(directory), key=natural_keys) if filename.endswith(".tif")]


def preprocess_image(image):
    """Convert to grayscale, blur, and threshold."""
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    blurred = cv2.GaussianBlur(gray, (5, 5), 0)
    adaptive_thresh = cv2.adaptiveThreshold(blurred, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 13, 2)
    return cv2.bitwise_not(adaptive_thresh)


def find_petri_dishes(image, min_radius, max_radius):
    """Identify petri dishes using HoughCircles."""
    circles = cv2.HoughCircles(image, cv2.HOUGH_GRADIENT, dp=1, minDist=50, param1=100, param2=30, minRadius=min_radius, maxRadius=max_radius)
    return [(x, y, r) for x, y, r in np.round(circles[0, :]).astype("int")] if circles is not None else []


def measure_density_and_variation(cropped_binary_image, cropped_gray_image):
    """Calculate the density and variation of the bacterial lawn."""
    return np.mean(cropped_gray_image), np.std(cropped_gray_image)


def draw_petri_dishes(image, petri_dishes):
    """Draw circles around petri dishes and label them."""
    output = image.copy()
    for index, (x, y, r) in enumerate(petri_dishes):
        cv2.circle(output, (x, y), r, (0, 255, 0), 2)
        cv2.putText(output, f'Plate {index + 1}', (x - r, y - r), cv2.FONT_HERSHEY_SIMPLEX, 3, (255, 255, 255), 2)
    return output


def export_to_csv(filename, header, data, endpoint_dir):
    """Export data to a CSV file."""
    with open(os.path.join(endpoint_dir, filename), 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(header)
        csv_writer.writerows(data)



def main(directory, endpoint_dir, min_radius, max_radius):
    if not os.path.exists(endpoint_dir):
        os.makedirs(endpoint_dir)

    images = load_images(directory)
    first_image = images[0]
    petri_dishes = find_petri_dishes(preprocess_image(first_image), min_radius, max_radius)

    image_with_rectangles = draw_petri_dishes(first_image, petri_dishes)
    plt.imshow(cv2.cvtColor(image_with_rectangles, cv2.COLOR_BGR2RGB))


    densities = []
    variations = []

    for image in images:
        binary_image = preprocess_image(image)
        gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        dish_densities = []
        dish_variations = []

        # Measure densities and variations for each petri dish
        for x, y, r in petri_dishes:
            cropped_gray_image = gray_image[y-r:y+r, x-r:x+r]
            cropped_binary_image = binary_image[y-r:y+r, x-r:x+r]
            density, variation = measure_density_and_variation(cropped_binary_image, cropped_gray_image)
            dish_densities.append(density)
            dish_variations.append(variation)

        densities.append(dish_densities)
        variations.append(dish_variations)

    # Visualize the density results
    time_points = range(len(densities))
    for dish_index in range(len(petri_dishes)):
        dish_densities = [density[dish_index] for density in densities]
        plt.plot(time_points, dish_densities, label=f'Dish {dish_index + 1}')

    plt.xlabel('Time (hours)')
    plt.ylabel('Bacterial Lawn Density')
    plt.legend()


    # Visualize the variation results
    time_points = range(len(variations))
    for dish_index in range(len(petri_dishes)):
        dish_variations = [variation[dish_index] for variation in variations]
        plt.plot(time_points, dish_variations, label=f'Dish {dish_index + 1}')

    plt.xlabel('Time (hours)')
    plt.ylabel('Bacterial Lawn Variation')
    plt.legend()


    cv2.imwrite(os.path.join(endpoint_dir, 'image_with_rectangles.tif'), image_with_rectangles)


    # Export density data to CSV
    header = ['Plate'] + [ i + 1 for i in range(len(time_points))]
    density_data = [[f'Dish {i + 1}', *row] for i, row in enumerate(zip(*densities))]
    variation_data = [[f'Dish {i + 1}', *row] for i, row in enumerate(zip(*variations))]

    export_to_csv('density_data.csv', header, density_data, endpoint_dir)
    export_to_csv('variation_data.csv', header, variation_data, endpoint_dir)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process petri dish images.')
    parser.add_argument('directory', type=str, help='Path to folder with images')
    parser.add_argument('endpoint_dir', type=str, help='Directory to save output files')
    parser.add_argument('min_radius', type=int, help='Minimum radius of petri dishes')
    parser.add_argument('max_radius', type=int, help='Maximum radius of petri dishes')

    args = parser.parse_args()
    main(args.directory, args.endpoint_dir, args.min_radius, args.max_radius)