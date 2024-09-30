#include <iostream>
#include <random>
#include <cassert>

#include <Eigen/Dense>

// from https://github.com/nothings/stb/tree/master
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <image_path>" << std::endl;
		return 1;
	}

	const char* input_image_path = argv[1];

	// Load the image using stb_image
	int width, height, channels;
	unsigned char* original_image = stbi_load(
		input_image_path, &width, &height, &channels, 3);  // Force load as RGB
	if (!original_image) {
		std::cerr << "Error: Could not load image " << input_image_path
				  << std::endl;
		return 1;
	}

	// Task 1: report size of image
	std::cout << "Image loaded: " << width << "x" << height << " with "
			  << channels << " channels." << std::endl;

	// Task 2: add a noise in the range [-50;50] to each pixel
	Eigen::MatrixXd noisy(height, width);

	std::random_device dev;
	std::mt19937 rnd{dev()};
	std::uniform_real_distribution<double> dist{-50.0, 50.0};

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			const int index = (i * width + j) * 3;	// 3 channels (RGB)
			const double original = static_cast<double>(original_image[index]);
			noisy(i, j) = static_cast<unsigned char>(
				std::clamp(original + dist(rnd), 0.0, 255.0));
		}
	}

	const std::string noisy_image_path = "noisy.png";
	if (stbi_write_png(noisy_image_path.c_str(), width, height, 1, noisy.data(),
					   width) == 0) {
		std::cerr << "Error: Could not save noisy image" << std::endl;

		return 1;
	}
	std::cout << "Noisy image saved to " << noisy_image_path << std::endl;

	// Prepare Eigen matrices for each RGB channel
	Eigen::MatrixXd red(height, width);
	Eigen::MatrixXd green(height, width);
	Eigen::MatrixXd blue(height, width);

	// Fill the matrices with image data
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			const int index = (i * width + j) * 3;	// 3 channels (RGB)
			red(i, j) = static_cast<double>(original_image[index]) / 255.0;
			green(i, j) =
				static_cast<double>(original_image[index + 1]) / 255.0;
			blue(i, j) = static_cast<double>(original_image[index + 2]) / 255.0;
		}
	}

	// Free memory!!!
	stbi_image_free(original_image);

	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic,
				  Eigen::RowMajor>
		grayscale_image(height, width);

	// Save the grayscale image using stb_image_write
	const std::string output_image_path = "output_grayscale.png";
	if (stbi_write_png(output_image_path.c_str(), width, height, 1,
					   grayscale_image.data(), width) == 0) {
		std::cerr << "Error: Could not save grayscale image" << std::endl;

		return 1;
	}
	std::cout << "Grayscale image saved to " << output_image_path << std::endl;

	return 0;
}
