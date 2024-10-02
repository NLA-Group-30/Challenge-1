#include <iostream>
#include <random>
#include <cassert>

#include <Eigen/Dense>

// from https://github.com/nothings/stb/tree/master
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

void save_image(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
									Eigen::RowMajor>& m,
				const std::string& file_name) {
	// convert original image in one made of bytes instead of doubles
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic,
				  Eigen::RowMajor>
		tmp(m.rows(), m.cols());
	tmp = m.unaryExpr([](double val) -> unsigned char {
		return static_cast<unsigned char>(val);
	});
	// Save the image
	if (stbi_write_png(file_name.c_str(), tmp.cols(), tmp.rows(), 1, tmp.data(),
					   tmp.cols()) == 0) {
		std::cerr << " \u001b[31mERROR\u001b[0m: Could not save image to "
				  << file_name << std::endl;
		return;
	}
	std::cout << "Image saved to " << file_name << std::endl;
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <image_path>" << std::endl;
		return 1;
	}

	const char* input_image_path = argv[1];

	// Load the image using stb_image
	int width;
	int height;
	int channels;
	unsigned char* original_image =
		stbi_load(input_image_path, &width, &height, &channels, 1);
	if (!original_image) {
		std::cerr << "Error: Could not load image " << input_image_path
				  << std::endl;
		return 1;
	}

	// Task 1: report size of image
	std::cout << "Image loaded: " << width << "x" << height << " with "
			  << channels << " channels." << std::endl;

	// Create a Eigen copy of the original matrix
	Eigen::MatrixXd original_matrix(height, width);
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			const int index = (i * width + j) * channels;
			original_matrix(i, j) = static_cast<double>(original_image[index]);
		}
	}

	// we do not need the original image data
	stbi_image_free(original_image);

	// Task 2: add a noise in the range [-50;50] to each pixel
	Eigen::MatrixXd noisy(height, width);

	std::random_device dev;
	std::mt19937 rnd{dev()};
	std::uniform_real_distribution<double> dist{-50.0, 50.0};

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			const double original = original_matrix(i, j);
			const double x = std::clamp(original + dist(rnd), 0.0, 255.0);
			noisy(i, j) = x;
		}
	}

	save_image(noisy, "noisy.png");

	// Task 3: reshape the original matrix in the vector v and the noisy one in
	// the vector w. Make sure that they both have the shape m*n. Report the
	// euclidean norm of v.
	Eigen::VectorXd v =
		Eigen::Map<Eigen::VectorXd>(original_matrix.data(), height * width);
	Eigen::VectorXd w =
		Eigen::Map<Eigen::VectorXd>(noisy.data(), height * width);
	std::cout << " size of v = " << v.size() << std::endl;
	std::cout << " size of w = " << w.size() << std::endl;
	assert(v.size() == original_matrix.cols() * original_matrix.rows());
	assert(w.size() == original_matrix.cols() * original_matrix.rows());

	double norm{0};
	for (int i = 0; i < v.size(); i++) {
		norm += v(i) * v(i);
	}
	norm = std::sqrt(norm);
	std::cout << " norm of v = " << norm << std::endl;

	return 0;
}
