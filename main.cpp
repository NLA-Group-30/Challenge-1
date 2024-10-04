#include <iostream>
#include <random>
#include <cassert>

#include <Eigen/Dense>
#include <Eigen/Sparse>

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
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
		original_matrix(height, width);
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
	Eigen::RowVectorXd v =
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

	// Task 4: Write the convolution operation corresponding to the smoothing
	// kernel Hav2 as a matrix-vector multiplication between a matrix A1 having
	// size (m*n)x(m*n) and the image vector. Report the number of non-zero
	// entries in A1.
	Eigen::SparseMatrix<double, Eigen::RowMajor> A1(v.size(), v.size());
	std::vector<Eigen::Triplet<double>> triplets;
	// we know the exact number of non-zero values we will need
	const size_t expected_entries =
		// all elements that are not on the border
		9 * (width - 2) * (height - 2) +
		// all elements that are on the border but not in a corner
		6 * (2 * width + 2 * height - 8) +
		// all 4 corners
		4 * 4;
	triplets.reserve(expected_entries);
	for (int pixel_index{0}; pixel_index < v.size(); pixel_index++) {
		const int row = pixel_index / width;
		const int col = pixel_index % width;

		const double x = 1.0 / 9.0;

		const int northwest_index = (row - 1) * width + (col - 1);
		const int north_index = (row - 1) * width + col;
		const int northeast_index = (row - 1) * width + (col + 1);

		const int west_index = row * width + (col - 1);
		const int center_index = row * width + col;
		const int east_index = row * width + (col + 1);

		const int southwest_index = (row + 1) * width + (col - 1);
		const int south_index = (row + 1) * width + col;
		const int southeast_index = (row + 1) * width + (col + 1);

		if (row > 0) {
			if (col > 0) {
				triplets.push_back({pixel_index, northwest_index, x});
			}
			triplets.push_back({pixel_index, north_index, x});
			if (col < width - 1) {
				triplets.push_back({pixel_index, northeast_index, x});
			}
		}

		if (col > 0) {
			triplets.push_back({pixel_index, west_index, x});
		}
		triplets.push_back({pixel_index, center_index, x});
		if (col < width - 1) {
			triplets.push_back({pixel_index, east_index, x});
		}

		if (row < height - 1) {
			if (col > 0) {
				triplets.push_back({pixel_index, southwest_index, x});
			}
			triplets.push_back({pixel_index, south_index, x});
			if (col < width - 1) {
				triplets.push_back({pixel_index, southeast_index, x});
			}
		}
	}
	assert(expected_entries == triplets.size());
	A1.setFromTriplets(triplets.begin(), triplets.end());
	std::cout << "Number of non-zero entries: " << A1.nonZeros() << std::endl;

	std::cout << " A1: " << A1.rows() << " x " << A1.cols() << std::endl;
	std::cout << " w : " << w.rows() << " x " << w.cols() << std::endl;

	// Task 5: Apply the previous smoothing filter to the noisy image by
	// performing the matrix-vector multiplication A1*w. Export and upload the
	// resulting image.
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
		smooth_image = A1 * w;
	assert(smooth_image.rows() == w.rows() && smooth_image.cols() == w.cols());
	smooth_image = smooth_image.reshaped(height, width);
	save_image(smooth_image, "smooth.png");

	return 0;
}
