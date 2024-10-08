#include <cassert>
#include <iostream>
#include <random>
#include <sstream>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

// from https://github.com/nothings/stb/tree/master
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

void save_image(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& m,
				const std::string& file_name) {
	// convert original image in one made of bytes instead of doubles
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmp(m.rows(), m.cols());
	tmp = m.unaryExpr([](const double val) -> unsigned char { return static_cast<unsigned char>(val); });
	// Save the image
	if (stbi_write_png(file_name.c_str(), tmp.cols(), tmp.rows(), 1, tmp.data(), tmp.cols()) == 0) {
		std::cerr << " \u001b[31mERROR\u001b[0m: Could not save image to " << file_name << std::endl;
		return;
	}
	std::cout << "Image saved to " << file_name << std::endl;
}

Eigen::SparseMatrix<double, Eigen::RowMajor> convolution(const Eigen::Matrix3d& convolution_filter, const int width,
														 const int height) {
	Eigen::SparseMatrix<double, Eigen::RowMajor> spMat(width * height, width * height);
	std::vector<Eigen::Triplet<double>> A1_triplets;
	for (int pixel_index{0}; pixel_index < width * height; pixel_index++) {
		const int row = pixel_index / width;
		const int col = pixel_index % width;

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
				if (convolution_filter(0, 0) != 0.0)
					A1_triplets.push_back({pixel_index, northwest_index, convolution_filter(0, 0)});
			}
			if (convolution_filter(0, 1) != 0.0)
				A1_triplets.push_back({pixel_index, north_index, convolution_filter(0, 1)});
			if (col < width - 1) {
				if (convolution_filter(0, 2) != 0.0)
					A1_triplets.push_back({pixel_index, northeast_index, convolution_filter(0, 2)});
			}
		}

		if (col > 0) {
			if (convolution_filter(1, 0) != 0.0)
				A1_triplets.push_back({pixel_index, west_index, convolution_filter(1, 0)});
		}
		if (convolution_filter(1, 1) != 0.0)
			A1_triplets.push_back({pixel_index, center_index, convolution_filter(1, 1)});
		if (col < width - 1) {
			if (convolution_filter(1, 2) != 0.0)
				A1_triplets.push_back({pixel_index, east_index, convolution_filter(1, 2)});
		}

		if (row < height - 1) {
			if (col > 0) {
				if (convolution_filter(2, 0) != 0.0)
					A1_triplets.push_back({pixel_index, southwest_index, convolution_filter(2, 0)});
			}
			if (convolution_filter(2, 1) != 0.0)
				A1_triplets.push_back({pixel_index, south_index, convolution_filter(2, 1)});
			if (col < width - 1) {
				if (convolution_filter(2, 2) != 0.0)
					A1_triplets.push_back({pixel_index, southeast_index, convolution_filter(2, 2)});
			}
		}
	}
	spMat.setFromTriplets(A1_triplets.begin(), A1_triplets.end());
	return spMat;
}

bool is_symmetric(const Eigen::SparseMatrix<double, Eigen::RowMajor>& m) {
	if (m.rows() != m.cols()) {
		// rectangular matrices cannot be symmetric
		return false;
	}

	const Eigen::SparseMatrix<double, Eigen::RowMajor> mT(m.transpose());
	return (m - mT).norm() == 0.0;
}

void save_vector(const char* filename, const Eigen::VectorXd& vector) {
	const long n = vector.size();
	FILE* out = fopen(filename, "w");
	fprintf(out, "%%%%MatrixMarket vector coordinate real general\n");
	fprintf(out, "%ld\n", n);
	for (int i = 0; i < n; i++) {
		fprintf(out, "%d %f\n", i + 1, vector(i));
	}
	fclose(out);
}

Eigen::VectorXd read_vector(const char* filename) {
	std::ifstream input_file(filename);
	if (!input_file.is_open()) {
		std::cerr << "Error opening the file " << filename << "." << std::endl;
		std::exit(-1);
	}

	std::string line;
	long int num_elements = 0;

	// skip the first line as we know it is always a comment
	std::getline(input_file, line);

	if (std::getline(input_file, line)) {
		std::istringstream iss(line);
		if (!(iss >> num_elements)) {
			std::cerr << "Error reading the number of elements." << std::endl;
			std::exit(-1);
		}
		if (num_elements < 1) {
			std::cerr << "Error: wrong number of elements: " << num_elements << std::endl;
			std::exit(-1);
		}
	} else {
		std::cerr << "Error reading from the file " << filename << "." << std::endl;
		std::exit(-1);
	}

	Eigen::VectorXd x(num_elements);

	while (std::getline(input_file, line)) {
		std::istringstream iss(line);
		long int index;
		double value;

		if (!(iss >> index >> value)) {
			std::cerr << "Error reading index and value." << std::endl;
			std::exit(-1);
		}

		// LIS works with 1-based indices
		if (index < 1 || index > num_elements) {
			std::cerr << "Index out of bounds: " << index << std::endl;
			std::exit(-1);
		}

		x(index - 1) = value;
	}
	input_file.close();

	return x;
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
	unsigned char* original_image = stbi_load(input_image_path, &width, &height, &channels, 1);
	if (!original_image) {
		std::cerr << "Error: Could not load image " << input_image_path << std::endl;
		return 1;
	}

	// Task 1: report size of image
	std::cout << "Image loaded: " << width << "x" << height << " with " << channels << " channels." << std::endl;

	// Create a Eigen copy of the original matrix
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> original_matrix(height, width);
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			const int index = (i * width + j) * channels;
			original_matrix(i, j) = static_cast<double>(original_image[index]);
		}
	}

	// we do not need the original image data
	stbi_image_free(original_image);

	// Task 2: add a noise in the range [-50;50] to each pixel
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> noisy(height, width);

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
	Eigen::VectorXd v = Eigen::Map<Eigen::VectorXd>(original_matrix.data(), height * width);
	Eigen::VectorXd w = Eigen::Map<Eigen::VectorXd>(noisy.data(), height * width);
	std::cout << " size of v = " << v.size() << std::endl;
	std::cout << " size of w = " << w.size() << std::endl;
	assert(v.size() == original_matrix.cols() * original_matrix.rows());
	assert(w.size() == original_matrix.cols() * original_matrix.rows());

	std::cout << " norm of v = " << v.norm() << std::endl;

	// Task 4: Write the convolution operation corresponding to the smoothing
	// kernel Hav2 as a matrix-vector multiplication between a matrix A1 having
	// size (m*n)x(m*n) and the image vector. Report the number of non-zero
	// entries in A1.
	Eigen::Matrix3d Hav2;
	const double one_ninth = 1.0 / 9.0;
	Hav2 << one_ninth, one_ninth, one_ninth, one_ninth, one_ninth, one_ninth, one_ninth, one_ninth, one_ninth;
	Eigen::SparseMatrix<double, Eigen::RowMajor> A1 = convolution(Hav2, width, height);
	std::cout << "Number of A1 non-zero entries: " << A1.nonZeros() << std::endl;

	// Task 5: Apply the previous smoothing filter to the noisy image by
	// performing the matrix-vector multiplication A1*w. Export and upload the
	// resulting image.
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> smooth_image = A1 * w;
	assert(smooth_image.rows() == w.rows() && smooth_image.cols() == w.cols());
	smooth_image = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
		smooth_image.data(), height, width);
	save_image(smooth_image, "smooth.png");

	// Task 6: Write the convolution operation corresponding to the sharpening
	// kernel Hsh2 as a matrix-vector multiplication by a matrix A2 having size
	// (m*n)x(m*n). Report the number of non-zero entries in A2. Is A2
	// symmetric?
	Eigen::Matrix3d Hsh2;
	Hsh2 << 0.0, -3.0, 0.0, -1.0, 9.0, -3.0, 0.0, -1.0, 0.0;
	Eigen::SparseMatrix<double, Eigen::RowMajor> A2 = convolution(Hsh2, width, height);
	std::cout << "Number of A2 non-zero entries: " << A2.nonZeros() << std::endl;
	std::cout << "A2 is " << (is_symmetric(A2) ? "" : "not ") << "symmetric" << std::endl;

	// Task 7: Apply the previous sharpening filter to the original image by
	// performing the matrix-vector multiplication A2*v. Export and upload the
	// resulting image.
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> sharp_image = A2 * v;
	assert(sharp_image.rows() == v.rows() && sharp_image.cols() == v.cols());
	sharp_image = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(sharp_image.data(),
																									 height, width);
	sharp_image = sharp_image.unaryExpr([](const double val) -> double { return std::clamp(val, 0.0, 255.0); });
	save_image(sharp_image, "sharp.png");

	// Task 8: Export the Eigen matrix A2 and vector w in the .mtx format. Using
	// a suitable iterative solver and preconditioner technique available in the
	// LIS library compute the approximate solution to the linear system A2*x =
	// w prescribing a tolerance of 1e−9. Report here the iteration count and
	// the final residual.
	Eigen::saveMarket(A2, "A2.mtx");
	std::cout << "Matrix A2 saved to A2.mtx" << std::endl;
	save_vector("w.mtx", w);
	std::cout << "Vector w saved to w.mtx" << std::endl;

	system("./task_8.x A2.mtx w.mtx sol.txt hist.txt");

	// Task 9: Import the previous approximate solution vector x in Eigen and
	// then convert it into a .png image. Upload the resulting file here.
	Eigen::VectorXd x = read_vector("sol.txt");
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> solution =
		Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(x.data(), height, width);
	save_image(solution, "solution.png");

	// Task 10: Write the convolution operation corresponding to the detection
	// kernel Hlap as a matrix-vector multiplication by a matrix A3 having size
	// (m*n)x(m*n). Is matrix A3 symmetric?
	Eigen::Matrix3d Hlap;
	Hlap << 0.0, -1.0, 0.0, -1.0, 4.0, -1.0, 0.0, -1.0, 0.0;
	Eigen::SparseMatrix<double, Eigen::RowMajor> A3 = convolution(Hlap, width, height);
	std::cout << "Number of A3 non-zero entries: " << A3.nonZeros() << std::endl;
	std::cout << "A3 is " << (is_symmetric(A3) ? "" : "not ") << "symmetric" << std::endl;

	// Task 11: Apply the previous edge detection filter to the original image
	// by performing the matrix-vector multiplication A3*v. Export and upload
	// the resulting image.
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> edge_image = A3 * v;
	assert(edge_image.rows() == v.rows() && edge_image.cols() == v.cols());
	edge_image = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(edge_image.data(),
																									height, width);
	edge_image = edge_image.unaryExpr([](const double val) -> double { return std::clamp(val, 0.0, 255.0); });
	save_image(edge_image, "edge.png");

	// Task 12-13: Using a suitable iterative solver available in the Eigen library compute the approximate solution of
	// the linear system (I+A3)y = w, where I denotes the identity matrix, prescribing a tolerance of 10−10. Report here
	// the iteration count and the final residual. Convert the image stored in the vector y into a .png image and upload
	// it.
	Eigen::SparseMatrix<double, Eigen::RowMajor> I(A3.rows(), A3.cols());
	I.setIdentity();
	Eigen::VectorXd y(A3.rows());
	const double tol = 1.0e-10;
	const int maxit = 1000;

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>> cg;
	cg.setMaxIterations(maxit);
	cg.setTolerance(tol);
	cg.compute(I + A3);
	y = cg.solve(w);
	std::cout << "Eigen native CG" << std::endl;
	std::cout << "#iterations: " << cg.iterations() << std::endl;
	std::cout << "relative residual: " << cg.error() << std::endl;

	// Task 13: Convertire l'immagine memorizzata nel vettore y in un'immagine .png e caricarla.
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> output_image =
		Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(y.data(), height, width);
	output_image = output_image.unaryExpr([](const double val) -> double { return std::clamp(val, 0.0, 255.0); });
	save_image(output_image, "y.png");

	return 0;
}
