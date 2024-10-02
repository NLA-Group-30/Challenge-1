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
        input_image_path, &width, &height, &channels, 1);  // Force load as RGB
    if (!original_image) {
        std::cerr << "Error: Could not load image " << input_image_path
                  << std::endl;
        return 1;
    }

    // Task 1: report size of image
    std::cout << "Image loaded: " << width << "x" << height << " with "
              << channels << " channels." << std::endl;

    // Create Eigen matrices to store original and noisy image data
    Eigen::MatrixXd original_matrix(height, width);
    Eigen::MatrixXd noisy_matrix(height, width);

    // Task 2: add a noise in the range [-50;50] to each pixel
    std::random_device dev;
    std::mt19937 rnd{dev()};
    std::uniform_real_distribution<double> dist{-50.0, 50.0};
    double norm = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            const int index = (i * width + j);    // 3 channels (RGB)
            const double original = static_cast<double>(original_image[index]);
            
            // Store original grayscale value
            original_matrix(i, j) = original;
            norm = norm + original*original;
            // Calculate and store noisy value
            double noisy_value = std::clamp(original + dist(rnd), 0.0, 255.0);
            noisy_matrix(i, j) = noisy_value;


        }
    }

    // Reshape matrices to vectors
    Eigen::VectorXd original_vector = Eigen::Map<Eigen::VectorXd>(original_matrix.data(), height * width);
    Eigen::VectorXd noisy_vector = Eigen::Map<Eigen::VectorXd>(noisy_matrix.data(), height * width);

    // Verify that the vectors have the correct length and output norm
    assert(original_vector.size() == height * width);
    assert(noisy_vector.size() == height * width);
    norm=sqrt(norm);
    std::cout << "euclidean norm of original image:" << norm << std::endl;

    //output noisy image
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> noisy_image(height, width);
    noisy_image = noisy_matrix.unaryExpr([](double val) -> unsigned char {
    return static_cast<unsigned char>(val);
     });
    // Save the noisy image
    const std::string noisy_image_path = "noisy.png";
    if (stbi_write_png(noisy_image_path.c_str(), width, height, 1, noisy_image.data(),
                       width) == 0) {
        std::cerr << "Error: Could not save noisy image" << std::endl;
        return 1;
    }
    std::cout << "Noisy image saved to " << noisy_image_path << std::endl;

    stbi_image_free(original_image);

    return 0;
}