#include <opencv2/opencv.hpp>
#include <iostream>

using namespace cv;
using namespace std;

Mat image(Mat inputImage)
{
	double T = 0, B = 0, C = 0, D = 0, nx = 0, ny = 0, compress_direction, compress_amount, xn = 0, xdn = 0, yn = 0, ydn = 0, xCenter = 0, yCenter = 0;
	float M = 0, m_numerator = 0, m_denominator = 0;

	//a1-1 a2-1 (3.4)
	for (int x = 0; x < inputImage.cols; ++x)
	{
		for (int y = 0; y < inputImage.rows; ++y)
		{
			xn += x * inputImage.at<uchar>(y, x);
			xdn += inputImage.at<uchar>(y, x);
			yn += y * inputImage.at<uchar>(y, x);
			ydn += inputImage.at<uchar>(y, x);
		}
	}
	xCenter = round(xn / xdn); 
	yCenter = round(yn / ydn);
	Point2f center(xCenter, yCenter);

	for (int x = 0; x < inputImage.cols; ++x)
	{
		for (int y = 0; y < inputImage.rows; ++y)
		{
			T = inputImage.at<uchar>(y, x);
			nx = x - xCenter; // x*
			ny = y - yCenter; // y*
			B += T * (pow(nx, 2) - pow(ny, 2));
			C += T * 2 * nx * ny;
			D += T * (pow(nx, 2) + pow(ny, 2));
		}
	}
	compress_amount = sqrt((D - sqrt(pow(C, 2) + pow(B, 2))) / (D + sqrt(pow(C, 2) + pow(B, 2)))); //mu
	compress_direction = 0.5 * atan(C / B); // tetta

	// 3.16
	for (int x = 0; x < inputImage.cols; ++x)
	{
		for (int y = 0; y < inputImage.rows; ++y)
		{
			m_numerator += inputImage.at<uchar>(y, x) * sqrt(pow(x - xCenter, 2) + pow(y - yCenter, 2));
			m_denominator += inputImage.at<uchar>(y, x);
		}
	}
	m_denominator *= 10;
	M = m_numerator / m_denominator;

	Mat rotate_plus = getRotationMatrix2D(center, compress_direction * (180 / CV_PI), 1);
	Mat rotate_minus = getRotationMatrix2D(center, -compress_direction * (180 / CV_PI), 1);
	Mat compress = getRotationMatrix2D(center, compress_direction, 1 / compress_amount);
	Mat scale = getRotationMatrix2D(center, 0, 1 / M);

	if (atan(C / B) < 1)
	{
		if (B < C)
		{
			compress.at<double>(1, 1) = 1;
			compress.at<double>(1, 2) = 1;
		}
		else
		{
			compress.at<double>(0, 0) = 1;
			compress.at<double>(0, 2) = 1;
		}
	}
	else
	{
		if (B > C)
		{
			compress.at<double>(1, 1) = 1;
			compress.at<double>(1, 2) = 1;
		}
		else
		{
			compress.at<double>(0, 0) = 1;
			compress.at<double>(0, 2) = 1;
		}
	}
	Size size(inputImage.cols, inputImage.rows);

	warpAffine(inputImage, inputImage, rotate_plus, size);
	warpAffine(inputImage, inputImage, compress, size); // a3
	warpAffine(inputImage, inputImage, rotate_minus, size);
	warpAffine(inputImage, inputImage, scale, size); // a5

	return inputImage;
}

Mat Polarcreation(Mat inputImage)
{
	Size size(inputImage.cols, inputImage.rows);
	double xn = 0, xdn = 0, yn = 0, ydn = 0, xCenter = 0, yCenter = 0;
	for (int x = 0; x < inputImage.cols; ++x)
	{
		for (int y = 0; y < inputImage.rows; ++y)
		{
			xn += x * inputImage.at<uchar>(y, x);
			xdn += inputImage.at<uchar>(y, x);
			yn += y * inputImage.at<uchar>(y, x);
			ydn += inputImage.at<uchar>(y, x);
		}
	}
	xCenter = round(xn / xdn); 
	yCenter = round(yn / ydn);
	Point2f center(xCenter, yCenter);
	warpPolar(inputImage, inputImage, size, center, 50, INTER_LINEAR + WARP_POLAR_LINEAR);
	return inputImage;
}

vector<double> qazwsx(Mat inputImage)
{
	Size size(inputImage.cols, inputImage.rows);
	vector<double> priznaki(64);
	vector<double> invariant(64);
	double azimuth = 0, summa = 0;
	for (int x = 1; x < inputImage.cols; ++x)
	{
		for (int y = 1; y < inputImage.rows; ++y)
		{
			summa += inputImage.at<uchar>(y, x);
			azimuth = atan(y/x);
			for (int n = 0; n < 64; ++n)
				if (azimuth >= (2 * CV_PI * n) / 64 && azimuth < (2 * CV_PI * (n + 1)) / 64)
					priznaki[n] += inputImage.at<uchar>(y, x);
		}
	}
	for (int i = 0; i < 64; ++i)
		priznaki[i] = priznaki[i] / summa;
	for (int i = 0; i < 64; ++i)
		for (int j = 0; j < 64; ++j)
			invariant[i] += abs(priznaki[j] * exp((-2 * CV_PI * j * i) / 64));
	return invariant;
}

int main()
{
	Mat inputImage1 = imread("D:/Pictures/Lab5/image_1_14.jpg", IMREAD_UNCHANGED);
	Mat inputImage2 = imread("D:/Pictures/Lab5/image_2_14.jpg", IMREAD_UNCHANGED);
	Mat inputImage3 = imread("D:/Pictures/Lab5/image_3_14.jpg", IMREAD_UNCHANGED);
	Mat inputImage4 = imread("D:/Pictures/Lab5/image_4_14.jpg", IMREAD_UNCHANGED);
	imshow("Input image1", inputImage1);
	imshow("Input image2", inputImage2);
	imshow("Input image3", inputImage3);
	imshow("Input image4", inputImage4);
	if (inputImage1.empty() || inputImage2.empty() || inputImage3.empty() || inputImage4.empty())
		return -1;

	imshow("Default image1", image(inputImage1));
	imshow("Default image2", image(inputImage2));
	imshow("Default image3", image(inputImage3));
	imshow("Default image4", image(inputImage4));

	inputImage1 = imread("D:/Pictures/Lab5/image_1_14ns.jpg", IMREAD_UNCHANGED);
	inputImage2 = imread("D:/Pictures/Lab5/image_2_14ns.jpg", IMREAD_UNCHANGED);
	inputImage3 = imread("D:/Pictures/Lab5/image_3_14ns.jpg", IMREAD_UNCHANGED);
	inputImage4 = imread("D:/Pictures/Lab5/image_4_14ns.jpg", IMREAD_UNCHANGED);
	if (inputImage1.empty() || inputImage2.empty() || inputImage3.empty() || inputImage4.empty())
		return -1;

	imshow("Default image1ns", Polarcreation(inputImage1));
	imshow("Default image2ns", Polarcreation(inputImage2));
	imshow("Default image3ns", Polarcreation(inputImage3));
	imshow("Default image4ns", Polarcreation(inputImage4));

	vector<double> invariant1(64);
	vector<double> invariant2(64);
	vector<double> invariant3(64);
	vector<double> invariant4(64);
	invariant1 = qazwsx(Polarcreation(inputImage1));
	invariant2 = qazwsx(Polarcreation(inputImage2));
	invariant3 = qazwsx(Polarcreation(inputImage3));
	invariant4 = qazwsx(Polarcreation(inputImage4));

	double summ1 = 0, summ2 = 0, summ3 = 0;
	for (int i = 0; i < 64; ++i)
	{
		summ1 += pow((invariant4[i] - invariant1[i]), 2);
		summ2 += pow((invariant4[i] - invariant2[i]), 2);
		summ3 += pow((invariant4[i] - invariant3[i]), 2);
	}
	summ1 = sqrt(summ1);
	summ2 = sqrt(summ2);
	summ3 = sqrt(summ3);

	if (summ1 < summ2)
	{
		if (summ1 < summ3)
			std::cout << "Included in 1 class" << std::endl;
		else
			std::cout << "Included in 3 class" << std::endl;
	}
	else
	{
		if (summ2 < summ3)
			std::cout << "Included in 2 class" << std::endl;
		else
			std::cout << "Included in 3 class" << std::endl;
	}

	waitKey(0);
    destroyAllWindows();
	return 0;
}