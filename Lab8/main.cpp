#include <stdio.h>
#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

const int length = 512;

double** mattodouble(Mat inputImage)
{
    double** ptr = new double* [inputImage.rows];
    for (int i = 0; i < inputImage.rows; ++i)
    {
        ptr[i] = new double [inputImage.cols];
        for (int j = 0; j < inputImage.cols; ++j)
            ptr[i][j] = inputImage.at<uchar>(i, j);
    }
    return ptr;
}

int main(int, char**)
{
    // Mat image = imread("D:/Pictures/Lab8/Wood.0000.jpg", IMREAD_GRAYSCALE);
    Mat image = imread("D:/Pictures/Lab8/Brick.0005.jpg", IMREAD_GRAYSCALE);
    double summ = 0, p_f[256] = {0}, maxB[3] = {0}, index[3] = {0},
    q1 = 0, q2 = 0, prev_q1 = 0,
    u = 0, u1 = 0, u2 = 0,
    sigma = 0, sigma1 = 0, sigma2 = 0, sigmaW = 0, sigmaB = 0;
    double** image_double;

    imshow("Input", image);

    // перевод Mat в double
    image_double = mattodouble(image);

    //частота появления пикселов
    for (int i = 0; i < length; ++i)
        for (int j = 0; j < length; ++j)
            p_f[int(image_double[i][j])] += 1;
    for (int i = 0; i < 256; ++i)
        p_f[i] /= length * length;

    // разделение на 2 кластера
    for (int i = 0; i < 256; i++)
    {
        if (i < 256 / 2)
            q1 += p_f[i];
        else
            q2 += p_f[i];
    }

    // центры кластеров
    for (int i = 0; i < 256; ++i)
    {
        if (i < 256 / 2)
            u1 += (i * p_f[i]) / q1;
        else
            u2 += (i * p_f[i]) / q2;
    }
    u = q1 * u1 + q2 * u2;
    
    // дисперсия распределения
    for (int i = 0; i < 256; ++i)
    {
        if (i < 256 / 2)
            sigma1 += (pow((i - u1), 2) * p_f[i]) / q1;
        else
            sigma2 += (pow((i - u2), 2) * p_f[i]) / q2;
    }

    // суммарная внутрикласстерная дисперсия
    sigmaW = q1 * sigma1 + q2 * sigma2;

    // межкласстерная дисперсия
    sigmaB = q1 * pow((u1 - u), 2) + q2 * pow(u2 - u, 2);

    // дисперсия яркости изображения
    sigma = sigmaW + sigmaB;

    // оптимальное пороговое значение
    q1 = p_f[0];
    q2 = 1 - q1;
    u1 = 1;
    u2 = (u - q1 * u1) / (1 - q1);
    sigmaB = q1 * q2 * pow((u1 - u2), 2);
    maxB[0] = sigmaB;
    for (int i = 1; i < 256; ++i)
    {
        prev_q1 = q1;
        q1 = prev_q1 + p_f[i];
        q2 = 1 - q1;
        u1 = (prev_q1 * u1 + i * p_f[i]) / q1;
        u2 = (u - q1 * u1) / (1 - q1);
        sigmaB = q1 * q2 * pow((u1 - u2), 2);
        if (i < 256 / 3)
        {
            if (sigmaB > maxB[0])
            {
                    maxB[0] = sigmaB;
                    index[0] = i;
            }
        }
        else if (i < (2 * 256) / 3)
        {
            if (sigmaB > maxB[1])
            {
                    maxB[1] = sigmaB;
                    index[1] = i;
            }
        }
        else
        {
            if (sigmaB > maxB[2])
            {
                    maxB[2] = sigmaB;
                    index[2] = i;
            }
        }
    }

    for (int i = 0; i < length; ++i)
        for (int j = 0; j < length; ++j)
        {
            if (image_double[i][j] <= index[0])
                image.at<uchar>(i, j) = 0;
            else if (image_double[i][j] <= index[1])
                image.at<uchar>(i, j) = 100;
            else if (image_double[i][j] <= index[2])
                image.at<uchar>(i, j) = 170;
            else
                image.at<uchar>(i, j) = 255;
        }

    imshow("Output", image);
    // imwrite("D:/Pictures/Lab8/Wood_Result_4_clasters.png", image);
    imwrite("D:/Pictures/Lab8/Brick_Result_4_clasters.png", image);
    waitKey(0);
    destroyAllWindows();
}
