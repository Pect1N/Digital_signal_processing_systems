#include <stdio.h>
#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

const int length = 256;
const int N = 20;

double*** Ndcreate(int length)
{
    double*** ptr = new double** [length * 5];
    for (int k = 0; k < length * 5; ++k)
    {
        ptr[k] = new double* [8];
        for (int i = 0; i < 8; ++i)
        {
            ptr[k][i] = new double [8];
            for (int j = 0; j < 8; ++j)
                ptr[k][i][j] = 0;
        }
    }
    return ptr;
}

void quanting (double** image, int length)
{
    for (int i = 0; i < length; ++i)
        for (int j = 0; j < length; ++j)
            for (int k = 0; k < 8; ++k)
                if (image[i][j] < 32 * (k + 1))
                {
                    image[i][j] = k;
                    break;
                }
}

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

void scatteringmerged(double** image, double*** Nd, int k)
{
    for (int f = 1; f < 6; ++f)
        for (int i = 0; i < length * k; ++i)
            for (int j = 0; j < length * k; ++j)
            {
                if (i - f >= 0 && j - f >= 0)
                    ++Nd[0 + 8 * (f - 1)][(int)image[i][j]][(int)image[i - f][j - f]];
                if (i - f >= 0)
                    ++Nd[1 + 8 * (f - 1)][(int)image[i][j]][(int)image[i - f][j]];
                if (i - f >= 0 && j + f < length * k)
                    ++Nd[2 + 8 * (f - 1)][(int)image[i][j]][(int)image[i - f][j + f]];

                if (j - f >= 0)
                    ++Nd[3 + 8 * (f - 1)][(int)image[i][j]][(int)image[i][j - f]];
                if (j + f < length * k)
                    ++Nd[4 + 8 * (f - 1)][(int)image[i][j]][(int)image[i][j + f]];

                if (i + f < length * k && j - f >= 0)
                    ++Nd[5 + 8 * (f - 1)][(int)image[i][j]][(int)image[i + f][j - f]];
                if (i + f < length * k)
                    ++Nd[6 + 8 * (f - 1)][(int)image[i][j]][(int)image[i + f][j]];
                if (i + f < length * k && j + f < length * k)
                    ++Nd[7 + 8 * (f - 1)][(int)image[i][j]][(int)image[i + f][j + f]];
            }
}

void normingmatrixmerged(double*** Nd, double*** p)
{
    double summ = 0;
    for (int f = 1; f < 6; ++f)
        for (int k = 0; k < 8; ++k)
        {
            for (int i = 0; i < 8; ++i)
                for (int j = 0; j < 8; ++j)
                    summ += Nd[k + 8 * (f - 1)][i][j];
            for (int i = 0; i < 8; ++i)
                for (int j = 0; j < 8; ++j)
                    p[k + 8 * (f - 1)][i][j] = Nd[k + 8 * (f - 1)][i][j] / summ;
        }
}

int optimalvectorX2(double*** p)
{
    int optimal = 0;
    double X2 = 0;
    double maxX = 0;
    double pa[40][8] = {0};
    double pb[40][8] = {0};
    for (int k = 0; k < 40; ++k)
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j)
            {
                pa[k][i] += p[k][i][j];
                pb[k][j] += p[k][i][j];
            }
    for (int k = 0; k < 40; ++k)
    {
        X2 = 0;
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j)
                X2 += pow(p[k][i][j], 2) / (pa[k][i] * pb[k][j]);
        X2 -= 1;
        if (X2 > maxX)
        {
            optimal = k;
            maxX = X2;
        }
    }
    return optimal;
}

void scattering(double** image, double*** Nd, int k, int n)
{
    for (int i = 0; i < length * k; ++i)
        for (int j = 0; j < length * k; ++j)
            if (j - 1 >= 0)
                ++Nd[n][(int)image[i][j]][(int)image[i][j - 1]];
}

void normingmatrix(double*** Nd, double*** p)
{
    double summ = 0;
    for (int k = 0; k < 4; ++k)
    {
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j)
                summ += Nd[k][i][j];
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j)
                p[k][i][j] = Nd[k][i][j] / summ;
    }
}

int main(int argc, char** argv )
{
    double** image1;
    double** image2;
    double** image3;
    double** image4;
    double** image5;
    double*** Ndmerged = Ndcreate(8); // 40 (число смещений (8 * 5) 8 (строки) 8 (столбцы)
    double*** pmerged = Ndcreate(8); // 40 (число смещений (8 * 5)) 8 (строки) 8 (столбцы)
    double*** Nd = Ndcreate(1); // 5 (число изображений) 8 (число строк) 8 (число столбцов)
    double*** p = Ndcreate(1); // 5 (число изображений) 8 (строки) 8 (столбцы)
    double mathwaiting[2][5] = {0}; // 2 (a, b) 5 (число изображений)
    double middlequad[2][5] = {0}; // 2 (a, b) 5 (число изображений)
    double texturefeatures[8][5] = {0}; // 8 (число признаков) 5 (число изображений)
    int optimal = 0; // номер оптимального смещения (в нашем случае (3) влево на 1)
    Mat image[5] = {
        imread("D:/Pictures/Lab7/Fabric.0016.png", IMREAD_GRAYSCALE),
        imread("D:/Pictures/Lab7/Fabric.0017.png", IMREAD_GRAYSCALE),
        imread("D:/Pictures/Lab7/Fabric.0018.png", IMREAD_GRAYSCALE),
        imread("D:/Pictures/Lab7/Fabric.0019.png", IMREAD_GRAYSCALE),
        imread("D:/Pictures/Lab7/merged.png", IMREAD_GRAYSCALE)};
    
    imshow("Input", image[4]);

    // перевод изображения в double массив
    image1 = mattodouble(image[0]);
    image2 = mattodouble(image[1]);
    image3 = mattodouble(image[2]);
    image4 = mattodouble(image[3]);
    image5 = mattodouble(image[4]);

    //квантование до 8-ми
    quanting(image1, length * 2);
    quanting(image2, length * 2);
    quanting(image3, length * 2);
    quanting(image4, length * 2);
    quanting(image5, length * 4);

    // определение лучшего смещения
    scatteringmerged(image5, Ndmerged, 4);
    normingmatrixmerged(Ndmerged, pmerged);
    optimal = optimalvectorX2(pmerged);
    Nd[4] = Ndmerged[optimal];

    // формирование матрицы рассеяния (Nd) для текстур
    scattering(image1, Nd, 2, 0);
    scattering(image2, Nd, 2, 1);
    scattering(image3, Nd, 2, 2);
    scattering(image4, Nd, 2, 3);

    // формирование нормированной матрицы рассеяния (p) для текстур
    normingmatrix(Nd, p);
    p[4] = pmerged[optimal];

    // мат.ожидание
    for (int k = 0; k < 5; ++k)
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j)
            {
                mathwaiting[0][k] += i * p[k][i][j];
                mathwaiting[1][k] += j * p[k][i][j];
            }

    // среднее квадратическое
    for (int k = 0; k < 5; ++k)
    {
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j)
            {
                middlequad[0][k] += pow(i, 2) * p[k][i][j] - pow(mathwaiting[0][k], 2);
                middlequad[1][k] += pow(j, 2) * p[k][i][j] - pow(mathwaiting[1][k], 2);
            }
        middlequad[0][k] = sqrt(middlequad[0][k]);
        middlequad[1][k] = sqrt(middlequad[1][k]);
    }

    // признаки текстур
    for (int k = 0; k < 5; ++k)
    {
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j)
            {
                texturefeatures[0][k] += (i - mathwaiting[0][k]) * (j - mathwaiting[1][k]) * p[k][i][j]; // BC
                texturefeatures[1][k] += pow(i - j, 2) * p[k][i][j]; // BI
                texturefeatures[2][k] += abs(i - j) * p[k][i][j]; // BV
                texturefeatures[3][k] += pow(p[k][i][j], 2); // BN
                if (p[k][i][j] > 0)
                    texturefeatures[4][k] += p[k][i][j] * log(p[k][i][j]); // BE
                texturefeatures[5][k] += p[k][i][j] / (1 + pow(i - j, 2)); // BD
                texturefeatures[6][k] += p[k][i][j] / (1 + abs(i - j)); // Bhom
                texturefeatures[7][k] += i * j * p[k][i][j] - mathwaiting[0][k] * mathwaiting[1][k]; // rC
            }
        texturefeatures[4][k] *= -1;
        texturefeatures[7][k] /= middlequad[0][k] * middlequad[1][k];
    }

    // сегментация (работаем с image5)
    for (int i = 0; i < 1024; ++i)
    {
        for (int j = 0; j < 1024; ++j)
        {
            double localfeatures[8] = {0};
            double localNd[8][8] = {0};
            double localp[8][8] = {0};
            double summ = 0;
            double localmwmq[2][2] = {0};
            // локальный Nd
            for (int m = (i - N / 2); m <= i + N / 2; ++m)
                for (int n = (j - N / 2); n <= (j + N / 2); ++n)
                    if (m >= 0 && n >= 0 && m < 1024 && n < 1024)
                        if (n - 1 >= 0)
                        {
                            ++localNd[(int)image5[m][n]][(int)image5[m][n - 1]];
                            ++summ;
                        }
            
            // локальный p
            for (int m = 0; m < 8; ++m)
                for (int n = 0; n < 8; ++n)
                    localp[m][n] = localNd[m][n] / summ;

            // локальное мат.ожидание
            for (int m = 0; m < 8; ++m)
                for (int n = 0; n < 8; ++n)
                {
                    localmwmq[0][0] += m * localp[m][n];
                    localmwmq[0][1] += n * localp[m][n];
                }
            
            // локальное среднее квадратическое
            for (int m = 0; m < 8; ++m)
                for (int n = 0; n < 8; ++n)
                {
                    localmwmq[1][0] += pow(m, 2) * localp[m][n] - pow(localmwmq[0][0], 2);
                    localmwmq[1][1] += pow(n, 2) * localp[m][n] - pow(localmwmq[0][1], 2);
                }
            localmwmq[1][0] = sqrt(localmwmq[1][0]);
            localmwmq[1][1] = sqrt(localmwmq[1][1]);

            // локальные признаки
            for (int m = 0; m < 8; ++m)
                for (int n = 0; n < 8; ++n)
                {
                    localfeatures[0] += (m - localmwmq[0][0]) * (n - localmwmq[0][1]) * localp[m][n];
                    localfeatures[1] += pow(n - m, 2) * localp[m][n];
                    localfeatures[2] += abs(n - m) * localp[m][n];
                    localfeatures[3] += pow(localp[m][n], 2);
                    if (localp[m][n] > 0)
                        localfeatures[4] += localp[m][n] * log(localp[m][n]);
                    localfeatures[5] += localp[m][n] / (1 + pow(m - n, 2));
                    localfeatures[6] += localp[m][n] / (1 + abs(m - n));
                    localfeatures[7] += m * n * localp[m][n] - localp[0][0] * localp[0][1];
                }
            localfeatures[4] *= -1;
            localfeatures[7] /= localp[1][0] * localp[1][1];

            double mineuclid = 1000000000;
            double eucliddist = 0;
            int texturenumber = 0;
            for (int n = 0; n < 8; ++n)
                for (int m = 0; m < 4; ++m)
                {
                    eucliddist = sqrt(pow(texturefeatures[n][m] - localfeatures[n], 2));
                    if (eucliddist < mineuclid)
                    {
                        mineuclid = eucliddist;
                        texturenumber = m;
                    }
                }
            switch (texturenumber)
            {
            case 0:
                image[4].at<uchar>(i, j) = 0;
                break;
            case 1:
                image[4].at<uchar>(i, j) = 100;
                break;
            case 2:
                image[4].at<uchar>(i, j) = 170;
                break;
            case 3:
                image[4].at<uchar>(i, j) = 255;
                break;
            }
        }
    }

    imwrite("D:/Pictures/Lab7/merged_N20.png", image[4]);
    imshow("Result", image[4]);

    waitKey(0);
    destroyAllWindows();
    return 0;
}