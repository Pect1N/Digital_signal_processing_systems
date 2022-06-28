#include <stdio.h>
#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

const int N = 15;
double save = 0;

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

double** localbright(double** inputImage, int a)
{
    double** ptr = new double* [a];
    for (int i = 0; i < a; ++i)
    {
        ptr[i] = new double [a];
        for (int j = 0; j < a; ++j)
            ptr[i][j] = inputImage[i][j];
    }
    int flag = 0;

    for (int i = 0; i < a; ++i)
        for (int j = 0; j < a; ++j)
        {
            save = 0;
            flag = 0;
            for (int m = (i - N / 2); m <= (i + N / 2); ++m)
                for (int n = (j - N / 2); n <= (j + N / 2); ++n)
                    if (m >= 0 && n >= 0 && m < a && n < a)
                    {
                        save = save + inputImage[m][n];
                        flag++;
                    }
            ptr[i][j] = inputImage[i][j] - save / flag;
        }
    return ptr;
}

double** filtration(double** inputImage, int a, int* L5, int* S5)
{
    int M[5][5];
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            M[i][j] = L5[i] * S5[j];

    double** ptr = new double* [a];
    for (int i = 0; i < a; ++i)
    {
        ptr[i] = new double [a];
        for (int j = 0; j < a; ++j)
            ptr[i][j] = inputImage[i][j];
    }

    for (int i = 0; i < a; ++i)
        for (int j = 0; j < a; ++j)
        {
            save = 0;
            for (int u = -2; u <= 2; ++u)
                for (int v = -2; v <= 2; ++v)
                    if ((i + u) >= 0 && (j + v) >= 0 && (i + u) < a && (j + v) < a)
                        save += inputImage[i + u][j + v] * M[u + 2][v + 2];
            ptr[i][j] = save;
        }
    return ptr;
}

double** energytexture(double** inputImage, int a)
{
    double** ptr = new double* [a];
    for (int i = 0; i < a; ++i)
    {
        ptr[i] = new double [a];
        for (int j = 0; j < a; ++j)
            ptr[i][j] = inputImage[i][j];
    }

    for (int i = 0; i < a; ++i)
        for (int j = 0; j < a; ++j)
        {
            save = 0;
            for (int r = -7; r <= 7; ++r)
                for (int c = -7; c <= 7; ++c)
                    if ((i + r) >= 0 && (j + c) >= 0 && (i + r) < a && (j + c) < a)
                        save += abs(inputImage[i + r][j + c]);
            ptr[i][j] = save;
        }
    return ptr;
}

double** tfeature(double** Ekl, double** Elk, int a)
{
    double** ptr = new double* [a];
    for (int i = 0; i < a; ++i)
    {
        ptr[i] = new double [a];
        for (int j = 0; j < a; ++j)
            ptr[i][j] = Ekl[i][j];
    }

    for (int i = 0; i < a; ++i)
        for (int j = 0; j < a; ++j)
            ptr[i][j] = (Ekl[i][j] + Elk[i][j]) / 2;
    return ptr;
}

double** mediate(double** inputImage)
{
    int a = 1024;
    double** ptr = new double* [a];
    for (int i = 0; i < a; ++i)
    {
        ptr[i] = new double [a];
        for (int j = 0; j < a; ++j)
            ptr[i][j] = inputImage[i][j];
    }
    int flag = 0;

    for (int i = 0; i < 1024; ++i)
        for (int j = 0; j < 1024; ++j)
        {
            save = 0;
            flag = 0;
            for (int m = (i - N / 2); m <= (i + N / 2); ++m)
                for (int n = (j - N / 2); n <= (j + N / 2); ++n)
                    if (m >= 0 && n >= 0 && m < 1024 && n < 1024)
                    {
                        save = save + inputImage[m][n];
                        flag++;
                    }
            ptr[i][j] = save / flag;
        }
    return ptr;
}

int main(int argc, char** argv)
{
    double** image1;
    double** image2;
    double** image3;
    double** image4;
    double** image5;
    double** Ffilter[5][16];
    double** Emap[5][16];
    double** E_texturefeature[5][9];
    double save = 0, mathwaiting[5][9], middlequad[5][9], distance[12][9], effectivefeatures[2][9];
    double a = 0, b = 512*512, max = 0;
    a = 1 / (b * (b - 1));
    int vectors[4][5] = {{1, 4, 6, 4, 1}, {-1, 0, 2, 0, -1}, {-1, -2, 0, 2, 1}, {1, -4, 6, -4, 1}}, counter = 0, 
    L5[5] = {1, 4, 6, 4, 1}, E5[5] = {-1, 0, 2, 0, -1}, S5[5] = {-1, -2, 0, 2, 1}, R5[5] = {1, -4, 6, -4, 1};
    int length = 512;
    Mat images[5];
    string paths[5] = 
        {
            "D:/Pictures/Lab6/Bark.0000.png",
            "D:/Pictures/Lab6/Brick.0004.png",
            "D:/Pictures/Lab6/Fabric.0014.png",
            "D:/Pictures/Lab6/Grass.0001.png",
            "D:/Pictures/Lab6/nmerged.png"
        };
    images[0] = imread(paths[0], IMREAD_GRAYSCALE);
	images[1] = imread(paths[1], IMREAD_GRAYSCALE);
	images[2] = imread(paths[2], IMREAD_GRAYSCALE);
	images[3] = imread(paths[3], IMREAD_GRAYSCALE);
    images[4] = imread(paths[4], IMREAD_GRAYSCALE);
    if (images[0].empty() || images[1].empty() || images[2].empty() || images[3].empty()|| images[4].empty())
		return -1;
    // imshow("Input image", images[0]);
	// imshow("Input image2", images[1]);
	// imshow("Input image3", images[2]);
	// imshow("Input image4", images[3]);
    imshow("Image", images[4]);

    image1 = mattodouble(images[0]);
    image2 = mattodouble(images[1]);
    image3 = mattodouble(images[2]);
    image4 = mattodouble(images[3]);
    image5 = mattodouble(images[4]);

    image1 = localbright(image1, 512);
    image2 = localbright(image2, 512);
    image3 = localbright(image3, 512);
    image4 = localbright(image4, 512);
    image5 = localbright(image5, 1024);
    // imshow("Localbright (1.1)", images[0]);
    // imshow("Localbright1 (1.1)", images[1]);
    // imshow("Localbright2 (1.1)", images[2]);
    // imshow("Localbright3 (1.1)", images[3]);
    // imshow("Localbright4 (1.1)", images[4]);

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        {
            Ffilter[0][i * 4 + j] = filtration(image1, 512, vectors[i], vectors[j]);
            Ffilter[1][i * 4 + j] = filtration(image2, 512, vectors[i], vectors[j]);
            Ffilter[2][i * 4 + j] = filtration(image3, 512, vectors[i], vectors[j]);
            Ffilter[3][i * 4 + j] = filtration(image4, 512, vectors[i], vectors[j]);
            Ffilter[4][i * 4 + j] = filtration(image5, 1024, vectors[i], vectors[j]);
        }
    // imshow("Lavs filtration (1.2)", Ffilter[4][0]);

    for (int i = 0; i < 5; ++i)
    {
        if (i == 4)
            length = 1024;
        else
            length = 512;
        for (int j = 0; j < 16; ++j)
            Emap[i][j] = energytexture(Ffilter[i][j], length);
    }
    // imshow("Energy texture (1.3)", Emap[4][1]);
    // imshow("Energy texture1 (1.3)", Emap[4][4]);

    for (int i = 0; i < 5; ++i)
    {
        if (i == 4)
            length = 1024;
        else
            length = 512;
        E_texturefeature[i][0] = tfeature(Emap[i][1], Emap[i][4], length);
        E_texturefeature[i][1] = tfeature(Emap[i][3], Emap[i][12], length);
        E_texturefeature[i][2] = tfeature(Emap[i][1], Emap[i][4], length);
        E_texturefeature[i][3] = tfeature(Emap[i][5], Emap[i][5], length);
        E_texturefeature[i][4] = tfeature(Emap[i][2], Emap[i][8], length);
        E_texturefeature[i][5] = tfeature(Emap[i][6], Emap[i][9], length);
        E_texturefeature[i][6] = tfeature(Emap[i][10], Emap[i][10], length);
        E_texturefeature[i][7] = tfeature(Emap[i][14], Emap[i][11], length);
        E_texturefeature[i][8] = tfeature(Emap[i][15], Emap[i][15], length);
    }
    // imshow("Texture features (1.4)", E_texturefeature[4][0]);

    // мат.ожидание и среднеквадратическое
    for (int i = 0; i < 5; ++i)
    {
        if (i == 4)
            length = 1024;
        else
            length = 512;
        for (int j = 0; j < 9; ++j)
        {
            for (int y = 0; y < length; ++y)
                for (int x = 0; x < length; ++x)
                    save += E_texturefeature[i][j][y][x];
            mathwaiting[i][j] = save / (length * length);
            save = 0;
            for (int y = 0; y < length; ++y)
                for (int x = 0; x < length; ++x)
                    save += pow(E_texturefeature[i][j][y][x] - mathwaiting[i][j], 2);
            middlequad[i][j] = sqrt(save / a);
            save = 0;
        }
    }

    // расстояние
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (i != j)
                for (int p = 0; p < 9; ++p)
                    distance[i * 4 + j - counter][p] = (mathwaiting[i][p] - mathwaiting[j][p]) / (middlequad[i][p] + middlequad[j][p]);
            else
                counter++;

    // определение лучших
    for (int i = 0; i < 9; ++i)
    {
        max = distance[0][i];
        for (int j = 0; j < 12; ++j)
            if (distance[j][i] > max)
                max = distance[j][i];
        effectivefeatures[0][i] = max;
        effectivefeatures[1][i] = i;
    }

    // сортировка по убыванию
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8 - i; ++j)
            if (effectivefeatures[0][j] < effectivefeatures[0][j + 1])
            {
                double tmp = 0, tmp1 = 0;
                tmp = effectivefeatures[0][j];
                effectivefeatures[0][j] = effectivefeatures[0][j + 1];
                effectivefeatures[0][j + 1] = tmp;
                tmp1 = effectivefeatures[1][j];
                effectivefeatures[1][j] = effectivefeatures[1][j + 1];
                effectivefeatures[1][j + 1] = tmp1;
            }

    // среднее значение
    for (int i = 0; i < 9; ++i)
        E_texturefeature[4][i] = mediate(E_texturefeature[4][i]);

    const int K = 9;
    for (int i = 0; i < 1024; i++)
        for (int j = 0; j < 1024; j++)
        {
            auto tmp = 0.l;
            for (int l = 0; l < K; l++)
            {
                double tmp2;
                int f = effectivefeatures[1][l];
                tmp2 = E_texturefeature[4][f][i][j] - mathwaiting[0][f];
                tmp += tmp2 * tmp2;
            }
            tmp = sqrt(tmp);
            auto k_min = 0;
            auto min_dist = tmp;
            for (int k = 1; k < 4; k++)
            {
                auto tmp = 0.l;
                for (int l = 0; l < K; l++)
                {
                    double tmp2;
                    int f = effectivefeatures[1][l];
                    tmp2 = E_texturefeature[4][f][i][j] - mathwaiting[k][f];
                    tmp += tmp2 * tmp2;
                }
                tmp = sqrt(tmp);
                if (tmp < min_dist)
                {
                    min_dist = tmp;
                    k_min = k;
                }
            }
            switch (k_min)
            {
            case 0:
                images[4].at<uchar>(i, j) = 0;
                break;
            case 1:
                images[4].at<uchar>(i, j) = 100;
                break;
            case 2:
                images[4].at<uchar>(i, j) = 170;
                break;
            case 3:
                images[4].at<uchar>(i, j) = 255;
                break;
            }
        }
    imshow("Result", images[4]);
    //imwrite("D:/Pictures/Lab6/result_N15.png", images[4]);

    waitKey(0);
    destroyAllWindows();
    return 0;
}