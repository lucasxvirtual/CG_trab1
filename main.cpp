/* 
 * File:   main.cpp
 * Author: jcoelho
 *
 * Created on September 11, 2016, 2:18 PM
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Image.h"
#include "Cluster.h"

using namespace std;



/**
 * TODO: converte uma imagem rgb em uma imagem Lab.
 * @param rgb - imagem rgb de entrada.
 * @return - imagem em Lab.
 */
Image convertImageFromRGB2Lab( const Image& rgb )
{
	int i, j;

	Image Lab = Image(rgb.getW(), rgb.getH());
	for (j = 0; j < rgb.getH(); j++) {

		for (i = 0; i < rgb.getW(); i++) {

			Pixel pixel = rgb.getPixel(i, j);
			pixel = rgb.rgbToXYZ(pixel);
			pixel = rgb.XYZToLab(pixel);
			Lab.setPixel(i, j, pixel);

		}

	}

	return Lab;

}



/**
 * TODO: converte uma imagem Lab em uma imagem em rgb.
 * @param Lab - imagem Lab de entrada.
 * @return - imagem em RGB.
 */
Image convertImageFromLAB2RGB( const Image& Lab )
{

	int i, j;

Image rgb = Image(Lab.getW(), Lab.getH());
for (j = 0; j < Lab.getH(); j++) {

	for (i = 0; i < Lab.getW(); i++) {

		Pixel pixel = Lab.getPixel(i, j);
		pixel = Lab.LabToXYZ(pixel);
		pixel = Lab.XYZTorgb(pixel);
		rgb.setPixel(i, j, pixel);

	}

}

return rgb;

}



/**
 * TODO: inicializa os clusters.
 * @param clusters - clusters que DEVEM ser alocados e inicializados.
 * @param Lab - imagem em Lab.
 * @param k - numero desejado de superpixels.
 * @return - numero de superpixels.
 */
int initializeClusters(Cluster*& clusters, Image& Lab, int k)
{

	int tamSuperPixelx = Lab.getW() / sqrt(k);
	int tamSuperPixely = Lab.getH() / sqrt(k);

	printf("tamanho do super pixel: %d x %d y\ntamanho width: %d\ntamanho height: %d\n", tamSuperPixelx, tamSuperPixely, Lab.getW(), Lab.getH());

	int i, j, l = 0;

	for (i = 0; i < sqrt(k); i++) {

		for (j = 0; j < sqrt(k); j++) {

			l++;

		}

	}

	clusters = new Cluster[l];

	l = 0;

	int a;
	int b;

	for (j = 0; j < sqrt(k); j++) {

		for (i = 0; i < sqrt(k); i++) {
			a = (tamSuperPixelx * i) + (tamSuperPixelx / 2);
			b = (tamSuperPixely * j) + (tamSuperPixely / 2);
			clusters[l] = Cluster(Lab.getPixel((tamSuperPixelx * i) + (tamSuperPixelx / 2), (tamSuperPixely * j) + (tamSuperPixely / 2)), a, b);
			l++;
		}

	}

	return l;

}


double distancia(Pixel& p, Pixel& c, int px, int py, int cx, int cy, double mc, double ms) {
	double dc = sqrt(pow(p[0] - c[0], 2) + pow(p[1] - c[1], 2) + pow(p[2] - c[2], 2));
	double ds = sqrt(pow(px - cx, 2) + pow(py - cy, 2));
	double dt = sqrt(pow(dc / mc, 2) + pow(ds / ms, 2));

	return dt;
}


/**
 * TODO: realiza o algoritmo de superpixels.
 * @param Lab - Imagem em lab.
 * @param clusters - clusters inicializados.
 * @param labels - labels inicializadas.
 * @param k - numero de superpixels.
 * @param M - compacidade.
 */
void performSuperPixelsAlgorithm(Image& Lab, Cluster* clusters, int *labels, int k, double M)
{

	int erro = 6;
	int i = 0;
	int j;
	int tamSuperPixelx = Lab.getW() / sqrt(k);
	int tamSuperPixely = Lab.getH() / sqrt(k);
	int *numPixelSuperPixel = new int[k];
	int *somaPixelx = new int[k];
	int *somaPixely = new int[k];
	int novox, novoy;
	int p;
	for (p = 0; p < k; p++) {
		somaPixelx[p] = 0;
		somaPixely[p] = 0;
		numPixelSuperPixel[p] = 0;
	}

	int janelax1, janelay1, janelax2, janelay2;

	while (i < 7) {

		for (j = 0; j <= k; j++) {

			janelax1 = clusters[j].getX() - tamSuperPixelx;
			if (janelax1 < 0)
				janelax1 = 0;

			janelax2 = clusters[j].getX() + tamSuperPixelx;
			if (janelax2 >= Lab.getW()) {
				janelax2 = Lab.getW() - 1;
			}

			janelay1 = clusters[j].getY() - tamSuperPixely;
			if (janelay1 < 0)
				janelay1 = 0;

			janelay2 = clusters[j].getY() + tamSuperPixely;
			if (janelay2 >= Lab.getH()) {
				janelay2 = Lab.getH() - 1;
			}

			int percorrejanelax;
			int percorrejanelay;

			for (percorrejanelax = janelax1; percorrejanelax <= janelax2; percorrejanelax++) {
				for (percorrejanelay = janelay1; percorrejanelay <= janelay2; percorrejanelay++) {
					Pixel p = Lab.getPixel(percorrejanelax, percorrejanelay);
					if (labels[Lab.computePosition(percorrejanelax, percorrejanelay)] == -1) {
						labels[Lab.computePosition(percorrejanelax, percorrejanelay)] = j;
						somaPixelx[j] += percorrejanelax;
						somaPixely[j] += percorrejanelay;
						numPixelSuperPixel[j]++;
					}
					else {
						double distancia1 = distancia(p, clusters[j].getPixel(), percorrejanelax, percorrejanelay, clusters[j].getX(), clusters[j].getY(), M, tamSuperPixelx*tamSuperPixely);
						double distancia2 = distancia(p, clusters[labels[Lab.computePosition(percorrejanelax, percorrejanelay)]].getPixel(), percorrejanelax, percorrejanelay, clusters[labels[Lab.computePosition(percorrejanelax, percorrejanelay)]].getX(), clusters[labels[Lab.computePosition(percorrejanelax, percorrejanelay)]].getY(), M, tamSuperPixelx*tamSuperPixely);

						if (distancia1 < distancia2) {
							somaPixelx[labels[Lab.computePosition(percorrejanelax, percorrejanelay)]] -= percorrejanelax;
							somaPixely[labels[Lab.computePosition(percorrejanelax, percorrejanelay)]] -= percorrejanelay;
							numPixelSuperPixel[labels[Lab.computePosition(percorrejanelax, percorrejanelay)]]--;
							labels[Lab.computePosition(percorrejanelax, percorrejanelay)] = j;
							somaPixelx[j] += percorrejanelax;
							somaPixely[j] += percorrejanelay;
							numPixelSuperPixel[j]++;
							
						}

					}
				}
			}

		}

		for (p = 0; p < k; p++){
			novox = somaPixelx[p] / numPixelSuperPixel[p];
			novoy = somaPixely[p] / numPixelSuperPixel[p];

			clusters[p].setPixel(Lab.getPixel(novox, novoy));
			clusters[p].setPosition(novox, novoy);
		}

		i++;

	}
	//printf("%d\n", k);
	//for (i = 0; i < Lab.getH()*Lab.getW(); i++) {
		//printf("%d\n", labels[i]);
	//}

}



void drawContoursAroundSegments( Image& rgb, int* labels, Pixel borderColor = Pixel( ) )
{
    int w = rgb.getW( );
    int h = rgb.getH( );

    const int dx8[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };
    const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1, 1 };

    int sz = w * h;
    vector<bool> istaken( sz, false );
    vector<int> contourx( sz );
    vector<int> contoury( sz );
    int mainindex( 0 );
    int cind( 0 );
    for (int j = 0; j < h; j++)
    {
        for (int k = 0; k < w; k++)
        {
            int np( 0 );
            for (int i = 0; i < 8; i++)
            {
                int x = k + dx8[i];
                int y = j + dy8[i];

                if (( x >= 0 && x < w ) && ( y >= 0 && y < h ))
                {
                    int index = y * w + x;

                    if (false == istaken[index])//comment this to obtain internal contours
                    {
                        if (labels[mainindex] != labels[index]) np++;
                    }
                }
            }
            if (np > 1)
            {
                contourx[cind] = k;
                contoury[cind] = j;
                istaken[mainindex] = true;
                //img[mainindex] = color;
                cind++;
            }
            mainindex++;
        }
    }

    int numboundpix = cind; //int(contourx.size());
    for (int j = 0; j < numboundpix; j++)
    {
        for (int n = 0; n < 8; n++)
        {
            int x = contourx[j] + dx8[n];
            int y = contoury[j] + dy8[n];
            if (( x >= 0 && x < w ) && ( y >= 0 && y < h ))
            {
                int ind = rgb.computePosition( x, y );
                if (!istaken[ind])
                {
                    rgb.setPixel( ind, borderColor );
                }
            }
        }
    }
}



void enforceLabelConnectivity( const int* labels, //input labels that need to be corrected to remove stray labels
                               const int width,
                               const int height,
                               int*& nlabels, //new labels
                               int& numlabels, //the number of labels changes in the end if segments are removed
                               const int& K ) //the number of superpixels desired by the user
{
    const int dx4[4] = { -1, 0, 1, 0 };
    const int dy4[4] = { 0, -1, 0, 1 };

    const int sz = width * height;
    const int SUPSZ = sz / K;

    for (int i = 0; i < sz; i++) nlabels[i] = -1;
    int label( 0 );
    int* xvec = new int[sz];
    int* yvec = new int[sz];
    int oindex( 0 );
    int adjlabel( 0 ); //adjacent label
    for (int j = 0; j < height; j++)
    {
        for (int k = 0; k < width; k++)
        {
            if (0 > nlabels[oindex])
            {
                nlabels[oindex] = label;
                //--------------------
                // Start a new segment
                //--------------------
                xvec[0] = k;
                yvec[0] = j;
                //-------------------------------------------------------
                // Quickly find an adjacent label for use later if needed
                //-------------------------------------------------------
                {
                    for (int n = 0; n < 4; n++)
                    {
                        int x = xvec[0] + dx4[n];
                        int y = yvec[0] + dy4[n];
                        if (( x >= 0 && x < width ) && ( y >= 0 && y < height ))
                        {
                            int nindex = y * width + x;
                            if (nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
                        }
                    }
                }

                int count( 1 );
                for (int c = 0; c < count; c++)
                {
                    for (int n = 0; n < 4; n++)
                    {
                        int x = xvec[c] + dx4[n];
                        int y = yvec[c] + dy4[n];

                        if (( x >= 0 && x < width ) && ( y >= 0 && y < height ))
                        {
                            int nindex = y * width + x;

                            if (0 > nlabels[nindex] && labels[oindex] == labels[nindex])
                            {
                                xvec[count] = x;
                                yvec[count] = y;
                                nlabels[nindex] = label;
                                count++;
                            }
                        }

                    }
                }
                //-------------------------------------------------------
                // If segment size is less then a limit, assign an
                // adjacent label found before, and decrement label count.
                //-------------------------------------------------------
                if (count <= SUPSZ >> 2)
                {
                    for (int c = 0; c < count; c++)
                    {
                        int ind = yvec[c] * width + xvec[c];
                        nlabels[ind] = adjlabel;
                    }
                    label--;
                }
                label++;
            }
            oindex++;
        }
    }
    numlabels = label;

    if (xvec) delete [] xvec;
    if (yvec) delete [] yvec;
}



Image SuperPixels( Image& rgb, int k, double M )
{
    //TODO: Converte a imagem para LAb.
	Image Lab = convertImageFromRGB2Lab(rgb);

    //TODO: Calcula o numero de pixels cada superpixel.

    //Todo: Inicializa os os clusters.
    Cluster* clusters;
    
    k = initializeClusters( clusters, Lab, k );


    //TODO: aloca e inicializa labels.
	int* labels = new int[rgb.getH()*rgb.getW()];
	int i = 0;
	
	for (i = 0; i < rgb.getH()*rgb.getW(); i++) {
		labels[i] = -1;
	}

    //TODO: Executa o algoritmo.
	performSuperPixelsAlgorithm(Lab, clusters, labels, k, 1);
    
        //int* nlabels = new int[rgb.getW()*rgb.getH()];
        //enforceLabelConnectivity( labels, rgb.getW(), rgb.getH(), nlabels, k, double(rgb.getW()*rgb.getH()) / double(sqrt(rgb.getW()*rgb.getH()/k) * rgb.getW()*rgb.getH() / k));
        //for (int i = 0; i < rgb.getW()*rgb.getH(); i++)
            //labels[i] = nlabels[i];

    //if (nlabels)
    //    delete [] nlabels;

    //TODO: define as novas cores dos pixels.


    //TODO: Converte a imagem de volta.
	Image rgb2 = convertImageFromLAB2RGB(Lab);

	//printf("%d\n", k);
	//for (i = 0; i < Lab.getH()*Lab.getW(); i++) {
		//printf("%d\n", labels[i]);
	//}

    //Desenha os contornos. Deve passar a imagem em rgb e o vetor de labels.
    drawContoursAroundSegments( rgb2, labels );

	return rgb2;
}



/*
 * 
 */
int main( int argc, char** argv )
{

    Image l;
    if (l.readBMP( "C:\\Users\\Lucas\\Documents\\Visual Studio 2015\\Projects\\CG\\Debug\\AB_ufv_0675.bmp" ))
    {
        printf( "Leitura executada com sucesso\n" );
    }
    
    Image m = SuperPixels( l, 512, 20 );
    
    if (m.writeBMP( "AB_ufv_06752.bmp" ))
    {
        printf( "Escrita executada com sucesso\n" );
    }

    return 0;
}

