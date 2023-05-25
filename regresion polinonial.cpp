#include <iostream>
#include <cmath>

using namespace std;

int tamano;

void eliminacionGauss(double **A, double B[], int n) {
 double inv;
 for (int k = 0; k < n; k++) {
 for (int i = k + 1; i < n; i++) {
 inv = A[i][k] / A[k][k];
 for (int j = k; j < n; j++) {
 A[i][j] = A[i][j] - inv*A[k][j];
 }
 B[i] = B[i] - inv*B[k];
 }
 }
}

void sustitucionAtras(double **A, double B[], int n, double C[]) {
 double suma;
 C[n - 1] = B[n - 1] / A[n - 1][n - 1];
 for (int i = n - 2; i >= 0; i--) {
 suma = 0;
 for (int j = i + 1; j < n; j++) {
 suma = suma + A[i][j] * C[j];
 }
 C[i] = (B[i] - suma) / A[i][i];
 }
}

void regresionPolinomial(double x[], double y[], int n, int m) {
 double sum_x = 0, sum_xy = 0;
 tamano = m + 1;
 double *solucion = new double[tamano];
 double **ecuaciones;

 ecuaciones = new double*[tamano];
 
 for (int i = 0; i < tamano; i++) {
 ecuaciones[i] = new double[tamano];
 }

 for (int i = 0; i < tamano; i++) {
 sum_xy = 0;

 for (int j = 0; j < n; j++)
 sum_xy += pow(x[j], i)*y[j];
 solucion[i] = sum_xy;

 for (int j = 0; j < tamano; j++) {
 sum_x = 0;
 if (i == 0 && j == 0)
 ecuaciones[i][j] = n;
 else {
 for (int h = 0; h < n; h++)
 sum_x += pow(x[h], (j + i));
 ecuaciones[i][j] = sum_x;
 }
 }
 }

 eliminacionGauss(ecuaciones, solucion, tamano);

 double *x_vector = new double[tamano];

 sustitucionAtras(ecuaciones, solucion, tamano, x_vector);

 cout << "La ecuacion es: ";
 for (int i = 0; i < tamano; i++) {
 cout << x_vector[i];
 if (i != 0) {
 cout << "x^" << i;
 }
 if (i != tamano - 1) {
 cout << " + ";
 }
 }
 cout << endl;

 double *e = new double[n];
 for (int i = 0; i < n; i++) {
 double y_calculada = 0;
 for (int j = tamano - 1; j >= 1; j--)
 y_calculada += x_vector[j] * (pow(x[i], j));
 y_calculada += x_vector[0];
 e[i] = pow(y[i] - y_calculada, 2);
 }

 double sum_y = solucion[0];

 double sr = 0;
 double st = 0;
 for (int i = 0; i < n; i++) {
 sr += e[i];
 st += pow(y[i] - (sum_y / n), 2);
 }

 double err = sqrt(sr / (n - tamano));

 double r2 = (st - sr) / st;
 double r = sqrt(r2);
 
 cout << "Error estandar de la estimacion: " << err << endl;
 cout << "Coeficiente de determinacion: " << r2 << endl;
 cout << "Coeficiente de correlacion: " << r << endl;
 cout << endl;

}

int main(){

 double x[] = {-3, -2, -1, 0, 1, 2, 3};
 double y[] = {7.5, 3, 0.5, 1, 3, 6, 14};

 regresionPolinomial(x, y, sizeof(x)/sizeof(x[0]), 4);

 return 0;
   
}