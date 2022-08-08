#include <iostream>
#include <fstream>
#include <complex>
#include <vector>

using namespace std;

double Pi = acos(-1);
double R = 5.0;
double K0 = 1;
//double K = K0 * 1.2;
//double K0 = 1;
double K = K0 * 1.5;


double f(double x) {
    //return 1;
    //return 1 + 31.0 / 3.0 * x * x;
	//return 1+4.0/3.0*x;
    return 1 - x * x / 3.0;
}

double u(double x) {
    return 5 * x * x + 1;
}

double Ker(double x, double s) {
    //return x * s * s - x;
    return x * s + x * x;
}

complex <double> Kernel(double k, double rho_1, double rho_2, double phi_1, double phi_2, bool eq = true) { //если rho_1 = rho_2 то true, иначе false
    complex <double> ed(0, 1.0);
    double l;
    if (eq)
        l = rho_1 * sqrt(2 - 2 * cos(phi_1 - phi_2));
    else
        l = sqrt(pow(rho_1, 2) + pow(rho_2, 2) - 2 * rho_1 * rho_2 * cos(phi_1 - phi_2));
    return exp(ed * k * l) / (4.0 * Pi * l);
}

complex <double> diffKernel(double k, double rho, double phi_1, double phi_2) {
    complex <double> ed(0, 1.0);
    double l = sqrt(2 - 2 * cos(phi_1 - phi_2));
    complex <double> tmp = ed * k * rho * l;
    return exp(tmp) * (tmp - 1.0) / (4.0 * Pi * pow(rho, 2) * l);
}

complex <double> fallWave(double k, double rho, double phi) {
    complex <double> ed(0, 1.0);
    return 5.0*exp(ed * k * rho * cos(phi));
}

complex <double> difffallWave(double k, double rho, double phi) {
    complex <double> ed(0, 1.0);
    complex <double> tmp = ed * k * cos(phi);
    return 5.0 * tmp * exp(tmp * rho);
}

complex <double> timeFunc(double t) {
    complex <double> ed(0, 1.0);
    double omega = K0 * 3 * pow(10, 8);
    return exp(ed * omega * t);
}

void CreateMatrixMemory(int I, int J, complex <double>**& A)
{
    int i1, i2;
    A = new complex <double>* [I];
    for (i1 = 0; i1 < I; i1++) {
        A[i1] = new complex <double>[J];
        for (i2 = 0; i2 < J; i2++) {
            A[i1][i2] = 0.0;
        }
    }
}

//Освобождает память (количество строк)
void DeleteMatrixMemory(int I, complex <double>** A)
{
    int i1;
    for (i1 = 0; i1 < I; i1++) {
        delete A[i1];
    }
    delete[]A;
}

void printMatrix(complex <double>** A, int N, int M) {
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < M; j++) {
            cout << A[i][j] << "  ";
        }
        cout << endl;
    }
}

void Gauss(complex <double>** Matrix, complex <double>* Vec, int Nm) {
    complex <double> ed(0, 1.0);
    complex <double> nul(0.0, 0.0);

    for (int k = 0; k < Nm; k++) {
        if (Matrix[k][k] != ed) {
            complex <double> T = Matrix[k][k];
            for (int j = k; j < Nm; j++) {//нормирование строки
                Matrix[k][j] = Matrix[k][j] / T;
            }
            Vec[k] = Vec[k] / T;
        }
        for (int i = 0; i < Nm; i++) { //проходим по столбцу
            if ((Matrix[i][k] != nul) & (i != k)) {
                complex <double> T = Matrix[i][k];
                Matrix[i][k] = 0;
                for (int j = k + 1; j < Nm; j++) { //проходим по двум строкам и вычитаем их
                    Matrix[i][j] -= Matrix[k][j] * T;
                }
                Vec[i] -= Vec[k] * T;
            }
        }
    }
}

void GenerateVTK(vector<vector<double>> data, string name_file = "out.vtk", double time = -1.0) {

    int points = data.size();
    cout << "Points " << points << endl;

    ofstream file2(name_file);
    file2 << "# vtk DataFile Version 2.0\n" <<
        "Cube example\n" <<
        "ASCII\n" <<
        "DATASET POLYDATA\n" <<
        "POINTS " << points << " float" << endl;

    for (size_t i = 0; i < points; i++)
    {
        file2 << data[i][0] << " " << data[i][1] << " " << data[i][2] << endl;
    }
    int side = sqrt(points);
    int CountOfPolygons = pow(side - 1, 2);
    file2 << "POLYGONS " << CountOfPolygons << " " << CountOfPolygons * 5 << endl;//5 - количество координат на полигон 

    for (size_t i = 0; i < points - side - 1; i++)
    {
        int wP = i;
        if ((wP + 1) % side != 0 || wP == 0)
            file2 << 4 << " " << wP << " " << wP + 1 << " " << wP + side + 1 << " " << wP + side << endl;
    }
    file2 << "POINT_DATA " << points << endl <<
        "SCALARS Magnitude float 1\n" <<
        "LOOKUP_TABLE default" << endl;

    for (size_t i = 0; i < points; i++)
    {
        if (time == -1)
            file2 << data[i][3] << endl;
        else
            file2 << abs(complex<double>(data[i][4], data[i][5]) * timeFunc(time)) << endl;
    }

    file2 <<
        "SCALARS real float 1\n" <<
        "LOOKUP_TABLE default" << endl;

    for (size_t i = 0; i < points; i++)
    {
        if (time == -1)
            file2 << data[i][4] << endl;
        else
            file2 << (complex<double>(data[i][4], data[i][5]) * timeFunc(time)).real() << endl;
    }

    file2 <<
        "SCALARS imag float 1\n" <<
        "LOOKUP_TABLE default" << endl;

    for (size_t i = 0; i < points; i++)
    {
        if (time == -1)
            file2 << data[i][5] << endl;
        else
            file2 << (complex<double>(data[i][4], data[i][5]) * timeFunc(time)).imag() << endl;
    }
    file2.close();
}

complex <double> Integr(double phi_beg, double phi_end, double koll, double K, double rho_1, double rho_2, bool eq = true) {
    int N_int = 3;
    double h_int = (phi_end - phi_beg) / double(N_int);

    complex <double> Sum = 0;
    for (size_t i = 0; i < N_int; i++)
    {
        double s = phi_beg + i * h_int;
        Sum += (Kernel(K, rho_1, rho_2, koll, s, eq) + Kernel(K, rho_1, rho_2, koll, s + h_int, eq)) / 2.0;
    }
    return Sum * h_int;
}

complex <double> diffIntegr(double phi_beg, double phi_end, double koll, double K, double rho) {
    int N_int = 3;
    double h_int = (phi_end - phi_beg) / double(N_int);

    complex <double> Sum = 0;
    for (size_t i = 0; i < N_int; i++)
    {
        double s = phi_beg + i * h_int;
        Sum += (diffKernel(K, rho, koll, s) + diffKernel(K, rho, koll, s + h_int)) / 2.0;
    }
    return Sum * h_int;
}

void ReadFile(const char Patch[], double** Matr, int N, int M) {
    ifstream file(Patch);
    if (!file.is_open()) {
        cout << "Read File " << Patch << " fail\n";
        exit(-3);
    }
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            file >> Matr[i][j];
        }
    }
}

void ReadFile(const char Patch[], double* Vec, int N) {
    ifstream file(Patch);
    if (!file.is_open()) {
        cout << "Read File " << Patch << " fail\n";
        exit(-3);
    }
    for (size_t i = 0; i < N; i++) {
        file >> Vec[i];
    }
}



double xyToPhi(double x, double y) {
    if (x == 0) {
        if (y > 0) return Pi / 2.0;
        else if (y < 0) return 3.0 * Pi / 2.0;
    }
    if (y == 0) {
        if (x > 0) return 0;
        else return Pi;
    }
    if (x > 0 && y > 0)         return atan(y / x);
    else if (x < 0 && y > 0)    return atan(y / x) + Pi;
    else if (x < 0 && y < 0)    return atan(y / x) + Pi;
    else if (x > 0 && y < 0)    return atan(y / x) + 2.0 * Pi;

    printf("X = %f, Y = %f", x, y);
    exit(-2);
}
int main() {
    double A = 0, B = 2 * Pi;
    int n = 200;
    int N = n * 2;
    double h = (B - A) / double(n);
    complex <double>** Am, * Vec = new complex <double>[N], * f_x = new complex <double>[N];

    complex<double>* alpha_vec = new complex <double>[N];
    complex<double>* beta_vec = new complex <double>[N];

    CreateMatrixMemory(N, N, Am);
    for (size_t i = 0; i < n; i++)  //точки коллокации
    {
        double phi_koll = A + i * h + h / 2.0;
        for (size_t j = 0; j < n; j++)  //координаты
        {
            double phi_beg = A + j * h;
            double phi_end = phi_beg + h;

            if (i == j) {
                Am[i + n][j] = -0.5;
                Am[i + n][j + n] = -0.5;
            }
            else {
                Am[i + n][j] = 0.0;
                Am[i + n][j + n] = 0.0;
            }

            Am[i][j] = Integr(phi_beg, phi_end, phi_koll, K0, R, R) * R;       //Якобиан (без него наверное даже лучше?)
            Am[i][j + n] = -Integr(phi_beg, phi_end, phi_koll, K, R, R) * R;
            Am[i + n][j] += diffIntegr(phi_beg, phi_end, phi_koll, K0, R) * R;
            Am[i + n][j + n] -= diffIntegr(phi_beg, phi_end, phi_koll, K, R) * R;
            
            //Am[i][j] -= lambda * Integr(A + j * h, A + (j + 1.0) * h, x_i);
        }
        Vec[i] = fallWave(K0, R, phi_koll);
        Vec[i + n] = difffallWave(K0, R, phi_koll);
    }
    //for (size_t i = 0; i < N; i++)
    //{
    //    for (size_t j = 0; j < N; j++)
    //    {
    //        cout << Am[i][j] << " ";
    //    }
    //    cout << endl;
    //}
    Gauss(Am, Vec, N);
    ofstream alpha("alpha.txt");
    ofstream beta("beta.txt");

    alpha << "X Y Z F real imag" << endl;
    beta << "X Y Z F real imag" << endl;

    for (size_t i = 0; i < n; i++)
    {
        
        double phi = A + i * h + h / 2.0;
        double x = R * cos(phi);
        double y = R * sin(phi);
        alpha_vec[i] = Vec[i];
        beta_vec[i] = Vec[i + n];
        alpha << x << " " << y << " 0 " << abs(Vec[i]) <<" "<< Vec[i].real() << " " << Vec[i].imag() << endl;
        beta << x << " " << y << " 0 " << abs(Vec[i + n]) << " " << Vec[i+n].real() << " " << Vec[i + n].imag() << endl;
        //file << x << " " << Vec[i] << " " << u(x) << endl;
        cout << phi << " " << Vec[i] << endl;
    }

    ofstream map_int("map.txt");
    map_int << "X Y Z abs real imag" << endl;
    //for (double x = -10; x <= 10; x+=0.1)
    //{

    //    for (double y = -10; y <= 10; y += 0.1)
    //    {
    //        double new_rho = sqrt(pow(x, 2) + pow(y, 2));
    //        double new_phi = atan(y / x);
    //        if (x > -0.001 && x<0.001 || y>-0.001 && y < 0.001) {}
    //        else {
    //            
    //            complex<double> Intens(0, 0);

    //            if(new_rho<=R)
    //                for (size_t j = 0; j < n; j++)  //координаты
    //                {
    //                    double phi_beg = A + j * h;
    //                    double phi_end = phi_beg + h;

    //                    Intens += Integr(phi_beg, phi_end, new_phi, K0, R, new_rho, false) * alpha_vec[j] * R + fallWave(K0, new_rho, new_phi);

    //                }
    //            else
    //                for (size_t j = 0; j < n; j++)  //координаты
    //                {
    //                    double phi_beg = A + j * h;
    //                    double phi_end = phi_beg + h;

    //                    Intens += Integr(phi_beg, phi_end, new_phi, K, R, new_rho, false) * beta_vec[j] * R;

    //                }
    //            map_int << x << " " << y << " 0 " << abs(Intens) << " " << Intens.real() << " " << Intens.imag() << endl;
    //        }
    //    }
    //}
    
    ////////////////////////////
    //for (double new_rho = 5.1; new_rho < 20; new_rho += 0.1) {
    //    for (size_t j = 0; j < n; j++)  //координаты
    //    {
    //        double phi_beg = A + j * h;
    //        double phi_end = phi_beg + h;
    //        double phi = A + j * h + h / 2.0;
    //        double x = new_rho * cos(phi);
    //        double y = new_rho * sin(phi);
    //        complex<double> Intens = Integr(phi_beg, phi_end, phi, K0, R, new_rho, false) * alpha_vec[j] * R + fallWave(K0, new_rho, phi);
    //        map_int << x << " " << y << " 0 " << abs(Intens) << " " << Intens.real() << " " << Intens.imag() << endl;
    //    }
    //}
    //////////////////////

    //ofstream file2 = GenerateVTK(n);

    //for (double new_rho = 5.1; new_rho < 10; new_rho += 0.1) {    //самый удачный
    //    for (size_t j = 0; j < n; j++)  //координаты
    //    {
    //        
    //        double phi = A + j * h + h / 2.0;
    //        double x = new_rho * cos(phi);
    //        double y = new_rho * sin(phi);
    //        complex<double> Intens = 0;
    //        for (int ind = 0; ind < n; ind++) {
    //            double phi_beg = A + ind * h;
    //            double phi_end = phi_beg + h;
    //            Intens += Integr(phi_beg, phi_end, phi, K0, R, new_rho, false) * alpha_vec[ind] * R;
    //        }
    //        Intens += fallWave(K0, new_rho, phi);
    //        map_int << x << " " << y << " 0 " << abs(Intens) << " " << Intens.real() << " " << Intens.imag() << endl;
    //    }
    //}

    
    vector <vector<double>> data;
    double delta = 0.1;
    double edge = 20;
    for (double x = -edge; x <= edge; x += delta) {
        for (double y = -edge; y <= edge; y+= delta)  //координаты
        {
            if ((x < -0.01 || x>0.01) && (y < -0.01 || y>0.01)) {
                
                double phi = xyToPhi(x, y);
                double new_rho = sqrt(pow(x, 2) + pow(y, 2));
                if (new_rho > R) {
                    complex<double> Intens = 0;
                    for (int ind = 0; ind < n; ind++) {
                        double phi_beg = A + ind * h;
                        double phi_end = phi_beg + h;
                        Intens += Integr(phi_beg, phi_end, phi, K0, R, new_rho, false) * alpha_vec[ind] * R;
                    }
                    Intens += fallWave(K0, new_rho, phi);
                    map_int << x << " " << y << " 0 " << abs(Intens) << " " << Intens.real() << " " << Intens.imag() << endl;
                    data.push_back({ x,y,0.0, abs(Intens), Intens.real(), Intens.imag() });
                    //file2 << x << " " << y << " 0 " << abs(Intens) << " " << Intens.real() << " " << Intens.imag() << endl;
                }
                else {
                    complex<double> Intens = 0;
                    for (int ind = 0; ind < n; ind++) {
                        double phi_beg = A + ind * h;
                        double phi_end = phi_beg + h;
                        Intens += Integr(phi_beg, phi_end, phi, K, R, new_rho, false) * beta_vec[ind] * R;
                    }
                    //Intens += fallWave(K0, new_rho, phi);
                    map_int << x << " " << y << " 0 " << abs(Intens) << " " << Intens.real() << " " << Intens.imag() << endl;
                    data.push_back({ x,y,0.0, abs(Intens), Intens.real(), Intens.imag() });
                }
            }
        }
    }
    map_int.close();
    
    GenerateVTK(data);

    int count = data.size();
    double time_one_oscl = 2 * Pi / 3 / pow(10, 8);
    int i = 1;
    for (double time = 0; time < 1*time_one_oscl; time+=time_one_oscl/10)
    {
        string name = "time" + to_string(i) + ".vtk";
        GenerateVTK(data, name, time);
        i++;
    }
	return 0;
}
