using namespace std;

//nikita.orekhov@gmail.com

#include <iostream>
#include <random>
#include <ctime>
#include <fstream>


long double times1 = 0;//TODELETE

long int const c =  10000; //particles count

//data for changing in calculation process (for one dt)
double x[c];
double y[c];
double z[c];
double Vx[c];
double Vy[c];
double Vz[c];

//constant data, it refresh to first data after calculating (after each dt)
double _x[c];
double _y[c];
double _z[c];
double _Vx[c];
double _Vy[c];
double _Vz[c];

double dt = 0.01; //seconds (the accuracy of integration)
double p = 0.000000001; //p^3 = initial volume for c particles
double T = 293;  //average temperature
double epsilom = 4.2 / (6 * pow(10, 23));
double sigma = 3 * pow(10, -10);
double m = 16 * 1.66 * pow(10, -27); //particle's mass (kg)
double k = 1.38 * pow(10, -23); //Boltzmann constant

//counts for accelerated calculation
double zz = pow(dt, 2) / 2 / m;
double xx = dt / m;

ofstream fout("output.txt");


//initializing all molecules with average speed
//average molecula's kinetic energy mV^2/2 = 3/2kT
void initial() {
	double V = pow((3 * k * T / (m)), 0.5); //TODO
	for (int i = 0; i < c; i++) {
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<> dist0(1,2);
		uniform_int_distribution<> dist1(0,1000);
		double a = dist1(gen);
		uniform_int_distribution<> dist2(0,1000 - a);
		double b = dist2(gen);
		double c = (1000 - a - b);
		a = a / 1000;
		b = b / 1000;
		c = c / 1000;
		Vx[i] = pow(-1,dist0(gen)) * pow(a,0.5) * V;
		Vy[i] = pow(-1,dist0(gen)) * pow(b,0.5) * V;
		Vz[i] = pow(-1,dist0(gen)) * pow(c,0.5) * V;

		_Vx[i] = Vx[i];
		_Vy[i] = Vy[i];
		_Vz[i] = Vz[i];

		uniform_int_distribution<> dist3(0, pow(10,9));
		x[i] = dist3(gen) / pow(10, 9) * p;
		y[i] = dist3(gen) / pow(10, 9) * p;
		z[i] = dist3(gen) / pow(10, 9) * p;

		_x[i] = x[i];
		_y[i] = y[i];
		_z[i] = z[i];
	}
	cout << "initial counting finished" << endl;
}



//particle j acts on particle i
void forceij(int j, int i, double* force) {
	std::clock_t start;
	double duration;
	start = std::clock();

	//rij = ri - rj;
	double rij = pow((pow((_x[i] - _x[j]), 2) + pow((_y[i] - _y[j]), 2) + pow((_z[i] - _z[j]), 2)), 0.5);
	//du/dr
	double dUdr = 4 * epsilom / sigma * (12 * pow(sigma / rij, 13) - 6 * pow(sigma / rij, 7));

	force[0] = (_x[i] - _x[j])*dUdr;
	force[1] = (_y[i] - _y[j])*dUdr;
	force[2] = (_z[i] - _z[j])*dUdr;

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	times1 += duration;
}


void update() {
	cout << "start update" << endl;

	std::clock_t start;
	double duration;
	start = std::clock();

	for (int i = 0; i < c; i++) {
		double f[3];
		x[i] = _x[i] + _Vx[i] * dt;
		y[i] = _y[i] + _Vy[i] * dt;
		z[i] = _z[i] + _Vz[i] * dt;
		for (int j = i + 1; j < c; j++) {
			forceij(j, i, f);

			x[i]  += f[0] * zz;
			y[i]  += f[1] * zz;
			z[i]  += f[2] * zz;
			Vx[i] += f[0] * xx;
			Vy[i] += f[1] * xx;
			Vz[i] += f[2] * xx;

			x[j]  -= f[0] * zz;
			y[j]  -= f[1] * zz;
			z[j]  -= f[2] * zz;
			Vx[j] -= f[0] * xx;
			Vy[j] -= f[1] * xx;
			Vz[j] -= f[2] * xx;
		}
	}

	//refresh constant data
	for (int i = 0; i < c; i++) {
		_x[i] = x[i];
		_y[i] = y[i];
		_z[i] = z[i];
		_Vx[i] = Vx[i];
		_Vy[i] = Vy[i];
		_Vz[i] = Vz[i];
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	cout << "finished update " << duration << endl;
}



int main() {
	cout << "start here" << endl;
	fout << "started" << "\n";

	initial();
	cout << x[0] << endl;
	cout << y[0] << endl;
	cout << z[0] << endl;
	cout << endl;
	cout << Vx[0] << endl;
	cout << Vy[0] << endl;
	cout << Vz[0] << endl;
	cout << endl;

	fout << "particle 0" << endl;
	fout << "x y z" << endl;
	fout << x[0] << endl;
	fout << y[0] << endl;
	fout << z[0] << endl;
	fout << "Vx Vy Vz" << endl;
	fout << Vx[0] << endl;
	fout << Vy[0] << endl;
	fout << Vz[0] << endl;



	for (int k = 0; k < 10; k++) {
		//	update();
		cout << "k = " << k << endl;
		std::clock_t start;
		double duration;
		start = std::clock();

		for (int i = 0; i < c; i++) {
			double f[3];
			x[i] = _x[i] + _Vx[i] * dt;
			y[i] = _y[i] + _Vy[i] * dt;
			z[i] = _z[i] + _Vz[i] * dt;
			for (int j = i + 1; j < c; j++) {

				//start calculate force
				//rij = ri - rj;
				double rij = pow((pow((_x[i] - _x[j]), 2) + pow((_y[i] - _y[j]), 2) + pow((_z[i] - _z[j]), 2)), 0.5);
				//du/dr
				double dUdr = 4 * epsilom / sigma * (12 * pow(sigma / rij, 13) - 6 * pow(sigma / rij, 7));

				f[0] = (_x[i] - _x[j])*dUdr;
				f[1] = (_y[i] - _y[j])*dUdr;
				f[2] = (_z[i] - _z[j])*dUdr;
//				cout << "force[0] " << f[0] << endl;
				//finish calculate force

				x[i] += f[0] * zz;
				y[i] += f[1] * zz;
				z[i] += f[2] * zz;
				Vx[i] += f[0] * xx;
				Vy[i] += f[1] * xx;
				Vz[i] += f[2] * xx;

				x[j] -= f[0] * zz;
				y[j] -= f[1] * zz;
				z[j] -= f[2] * zz;
				Vx[j] -= f[0] * xx;
				Vy[j] -= f[1] * xx;
				Vz[j] -= f[2] * xx;
			}
		}

		//refresh constant data
		for (int i = 0; i < c; i++) {
			_x[i] = x[i];
			_y[i] = y[i];
			_z[i] = z[i];
			_Vx[i] = Vx[i];
			_Vy[i] = Vy[i];
			_Vz[i] = Vz[i];
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	//	cout << "finished update " << duration << endl;
	}






	cout << endl;
	cout << x[0] << endl;
	cout << y[0] << endl;
	cout << z[0] << endl;
	cout << endl;
	cout << Vx[0] << endl;
	cout << Vy[0] << endl;
	cout << Vz[0] << endl;

	fout << "x y z" << endl;
	fout << x[0] << endl;
	fout << y[0] << endl;
	fout << z[0] << endl;
	fout << "Vx Vy Vz" << endl;
	fout << Vx[0] << endl;
	fout << Vy[0] << endl;
	fout << Vz[0] << endl;

	fout.close();

	int m;
	cin >> m;

	return 0;
}