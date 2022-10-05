// RsHomoSepSearch.cpp : This file contains the 'main' function. Program execution begins and ends there.
// by N Mizumoto

// summary
// a separated leader and follower search for each other
// movement from experimental data of termites (Follower/Leader, FF/FM/MM)

// 1 step = 0.2 sec
// CRW (resample every 0.2 sec)

// The simulation program was implemented in Microsoft Visual Studio C++ 2017 (Windows 10).

//// libraries ////
#include <iostream>					// cout; cin
#include <sstream>					// to_string
#include <fstream>					// writing
#include <string>
using namespace std;

///*// opencv (for plot)
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace cv;
//*/

// Random sampling
#include <random>
#include <time.h>
std::uniform_real_distribution<double> distribution(0.0, 1.0);
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
auto rnd = [&] { return distribution(engine); };

#define pi 3.1415926535
#define nipi 3.1415926535*2

//// Paremeters ////
// movement parameters (order is "female after, female before, male after, male before")
// order: FF-F, FM-F, MM-F, FF-L, FM-L, MM-L
double speed[6] = { 9.646121, 10.534230, 12.828930,  2.254564,  3.949095,  4.953810 };
double rho[6] = { 0.7886909, 0.7567764, 0.7500900, 0.6007323, 0.6531174, 0.6274088 };

#define steps 300			// number of time steps (*0.2 = seconds)
#define detection 7			// the distance between two becomes smaller than É”, they encounter (sensitivity analysis: 5.6, 6.3, 7.7, 8.4)
#define d 20.1192//16.766			// female and male are initially separated by the distance d for conditions after separation (sensitivity analysis: 13.4128, 15.0894, 18.4426, 20.1192)
#define study_pairs 100000	// the number of simulated pairs to measure encounter efficiency
//// variables ////
double LocX[2];			// location (x-axis)
double LocY[2];			// location (y-axis)
double Angle[2];		// direction of movement
vector<int> MovementPattern(6);		// put the sets of movement pattern parameters
std::string FileName;	// result name


//// functions ////
// random angle
double dir;
double allocate() {
	dir = rnd() * nipi;
	return(dir);
}

// Wrapped Cauchy distribution (Bartumeus and Levin 2008)
double U, s, Y, Z;
double WrappedCauchy(double rho) {
	U = 1 - rnd();
	Y = 2 * atan(((1 - rho) / (1 + rho)) * tan(pi * (U - 0.5)));
	return(Y);
}

double Rn;
// change position (1:FF-F, 2:FM-F, 3:MM-F, 4:FF-L, 5:FM-L, 6:MM-L)
int change_position(int& j, int pattern) {
	LocX[j] += cos(Angle[j]) * 0.2 * speed[pattern - 1];
	LocY[j] += sin(Angle[j]) * 0.2 * speed[pattern - 1];
	Angle[j] += WrappedCauchy(rho[pattern - 1]);
	return(0);
}


// Reunion search (after separation)
double xdis, ydis, sepdis;
int i, j, simtime, endtime;
double Sep_search(vector<int>& MovementPattern, vector<int>& res, int p) {
	for (i = 0; i < study_pairs; i++) {
		printf("female:%d, male:%d, rep:%d        \r", MovementPattern[0], MovementPattern[1], i + 1);
		// initialization
		LocX[0] = 0;
		LocY[0] = 0;
		LocX[1] = d;
		LocY[1] = 0;

		// determine direction of individuals at the initial condition
		for (j = 0; j < 2; j++) {
			Rn = rnd();
			Angle[j] = allocate();
		}
		endtime = steps;	// record the encountered time for result output
		for (simtime = 0; simtime < steps + 1; simtime++) {
			// change in position
			for (j = 0; j < 2; j++) {
				change_position(j, MovementPattern[j]);
			}
			// encounter determination
			xdis = LocX[1] - LocX[0];
			ydis = LocY[1] - LocY[0];
			sepdis = sqrt(xdis * xdis + ydis * ydis);
			if (sepdis <= detection) {
				endtime = simtime + 1;
				break;
			}
			else if (sepdis > (steps - simtime) / 2) {
				break;
			}
			/*///*
			cv::Mat Img(cv::Size(500, 500), CV_8UC3, cv::Scalar(255, 255, 255));
			cv::namedWindow("search image", cv::WINDOW_AUTOSIZE);
			cv::circle(Img, cv::Point(LocX[0] + (500 / 2), LocY[0] + (500 / 2)), 7, cv::Scalar(0, 0, 200), -1, CV_AA);
			cv::circle(Img, cv::Point(LocX[1] + (500 / 2), LocY[1] + (500 / 2)), 7, cv::Scalar(200, 0, 0), -1, CV_AA);
			cv::imshow("search image", Img);
			cv::waitKey(1);
			//*///*/
		}

		for (j = endtime; j < steps; j++) {
			res[j] ++;
		}
	}
	return(0);
}


//// main simulations ////
int OK, search_mode, emp;
int main() {

	// chech the setting
	cout << "Encounter simulation" << endl;
	cout << "Analysisng steps: 300steps = 60sec" << endl;
	cout << "Detection range:" << detection << " mm" << endl;
	cout << "The separated distance:" << d << "mm" << endl;
	cout << "The encounter probability will be measured from 100,000 pairs" << endl;
	cout << "This 100,000 simulations will be repeated for 10 times" << endl;
	cout << "number indicates : 1 : FF - F, 2 : FM - F, 3 : MM - F, 4 : FF - L, 5 : FM - L, 6 : MM - L" << endl;
	cout << "Are you OK?" << endl;
	cout << "Yes: 1, No: 2" << endl;
	cin >> OK;
	if (OK == 2) { return(2); }

	// output data
	vector<int> res(steps);
	FileName = "SimRes";
	FileName += "_sec";
	FileName += std::to_string(int(steps * 0.2));
	FileName += ".csv";
	std::ofstream ofs(FileName);
	ofs << "Leader,Follower,sim_iter";
	for (i = 0; i < steps; i++) {
		ofs << "," << i + 1;
	}
	ofs << endl;

	// the whole simulation will be repeated for 10 times
	for (int sim_i = 0; sim_i < 10; sim_i++) {
		cout << "simulation iter:" << sim_i << endl;
		// simulations are performed in order of
		// FF-Leader-Follower, FF-Leader-Leader
		// FM-Leader-Follower
		// MM-Leader-Follower, MM-Follower-Follower
		// number indicates: 1:FF-F, 2:FM-F, 3:MM-F, 4:FF-L, 5:FM-L, 6:MM-L
		int Ind1MovementPattern[5] = { 1,4,2,3,3 };	// movement patterns for males
		int Ind2MovementPattern[5] = { 4,4,5,6,3 };	// movement patterns for females
		for (int p = 0; p < 5; p++) {
			// run simulations
			// put movement patterns
			MovementPattern[0] = Ind1MovementPattern[p];
			MovementPattern[1] = Ind2MovementPattern[p];
			fill(res.begin(), res.end(), 0);
			// simulation function for random search (movement patterns, replications, num of steps, container for results, separated distance d)
			Sep_search(MovementPattern, res, p);

			// for result outputs
			cout << endl;
			for (j = 0; j < 2; j++) {
				switch (MovementPattern[j]) {
				case 1:
					ofs << "FF_Follower" << ",";
					break;
				case 2:
					ofs << "FM_Follower" << ",";
					break;
				case 3:
					ofs << "MM_Follower" << ",";
					break;
				case 4:
					ofs << "FF_Leader" << ",";
					break;
				case 5:
					ofs << "FM_Leader" << ",";
					break;
				case 6:
					ofs << "MM_Leader" << ",";
					break;
				}

			}

			ofs << sim_i << ",";
			
			for (i = 0; i < steps - 1; i++) {
				ofs << res[i] << ",";
			}
			ofs << res[steps - 1] << endl;

		}
	}
}