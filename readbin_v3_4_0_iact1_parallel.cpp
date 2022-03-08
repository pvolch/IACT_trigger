//от 3.4 отличается отсутствием приведения спектра
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <utility>
#include <chrono>
#include <random>
#include <vector>
#include <algorithm>
#include <map>
#include <omp.h>
#include "Math/Interpolator.h"
using namespace std;
#define Pi 3.14159265359
FILE *f5;
double inf = numeric_limits<double>::infinity();
int32_t N_run, i_scat, i_tel, N_b_f, N_shower, N_phe, read32[4], history, roco[2];
double buf_header_angles[5], phe[2], tim_min, tim_max, amp;
double print_buf_head[15], readdubble[20], rowcol[560][5], integral35, neigh_clean[560][13], hold, trig_time, amp1, amp2;
string line;
int neigh[560][10], trigger, trig_clust, trig_hold, grid_max, cluster, real_numb[560][2];
double t_grid = 1; //время между записями во временной сетке
double t_sign = 20; //время в течении которого считаем что сигнал от фотоэлектрона значим
double trig_amp = 8.2; //превышение данной амплитуды дает холд в подтягиваемом кластере и окно 80 нс и 15 нс для холдов и второго пикселя в этом кластере соответственно
double t_hold = 160; //длинна сигнала HOLD
double integrate_window = 48; //время, в течении которого суммируем фотоэлектроны
double trig_window = 15; //окно триггера
double t_after_before = 100; //отступ от первого черенковского фэ и после последнего черенковского фэ
int number_of_pixels = 560, ev = 91, pix1, pix2;
int numb_of_clusters = 22; //потому что массив от нуля а у нас в файле запись от 1
double delta_e_range = 10;
int list_numb = 0;
int list_scat = 0;
int hidmir[5][2] = {{3, -1},{3, 0},{3, 1},{3, 2},{2, 0}};
typedef pair <int, int> type_key;
map <type_key, int> Mapa;
type_key key[560];
double tet_center, fi_center, tet_source, fi_source;
ROOT::Math::Interpolator inter (73, ROOT::Math::Interpolation::kAKIMA_PERIODIC);

//double fy = 0, fb = 1, fa = 6., fx = 0, fc = 2.03, fd = 2.5;
//double fymx = 0;
//map <pair, int> pmt_coord;
//vector<Photoelectron> pe_focal;
//map<int32_t, int32_t> pixN;
//vector<int32_t> a_hs;
//vector<double> t_hs;
//vector<double> dt_hs;
int32_t a_final;
int mirrow, mircol, bad_mirror, miss_phe, evch = 0;

void readSlowPulse(){
	ifstream file0("slow_pulse.txt");
	if (!file0.is_open()) {
		cout << "Файл slow_pulse.txt не найден" << endl;
	}
	else {
		double pulse[73][3];
		vector<double> vector_x;
		vector<double> vector_y;
		for(int ii = 0; ii < 73; ii++) {
			getline(file0, line);
			stringstream ist(line);
			for(int i = 0; i < 3; i++) {
				ist >> pulse[ii][i];
				cout << pulse[ii][i] << "\t";
			}
			vector_x.push_back(pulse[ii][0] + 10);
			vector_y.push_back(pulse[ii][2]/868);
			cout << endl;
		}
		inter.SetData(vector_x, vector_y);
	}
	file0.close();
}

int mod(int n, int d)
{
	int result = n % d;
	if ((result * d) < 0)
		result += d;
	return result;
}
void readRowColClust(){
	ifstream file0("RowColClust560_fix_375.txt");
	if (!file0.is_open()) {
		cout << "Файл не найден" << endl;
	}
	else {
		for(int ii = 0; ii < number_of_pixels; ii++) {
			getline(file0, line);
			stringstream ist(line);
			for(int i = 0; i < 5; i++) {
				ist >> rowcol[ii][i];
				//cout << rowcol[ii][i] << "\t";
			}
			rowcol[ii][2] = rowcol[ii][2] - 1;
			cout << ii << "\t" << rowcol[ii][2] << endl;
			key[ii] = make_pair(rowcol[ii][0], rowcol[ii][1]);
			Mapa[key[ii]] = ii;
		}
	}
	file0.close();
}

void readNeighbours(){
	ifstream file0("neighbour_cluser_pairs_iact1");
	if (!file0.is_open()) {
		cout << "Файл не найден" << endl;
	}
	else {
		for(int ii = 0; ii < number_of_pixels; ii++) {
			getline(file0, line);
			stringstream ist(line);
			for(int i = 0; i < 9; i++) {
				ist >> neigh[ii][i];
				/*if(neigh[ii][i] == -1){
				   neigh[ii][i] = ii;
				   }*/
			}
		}
	}
	file0.close();
}

void readNeighbours_clean(){
	ifstream file1("neighbour_clean_corsica_iact1_375.txt");
	if (!file1.is_open()) {
		cout << "Файл не найден" << endl;
	}
	else {
		for(int ii = 0; ii < number_of_pixels; ii++) {
			getline(file1, line);
			stringstream ist(line);
			for(int i = 0; i < 13; i++) {
				ist >> neigh_clean[ii][i];
			}
		}
	}
	file1.close();
}

void readChNumbs_clean(){
	ifstream file4("neighbour_clean_experiment_iact1_375.txt");
	if (!file4.is_open()) {
		cout << "neighbour_clean Файл не найден" << endl;
	}
	else {
		for(int ii = 0; ii < number_of_pixels; ii++) {
			getline(file4, line);
			stringstream ist(line);
			for(int i = 0; i < 2; i++) {
				ist >> real_numb[ii][i];
			}
		}
	}
	file4.close();
}

vector<double> random_poisson(double time_begin, double time_end)
{
	double mean_ph = 3; //mean number of NSB photons
	double mean_ph_time = 35; //period with 2 NSB photons
	vector<double> result;
	// construct a trivial random generator engine from a time-based seed:
	long int seed = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator (seed);
	poisson_distribution<int> distribution (mean_ph*(time_end-(time_begin-500))/mean_ph_time);
	double numb_NSB_phe = distribution(generator);
	random_device myRandomDevice;
	unsigned seed2 = myRandomDevice();
	default_random_engine generator2 (seed2);
	uniform_real_distribution<double> distribution2(time_begin-500, time_end);
	for(int t = 0; t <= round(numb_NSB_phe); t++) {
		result.push_back(distribution2(generator2));
	}
	return result;
}

double amp_t(double t, double ampmx){
	double y = 0, b = 1, a = 6., c = 2.03, d = 2.5;
	double x = t - a;
	double ymx = b * exp(-pow((c/(2*d)),2))*exp(pow(c,2)/(2*pow(d,2)));
	if(t < a) {
		y=b*exp(-pow((x/c),2))*exp(-x/d)*ampmx/ymx;
	}
	else{
		y = b*exp(-x/d)*ampmx/ymx;
	}
	return y;
}

vector<float> amp_distr(string path)
{
	vector<float> result;
	ifstream file(path);
	if (!file.is_open()) {
		cout << "Файл распределения амлитуд не найден" << endl;
		exit(0);
	}
	else {
		getline(file, line);
		stringstream ist(line);
		while(ist) {
			float amp = 0;
			ist >> amp;
			if(amp > 0.6) {
				result.push_back(amp);
			}
		}
		file.close();
	}
	return result;
}

int main()
{
	readSlowPulse();
	readRowColClust();
	readNeighbours_clean();
	readNeighbours();
	readChNumbs_clean();
	//read32: N_shower, i_scat, i_tel, N_phe
	//readdubble: energy, zen, az, x_core, y_core, z_core, height_1_interaction, particle_type,
	//xmax, hmax, xtel, ytel, ztel, Xoffsetmm, Yoffsetmm, Thetatelrad, Phitelrad, alfarad, alfamaxrad,
	//T
	///k1/taiga_pool/CORSIKA_TAIGA/output/bpe551_30_da0.1_md5/taiga551_feb
	vector <float> vector_amp_phe;
	cout << "getting vector of phe amps" << endl;
	vector_amp_phe = amp_distr("./probablies8.txt");
	cout << "vector complitely obtaned" << endl;
	ofstream fout("/home/volchugov/corsika/data/zen_0-30/taiga631_iact1.txt");
	f5 = fopen("/k1/taiga_pool/CORSIKA_TAIGA/mc/bpe631_30_da5.0_md5/taiga631_feb", "rf");
	int event = 0;
	while(1) {
		//cout << event << endl;
		miss_phe = 0;
		tim_min = inf;
		tim_max = -inf;
		fread(&read32, sizeof(int32_t), 4, f5);         //read32[0] number shower, read32[1] - shower scattering ID

		//read32[2] - tel ID
		if(feof(f5) != 0) {
			cout << "Файл закончился или пустой" << endl;
			break;
		}
		if(read32[3] < 0) {
			cout << "Wrong number photoelectrons: " << read32[3] << endl;
			break;
		}
		//cout << read32[0] << "\t" << read32[1] << "\t" << read32[2] << "\t" << list_numb << "\t" << list_scat << endl;
		if(read32[2] == 0 && read32[3] >= 0) {
			vector<vector<double> > vector_phe_time(number_of_pixels);
			vector<vector<double> > vector_amp(number_of_pixels);
			vector<vector<double> > vector_amp_MC(number_of_pixels);
			vector<vector<double> > vector_integral35(number_of_pixels);
			vector<vector<double> > vector_integral_hold(numb_of_clusters);
			vector<vector<double> > vector_time10(number_of_pixels);
			vector<vector<double> > vector_prelim_hold_time(numb_of_clusters);
			vector<vector<int> > vector_prelim_hold_pix(numb_of_clusters);
			fread(&readdubble, sizeof(double), 20, f5);         // readdubble[0]/1E12 - energy, readdubble[1]*180/Pi - Thetarad
			//readdubble[2] - Phirad, readdubble[3], readdubble[4], readdubble[5] - XYZ, readdubble[6]/1E3 - higth(m),
			//(readdubble[12] - Ztel, readdubble[19] - time(sec)
			//double R = pow(pow((readdubble[10] + readdubble[13] - readdubble[3])/1E3, 2) + pow((readdubble[11] + readdubble[14] - readdubble[4])/1E3, 2), 0.5);
			//cout << event << "\t" << readdubble[10] << "\t" << readdubble[11] << "\t" << readdubble[12] << "\t" << readdubble[13] << "\t" << readdubble[14] << "\t" << readdubble[3] << "\t" << readdubble[4] << "\t" << readdubble[5] << endl;
			//cout << event << "\t" << read32[0] << "\t" << read32[1] << "\t" << read32[2] << "\t" << readdubble[0]/1E12  << endl;
			//cout << read32[0] << "\t" << read32[1] << "\t" << read32[2] << "\t" << read32[3] << "\t" << readdubble[7] << endl;
			tet_center = readdubble[15];
			fi_center =  readdubble[16];
			tet_source = readdubble[1];
			fi_source =  readdubble[2];
			/*double xToSource=1*sin(tet_source)*cos(fi_source);
			   double yToSource=1*sin(tet_source)*sin(fi_source);
			   double zToSource=1*cos(tet_source);

			   double xToCenter=1*sin(tet_center)*cos(-fi_center);
			   double yToCenter=1*sin(tet_center)*sin(-fi_center);
			   double zToCenter=1*cos(tet_center);

			   double thetamir_c= (180./Pi)*acos((xToSource*xToCenter+yToSource*yToCenter+zToSource*zToCenter));
			   double axis[3] = {-sin(fi_center), cos(fi_center), 0};

			   double M[3][3]={{(1-cos(tet_center))*axis[0]*axis[0]+cos(tet_center),    (1-cos(tet_center))*axis[1]*axis[0]+sin(tet_center)*axis[2],  (1-cos(tet_center))*axis[2]*axis[0]-sin(tet_center)*axis[1]},
			        {(1-cos(tet_center))*axis[0]*axis[1]-sin(tet_center)*axis[2], (1-cos(tet_center))*axis[1]*axis[1]+cos(tet_center),     (1-cos(tet_center))*axis[2]*axis[1]+sin(tet_center)*axis[0]},
			        {(1-cos(tet_center))*axis[0]*axis[2]+sin(tet_center)*axis[1],  (1-cos(tet_center))*axis[1]*axis[2]- sin(tet_center)*axis[0], (1-cos(tet_center))*axis[2]*axis[2]+cos(tet_center)}};

			   double Ptheta[3]={axis[1]*sin(tet_center), -axis[0]*sin(tet_center), cos(tet_center)};
			   double XYZmir_c[3];
			   for(int i=0; i<3; i++) {
			        xToSource=1*sin(tet_source)*cos(-fi_source);
			        yToSource=1*sin(tet_source)*sin(-fi_source);
			        zToSource=1*cos(tet_source);
			        XYZmir_c[i] = M[i][0]*xToSource + M[i][1]*yToSource + M[i][2]*zToSource;
			   }
			   double phimiratan2_c = (180./Pi)*atan2((XYZmir_c[1]),(XYZmir_c[0]));
			   double XToSourcecamatan2=(180./Pi)*tan(thetamir_c*Pi/180.)*cos((phimiratan2_c + 180)*Pi/180.);
			   double YToSourcecamatan2=(180./Pi)*tan(thetamir_c*Pi/180.)*sin((phimiratan2_c + 180)*Pi/180.);

			   double XToSourcecamatan2_cm = (XToSourcecamatan2/(180./(475.*Pi)));
			   double YToSourcecamatan2_cm = (YToSourcecamatan2/(180./(475.*Pi)));*/
			tim_max = -inf;
			tim_min = inf;
			for(int phe_i = 0; phe_i < read32[3]; phe_i++) {
				bad_mirror = 0;
				phe[0] = 0;
				fread(&history, sizeof(int32_t), 1, f5);
				fread(&phe, sizeof(double), 2, f5);                 //time, ns and wavelength, nm
				if(phe[1] < 100 || phe[1] > 800) {
					cout << "Wrong wl: " << phe[1] << endl;
					break;
				}
				if(phe[0]*1E9 > tim_max) { tim_max = phe[0]*1E9;}
				if(phe[0]*1E9 < tim_min) { tim_min = phe[0]*1E9;}
				fread(&roco, sizeof(int32_t), 2, f5);
				mirrow = (int)floor((double)history / 256);
				mircol = mod(history,256) - 128;
				/*for(int i = 0; i < 5; i++){
				   if(hidmir[i][0] == mirrow && hidmir[i][1] == mircol){
				    //cout << "excluded\t" << phe_i << "\t" << history << "\t" << floor((double)history / 256) << "\t" << mod(history,256) - 128 << endl;
				    bad_mirror = 1;
				    miss_phe++;
				    break;
				   }
				   }*/
				//if(bad_mirror == 0){
				//cout << "\t\t" << phe_i << "\t" << history << "\t" << floor((double)history / 256) << "\t" << mod(history,256) - 128 << endl;
				map <type_key, int> :: iterator it;
				it = Mapa.find({roco[0], roco[1]});
				//cout << roco[0] << "\t" << roco[1] << "\t" << it->second << endl;
				if(it->second == 0 && roco[0] != -3 && roco[1] != 0) {}
				else{vector_phe_time[it->second].push_back(phe[0]*1E9);}
				//cout << event << "\t" << it->second << "\t" << phe[0]*10E9 << endl;
				//cout << "\t\t" << phe_i << "\t" << it->second << "\t" << phe[0]*10E9 << "\t" << phe[1] << "\t" << rowcol[it->second][3] << "\t" << rowcol[it->second][4] << endl;
				//cout << "\t\t" << it->second << endl;
				/*if(cluster != rowcol[i][2]){
				   cluster = rowcol[i][2];
				   }*/
				//  }
			}
			//cout << tim_max - tim_min << "\t" << tim_min << "\t" << tim_max << endl;
			//cout << vector_phe_time.size() << endl;;
			tim_min = tim_min - t_after_before;
			tim_max = tim_max + t_after_before;
			for(int i = 0; i < (int)vector_phe_time.size(); i++) {
				if(vector_phe_time[i].size() != 0 ) {
					vector <double> timesNSBphe = random_poisson(tim_min,tim_max);
					vector_phe_time[i].insert(vector_phe_time[i].end(), timesNSBphe.begin(), timesNSBphe.end());
					sort(vector_phe_time[i].begin(), vector_phe_time[i].end());
				}
				else{
					vector_phe_time[i] = random_poisson(tim_min,tim_max);
					sort(vector_phe_time[i].begin(), vector_phe_time[i].end());
				}
				for(int t = 0; t < (int)vector_phe_time[i].size(); t++) {
					vector_amp_MC[i].push_back(vector_amp_phe[rand()%vector_amp_phe.size()]);
					//vector_amp_MC[i].push_back(1.);
					//cout << i << "\t" << vector_amp_MC[i][t] << setw(4)<< "\t" << vector_phe_time[i][t] << "\t" << tim_min << "\t" << tim_max  << "\t" << vector_phe_time[i].size() << endl;
				}
			}
			trigger = 0;
			trig_time = -inf;
			trig_clust = -1;
			pix1 = -1;
			pix2 = -1;
			amp1 = -1;
			amp2 = -1;
			omp_set_num_threads(560);
			int i, ii;
			double t;
#pragma omp parallel shared(rowcol, vector_time10, vector_prelim_hold_time, vector_prelim_hold_pix, vector_phe_time, vector_amp_MC, t_sign, tim_min, tim_max, trig_amp) private(amp, t, ii, i)
			{
 #pragma omp for
				//fout1 << event << "\t" << (tim_max - tim_min) << endl;
				for(i = 0; i < (int)vector_phe_time.size(); i++) {                 //идем по сработавшим пикселям
					//cout << i << endl;
					if(vector_phe_time[i].size() != 0) {                 //если он не пустой:
						//grid_max = 0;
						//for(int ii = 0; ii < (int)vector_phe_time[i].size(); ii++){
						//cout << vector_phe_time[i][ii] << "\t" << vector_amp_MC[i][ii] << endl;
						for(t = tim_min; t <= tim_max; t=t+t_grid) {                 //идем по сетке от минимального в событии времени до максимального
							amp = 0;
							//#pragma omp for reduction(+: amp)
							for(ii = 0; ii < (int)vector_phe_time[i].size(); ii++) {                 //идем по конкретным значениям времён прихода фотоэлектронов
								if((t >= vector_phe_time[i][ii]) && (t <= (vector_phe_time[i][ii] + t_sign))) {                 //если время сетки идет после прихода фотоэлектрона(не посзже 30 ns):
									amp = amp + amp_t((t - vector_phe_time[i][ii]), vector_amp_MC[i][ii]);                 //считаем амплитуду в данном узле сетки от данного фэ и суммируем со вкладами от других фэ попавших в пиксель во время от t-30 до t
									//cout << "\t" << event << "\t" << i << "\t" << vector_amp_MC[i][ii] << "\t" << t - vector_phe_time[i][ii] << "\t" << amp << endl;
								}
								else if(vector_phe_time[i][ii] > t) {
									break;
								}
							}
							if(amp > trig_amp) {
						#pragma omp critical
								{
									vector_time10[i].push_back(t);
									//cout << i << "\t" <<  t << "\t" << amp << endl;
									vector_prelim_hold_time[(int)rowcol[i][2]].push_back(t);
									vector_prelim_hold_pix[(int)rowcol[i][2]].push_back(i);
								}
							}
						}
					}
				}
			}
			for(int i = 0; i < numb_of_clusters; i++) {
				sort(vector_prelim_hold_time[i].begin(), vector_prelim_hold_time[i].end());
				hold = -inf;
				for(int ii = 0; ii < vector_prelim_hold_time[i].size(); ii++) {
					if(vector_prelim_hold_time[i][ii] > (hold + t_hold)) {
						hold = vector_prelim_hold_time[i][ii];
						vector_integral_hold[i].push_back(hold);
						//cout << "\t" << hold[i] << endl;
					}
				}
				//if(event == ev){
				/*if(vector_integral_hold[i].size() > 0) {
				        for(int ii = 0; ii < vector_integral_hold[i].size(); ii++) {
				                cout << i << "\t" << vector_integral_hold[i][ii] << endl;
				        }
				   }*/
				//  }
			}
			for(int i = 0; i < (int)vector_phe_time.size(); i++) {
				if(vector_time10[i].size() != 0) {
					//cout << event << "\t" << i << "\t" << vector_time10[i].size() << endl;
					if(vector_integral_hold[(int)rowcol[i][2]].size() > 0) {
						for(int ii = 0; ii < neigh[i][2]; ii++) {
							//cout << i << "\tneigh\t" << neigh[i][3+ii]  << endl;
							if(neigh[i][1] == 1 && neigh[neigh[i][3+ii]][1] == 1) {
								//if(vector_integral_hold[(int)rowcol[neigh[i][3+ii]][2]].size() > 0){
								if(vector_time10[neigh[i][3+ii]].size() != 0) {
									for(int j = 0; j < vector_time10[i].size(); j++) {
										for(int jj = 0; jj < vector_time10[neigh[i][3+ii]].size(); jj++) {
											//cout << event << "\t" << i << "\t" << neigh[i][3+ii] << "\t" << vector_time10[i][j] << "\t" << vector_time10[neigh[i][3+ii]][jj] << endl;
											if(((vector_time10[i][j] - vector_time10[neigh[i][3+ii]][jj]) <= trig_window) && ((vector_time10[i][j] - vector_time10[neigh[i][3+ii]][jj]) >= 0)) {
												for(int ij = 0; ij < vector_integral_hold[(int)rowcol[i][2]].size(); ij++) {                                                                                                                                                                                                            //вероятно считать нужно по холдам neigh[i][3+ii] а не i, поскольку сосед сработал раньше. Однако холды соседа и i должны быть одними и теми же, ведь это один и тот-же кластер, так что можно и так оставить (или можно поменять vector_integral_hold[(int)rowcol[i][2]] на vector_integral_hold[(int)rowcol[neigh[i][3+ii]][2]]) но это один кластер, у них общие холды)
													//cout << vector_time10[i][j] << "\t" << vector_integral_hold[(int)rowcol[i][2]][ij] << endl;
													if(abs(vector_time10[neigh[i][3+ii]][jj] - vector_integral_hold[(int)rowcol[neigh[i][3+ii]][2]][ij]) < 0.1) {                 //а вот тут холды соседа сравниваются с превышениями порога этим же соседом. Раз сосед сработал первым - все правильно
														//cout << "trig" << endl;
														trigger = 1;
														trig_time = vector_time10[neigh[i][3+ii]][jj];                 //вероятно время триггера - это время первого превысевшего порог пикселя, так что тут должно быть время соседа
														trig_clust = rowcol[neigh[i][3+ii]][2];
														pix2 = i;                 //вероятно везеде нужно поменять местами значения 1 и 2, ведь первым был сосед
														pix1 = neigh[i][3+ii];
														amp2 = vector_time10[i][j];
														amp1 = vector_time10[neigh[i][3+ii]][jj];
														goto exit;
													}
												}
											}
											else if(((vector_time10[i][j] - vector_time10[neigh[i][3+ii]][jj]) >= -trig_window) && ((vector_time10[i][j] - vector_time10[neigh[i][3+ii]][jj]) <= 0)) {                 //тут первым срабатывает i, тк он меньше
												//cout << vector_time10[i][j] - vector_time10[neigh[i][3+ii]][jj] << endl;
												for(int ij = 0; ij < vector_integral_hold[(int)rowcol[neigh[i][3+ii]][2]].size(); ij++) {                 //нужно идти по холдам кластера с i а не его соседа, но они в одном кластере, так что норм
													//cout << vector_time10[neigh[i][3+ii]][jj] << "\t" << vector_integral_hold[(int)rowcol[neigh[i][3+ii]][2]][ij] << endl;
													if(abs(vector_time10[i][j] - vector_integral_hold[(int)rowcol[i][2]][ij]) < 0.1) {                 //тут и дальше - правильно
														//cout << "trig" << endl;
														trigger = 1;
														trig_time = vector_time10[i][j];
														trig_clust = rowcol[i][2];
														pix1 = i;
														pix2 = neigh[i][3+ii];
														amp1 = vector_time10[i][j];
														amp2 = vector_time10[neigh[i][3+ii]][jj];
														goto exit;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
exit:                   ;
			if(trigger == 1) {
				//cout << event  << "\t" << "trig" << endl;
				vector <double> vector_event_amp;
				vector <int> vector_event_numb;
				vector <double> vector_event_amp1;
				vector <int> vector_event_numb1;
				for(int i = 0; i < (int)vector_phe_time.size(); i++) {
					if(vector_phe_time[i].size() != 0) {
						if(vector_integral_hold[(int)rowcol[i][2]].size() > 0) {
							for(int ii = 0; ii < vector_integral_hold[(int)rowcol[i][2]].size(); ii++) {
								if(vector_integral_hold[(int)rowcol[i][2]][ii] >= trig_time - t_hold/2 && vector_integral_hold[(int)rowcol[i][2]][ii] <= trig_time + t_hold/2) {                 //т.к. через 80нс после триггера происходит считывание тех кластеров, у которых на этот момент есть холд
									amp = 0;
									amp1 = 0;
									for(int j = 0; j < (int)vector_phe_time[i].size(); j++) {                 //идем по конкретным значениям времён прихода фотоэлектронов
										if(vector_phe_time[i][j] > vector_integral_hold[(int)rowcol[i][2]][ii] - 600 && vector_phe_time[i][j] <= vector_integral_hold[(int)rowcol[i][2]][ii] + integrate_window) {
											amp = amp + vector_amp_MC[i][j]*inter.Eval(vector_integral_hold[(int)rowcol[i][2]][ii] + integrate_window - vector_phe_time[i][j]);
											//cout << i << "\t" << vector_phe_time[i][j] << "\t" << vector_integral_hold[(int)rowcol[i][2]][ii]+35 << "\t" <<  vector_integral_hold[(int)rowcol[i][2]][ii]+35 - vector_phe_time[i][j] << "\tamp:\t" << inter.Eval(vector_integral_hold[(int)rowcol[i][2]][ii] + 35 - vector_phe_time[i][j]) << "\t" << amp << endl;
										}
										if(vector_phe_time[i][j] >= vector_integral_hold[(int)rowcol[i][2]][ii] && vector_phe_time[i][j] <= (vector_integral_hold[(int)rowcol[i][2]][ii] + integrate_window)) {
											amp1 = amp1 + vector_amp_MC[i][j];
										}
									}
									vector_event_amp.push_back(amp);
									vector_event_numb.push_back(i);

									vector_event_amp1.push_back(amp1);
									vector_event_numb1.push_back(i);
									//cout << read32[0] << "\t" << read32[1] << "\t" << i << endl;
									break;
								}
							}
						}
					}
				}
				if(vector_event_numb.size() > 0) {
					cout << evch << "\t" << read32[2] << "\t" << read32[0] << "\t" << read32[1] << "\t" << vector_event_amp.size() << "\t" << readdubble[0]/1E12 << "\t" << (readdubble[10] + readdubble[13])/1E3 << "\t " << (readdubble[11] + readdubble[14])/1E3 << endl;
					//fout << read32[0] << "\t" << read32[1] << "\t" << vector_event_amp.size() << "\t" << readdubble[0]/1E12 << "\t" << XToSourcecamatan2_cm << "\t" << YToSourcecamatan2_cm << "\t" << (readdubble[10] + readdubble[13])/1E3 << "\t" << (readdubble[11] + readdubble[14])/1E3 << endl;
					fout << read32[0] << "\t" << read32[1] << "\t" << vector_event_amp.size() << "\t" << readdubble[0]/1E12 << "\t" << (readdubble[10] + readdubble[13])/1E3 << "\t" << (readdubble[11] + readdubble[14])/1E3 << "\t" << tet_center << "\t" << fi_center << "\t" << tet_source << "\t" << fi_source << "\t" << readdubble[8] << endl;
					//fout1 << read32[0] << "\t" << read32[1] << "\t" << trig_clust+1  << "\t" << vector_event_amp.size() << "\t" << readdubble[0]/1E12 << "\t" << R << "\t"<< readdubble[17]*180/Pi << endl;
					for(int i = 0; i < (int)vector_event_numb.size(); i++) {
						//cout << neigh_clean[vector_event_numb[i]][1] << "\t" << neigh_clean[vector_event_numb[i]][2] << "\t" << neigh_clean[vector_event_numb[i]][3] << "\t"<< neigh_clean[vector_event_numb[i]][4] << "\t" << vector_event_amp[i] << endl;
						fout << real_numb[vector_event_numb[i]][0] << "\t" << real_numb[vector_event_numb[i]][1] << "\t" << neigh_clean[vector_event_numb[i]][3] << "\t"<< neigh_clean[vector_event_numb[i]][4] << "\t" << vector_event_amp[i] << endl;
						//fout1 << real_numb[vector_event_numb1[i]][1] << "\t" << real_numb[vector_event_numb1[i]][2] << "\t" << neigh_clean[vector_event_numb1[i]][3] << "\t"<< neigh_clean[vector_event_numb1[i]][4] << "\t" << vector_event_amp1[i] << endl;
					}
				}
				for(int i = 0; i < (int)vector_event_numb.size(); i++) {
					vector_event_amp.clear();
					vector_event_numb.clear();
					vector_event_amp1.clear();
					vector_event_numb1.clear();
				}
				evch++;
			}
			for(int i = 0; i < (int)vector_phe_time.size(); i++) {
				vector_amp_MC[i].clear();
				vector_phe_time[i].clear();
				vector_amp[i].clear();
				vector_time10[i].clear();
				//test
				vector_integral35[i].clear();
			}
			for(int i = 0; i < (int)numb_of_clusters; i++) {
				vector_integral_hold[i].clear();
				vector_prelim_hold_time[i].clear();
				vector_prelim_hold_pix.clear();
			}
			//cout << "miss_phe " << miss_phe << endl;
			event++;
			/*if(event == ev + 2){
			   return 0;
			   }*/
		}
		else{
			fread(&readdubble, sizeof(double), 20, f5);
			for(int phe_i = 0; phe_i < read32[3]; phe_i++) {
				fread(&history, sizeof(int32_t), 1, f5);
				fread(&phe, sizeof(double), 2, f5);         //time, ns and wavelength, nm
				fread(&roco, sizeof(int32_t), 2, f5);
			}
		}
	}
	return 0;
}
