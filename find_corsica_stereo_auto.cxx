#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;
string line;
char path[100];
int numb, iscat, max_event, stereo, clust, number_of_pixels, ster[4]={0}, tel_ster2[4] = {0}, edge;
double energy, angle, size_e, Xc, Yc, con2, x_max_coord, y_max_coord, length, width, dist0, dist1, azwidth, miss, alpha0, alpha1, a_axis, b_axis, distance_t, x_cam, y_cam, x_ground, y_ground, tet_center, fi_center, tet_source, fi_source, xmax;
int gam = 100, edge_cut = -1;
char str_tel[200];
double edge1 = 0.0, edge2 = 0.0, size_cut = 100;
vector <int> iacts;
string input_file_name, out_path, out_add;

void read_param(){
	string word;
	int iact;
	ifstream file1("./clean.param");
	if (!file1.is_open()) {
		cout << "file neighbours is not found" << endl;
		exit(1);
	}
	getline(file1, line);
	istringstream ist(line);
	cout << line << endl;
	ist >> word >> edge1 >> edge2;
	getline(file1, line);
	istringstream ist1(line);
	ist1 >> word;
	while(ist1) {
		ist1 >> iact;
		if(ist1) {
			iacts.push_back(iact-1);
		}
	}
	getline(file1, line);
	istringstream ist2(line);
	ist2 >> word >> input_file_name;
	getline(file1, line);
	istringstream ist3(line);
	ist3 >> word >> out_path;
	getline(file1, line);
	istringstream ist4(line);
	ist4 >> word >> gam;
	if(gam != 0 && gam != 1) {
		cout << "invalid gamma parameter" << endl;
		exit(1);
	}
	getline(file1, line);
	istringstream ist5(line);
	ist5 >> word >> out_add;
	getline(file1, line);
	istringstream ist6(line);
	ist6 >> word >> edge_cut;
	getline(file1, line);
	istringstream ist7(line);
	ist7 >> word >> size_cut;
	file1.close();
}

int main()
{
	read_param();
	vector<vector<double> > vector_size(iacts.size());
	vector<vector<int> > event_number(iacts.size());
	vector<vector<int> > event_iscat(iacts.size());
	vector<vector<int> > vector_pix(iacts.size());
	vector<vector<double> > vector_energy(iacts.size());
	vector<vector<double> > vector_x_ground(iacts.size());
	vector<vector<double> > vector_y_ground(iacts.size());
	vector<vector<double> > vector_x_cam(iacts.size());
	vector<vector<double> > vector_y_cam(iacts.size());
	vector<vector<double> > vector_xc(iacts.size());
	vector<vector<double> > vector_yc(iacts.size());
	vector<vector<double> > vector_dist(iacts.size());
	vector<vector<double> > vector_edge(iacts.size());
	vector<vector<double> > vector_a_axis(iacts.size());
	vector<vector<double> > vector_b_axis(iacts.size());
	vector<vector<double> > vector_tet_center(iacts.size());
	vector<vector<double> > vector_fi_center(iacts.size());
	vector<vector<double> > vector_tet_source(iacts.size());
	vector<vector<double> > vector_fi_source(iacts.size());
	vector<vector<double> > vector_alpha(iacts.size());
	vector<vector<double> > vector_width(iacts.size());
	vector<vector<double> > vector_length(iacts.size());
	vector<vector<double> > vector_xmax(iacts.size());
	vector<vector<double> > vector_x_max_coord(iacts.size());
	vector<vector<double> > vector_y_max_coord(iacts.size());
	sprintf(str_tel, "%s%s_stereo_%dtel_%s_edge%d_size%.0f.txt", out_path.c_str(), input_file_name.c_str(), iacts.size(), out_add.c_str(), edge_cut, size_cut);
	cout << str_tel << endl;
	ofstream fout(str_tel);
	for(int tel = 0; tel < iacts.size(); tel++) {
		sprintf(str_tel, "%s%s_hillas_iact%d_%s.txt", out_path.c_str(), input_file_name.c_str(), iacts[tel]+1, out_add.c_str());
		ifstream file(str_tel);
		if (!file.is_open()) {
			cout << "Файл не найден" << endl;
		}
		else {
			while(!file.eof()) {
				getline(file, line);
				stringstream ist(line);
				ist >> numb >> numb >> iscat >> number_of_pixels >> energy >> size_e >> Xc >> Yc >> con2 >> length >> width >> dist0 >> dist1 >> azwidth >> miss >> alpha0 >> alpha1 >> a_axis >> b_axis >> x_max_coord >> y_max_coord >> edge >> x_cam >> y_cam >> x_ground >> y_ground >> tet_center >> fi_center >> tet_source >> fi_source >> xmax;
				if((edge == edge_cut || edge == (edge_cut - 1)) && size_e > size_cut) {
					event_number[tel].push_back(numb);
					event_iscat[tel].push_back(iscat);
					vector_pix[tel].push_back(number_of_pixels);
					vector_energy[tel].push_back(energy);
					vector_x_ground[tel].push_back(x_ground);
					vector_y_ground[tel].push_back(y_ground);

					//gamma
					if(gam == 1) {
						vector_x_cam[tel].push_back(x_cam);
						vector_y_cam[tel].push_back(y_cam);
						vector_tet_source[tel].push_back(tet_source);
						vector_fi_source[tel].push_back(fi_source);
					}
					else{
						//hadrons:
						vector_x_cam[tel].push_back(0);
						vector_y_cam[tel].push_back(0);
						vector_tet_source[tel].push_back(tet_center);
						vector_fi_source[tel].push_back(fi_center);
					}
					vector_xc[tel].push_back(Xc);
					vector_yc[tel].push_back(Yc);
					vector_size[tel].push_back(size_e);
					vector_dist[tel].push_back(dist1);
					vector_edge[tel].push_back(edge);
					vector_a_axis[tel].push_back(a_axis);
					vector_b_axis[tel].push_back(b_axis);
					vector_tet_center[tel].push_back(tet_center);
					vector_fi_center[tel].push_back(fi_center);
					vector_alpha[tel].push_back(alpha1);
					vector_width[tel].push_back(width);
					vector_length[tel].push_back(length);
					vector_xmax[tel].push_back(xmax);
					vector_x_max_coord[tel].push_back(x_max_coord);
					vector_y_max_coord[tel].push_back(y_max_coord);
					cout << number_of_pixels << "\t" << energy << "\t" << sqrt(pow(x_ground,2) + pow(y_ground,2)) << "\t" << size_e << "\t" << dist1  << "\t" << edge << endl;
				}
			}
		}
	}
	max_event = event_number[0][event_number[0].size()-1] + 1000; //первое событие 1 телескопа + 1000 событий(вдруг другие еще 1к настреляли)
	double init[4] = {0};
	//max_event = max_element(event_number[0].begin(), event_number[0].end());
	for(int i = 0; i < max_event; i++) {
		for(int j = 0; j < 10; j++) {
			stereo = 0;
			vector<vector <double> > telescopes(16);
			for(int tel = 0; tel < iacts.size(); tel++) {
				for(int ii = init[tel]; ii < event_number[tel].size(); ii++) {
					if(event_number[tel][ii] == i && event_iscat[tel][ii] == j) {
						stereo++;
						energy = vector_energy[tel][ii];
						tet_center = vector_tet_center[tel][ii];
						fi_center = vector_fi_center[tel][ii];
						tet_source = vector_tet_source[tel][ii];
						fi_source = vector_fi_source[tel][ii];
						x_cam = vector_x_cam[tel][ii];
						y_cam = vector_y_cam[tel][ii];
						xmax = vector_xmax[tel][ii];
						telescopes[0].push_back(iacts[tel]);
						telescopes[1].push_back(vector_pix[tel][ii]);
						telescopes[2].push_back(vector_size[tel][ii]);
						telescopes[3].push_back(vector_x_ground[tel][ii]);
						telescopes[4].push_back(vector_y_ground[tel][ii]);
						telescopes[5].push_back(vector_dist[tel][ii]);
						telescopes[6].push_back(vector_edge[tel][ii]);
						telescopes[7].push_back(vector_a_axis[tel][ii]);
						telescopes[8].push_back(vector_b_axis[tel][ii]);
						telescopes[9].push_back(vector_alpha[tel][ii]);
						telescopes[10].push_back(vector_width[tel][ii]);
						telescopes[11].push_back(vector_xc[tel][ii]);
						telescopes[12].push_back(vector_yc[tel][ii]);
						telescopes[13].push_back(vector_length[tel][ii]);
						telescopes[14].push_back(vector_x_max_coord[tel][ii]);
						telescopes[15].push_back(vector_y_max_coord[tel][ii]);
						init[tel] = ii;
						goto next_telescope;
					}
					if(event_number[tel][ii] > i) {
						goto next_telescope;
					}
				}
next_telescope:                 ;
			}
			if(stereo > 0 && energy > 0) {
				cout << i << "\t" << j << "\t" << telescopes[0].size() << "\t" << energy << "\t" << tet_center << "\t" << fi_center << "\t" << tet_source << "\t" << fi_source << "\t" << x_cam << "\t" << y_cam << "\t" << xmax << endl;
				fout << i << "\t" << j << "\t" << telescopes[0].size() << "\t" << energy << "\t" << tet_center << "\t" << fi_center << "\t" << tet_source << "\t" << fi_source << "\t" << x_cam << "\t" << y_cam << "\t" << xmax << endl;
				//cout << telescopes[0].size() << "\t" << telescopes.size() << endl;
				for(int jj = 0; jj < telescopes[0].size(); jj++) {
					for(int jjj = 0; jjj < telescopes.size(); jjj++) {
						cout << telescopes[jjj][jj] << "\t";
						fout << telescopes[jjj][jj] << "\t";
					}
					cout << endl;
					fout << endl;
				}
			}
			for(int jj = 0; jj < telescopes.size(); jj++) {
				telescopes[jj].clear();
			}
		}
	}
	cout << max_event << endl;
	return 0;
}
