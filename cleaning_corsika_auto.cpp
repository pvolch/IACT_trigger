#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <locale>
#include <cstdlib>
#define Pi 3.141592653589793238462643
using namespace std;
double x_pos[64][23], y_pos[64][23], bmp[64][23], sig[64][23] = {1};
int j, jj, jjj, f;
int kkk[6][64][23], pos[6][64][23];
int number_of_pixels, numb_clean = 0, gam = 100;
double edge1 = 0.0, edge2 = 0.0;
vector <int> iacts;
string line, input_file_name, out_path, out_add;
vector<vector<double> > edge_pix(2);
double tet_center, fi_center, tet_source, fi_source;

double * get_hillas(vector<vector<double> > vector_pixel, double x_cam, double y_cam, int gam){
	static double hillas[15];
	double Xc[3],Yc[3],Xc2[3],Yc2[3],XYc[3],sigx2[3], sigy2[3], sigxy[3], d[3],z[3],length[3],width[3],azwidth[3],U[3],V[3],bm[3],miss[3],dist[3],alpha[3], con2, con1, x_max_c[3], y_max_c[3], x_max_coord, y_max_coord, amp_max, amp2, event_size, a_axis[3], b_axis[3];
	//static double hillas[18]; //size, Xc[0],Yc[0], con2, length[0], width[0], dist[0], dist[1], dist[2], azwidth[0],
	//azwidth[1], azwidth[2], miss[0], miss[1], miss[2], alpha[0], alpha[1], alpha[2]
	//cout << "get hillas\t" << vector_pixel[0].size() << endl;
	amp2 = 0;
	amp_max = 0;
	x_max_coord = 0;
	y_max_coord = 0;
	for(int i = 0; i < 3; i++) {
	    x_max_c[i] = 0;
	    y_max_c[i] = 0;
	}
	con2 = 0;
	con1 = 0;
	event_size = 0;
	for(int i = 0; i < 2; i++) {
		Xc[i] = 0;
		Yc[i] = 0;
		Xc2[i]=0;
		Yc2[i]=0;
		XYc[i]=0;
		sigx2[i] = 0;
		sigy2[i] = 0;
		sigxy[i] = 0;
		d[i]=0;
		z[i]=0;
		length[i]=0;
		width[i] = 0;
		azwidth[i]=0;
		U[i]=0;
		V[i]=0;
		bm[i]=0;
		miss[i]=0;
		dist[i]=0;
		alpha[i]=0;
		a_axis[i] = 0;
		b_axis[i] = 0;
	}
	for(int k = 0; k < vector_pixel[4].size(); k++) {
		Xc[0] += (vector_pixel[4][k]*(vector_pixel[2][k]));
		Yc[0] += (vector_pixel[4][k]*(vector_pixel[3][k]));
		Xc2[0] += (vector_pixel[4][k]*pow(vector_pixel[2][k], 2));
		Yc2[0] += (vector_pixel[4][k]*pow(vector_pixel[3][k], 2));
		XYc[0] += (vector_pixel[4][k]*(vector_pixel[2][k])*(vector_pixel[3][k]));

		Xc[1] += (vector_pixel[4][k]*(vector_pixel[2][k] - x_cam));
		Yc[1] += (vector_pixel[4][k]*(vector_pixel[3][k] - y_cam));
		Xc2[1] += (vector_pixel[4][k]*pow(vector_pixel[2][k] - x_cam, 2));
		Yc2[1] += (vector_pixel[4][k]*pow(vector_pixel[3][k] - y_cam, 2));
		XYc[1] += (vector_pixel[4][k]*(vector_pixel[2][k] - x_cam)*(vector_pixel[3][k] - y_cam));

		if(amp_max <= vector_pixel[4][k]) {
			amp2 = amp_max;
			amp_max = vector_pixel[4][k];
			x_max_c[0] = x_max_c[1];
			y_max_c[0] = y_max_c[1];
			x_max_c[1] = x_max_c[2];
			y_max_c[1] = y_max_c[2];
			x_max_c[2] = vector_pixel[2][k];
			y_max_c[2] = vector_pixel[3][k];
		}
		else if(amp_max > vector_pixel[4][k] && vector_pixel[4][k] > amp2) {
			amp2 = vector_pixel[4][k];
		}

		event_size += vector_pixel[4][k];
	}
	
	x_max_coord = x_max_c[2];
	y_max_coord = y_max_c[2];
	con2 = (amp_max + amp2)/event_size;
	con1 = amp_max/event_size;

	for(int i = 0; i < 2; i++) {

		Xc[i] = Xc[i]/event_size;
		Yc[i] = Yc[i]/event_size;
		Xc2[i]=Xc2[i]/event_size;
		Yc2[i]=Yc2[i]/event_size;
		XYc[i]=XYc[i]/event_size;
		sigx2[i] = Xc2[i] - pow(Xc[i],2);
		sigy2[i] = Yc2[i] - pow(Yc[i],2);
		sigxy[i] = XYc[i] - (Xc[i])*(Yc[i]);
		d[i]=sigy2[i]-sigx2[i];
		a_axis[i] = (d[i] + sqrt(pow(d[i],2) + 4*pow(sigxy[i], 2)))/(2*sigxy[i]);
		b_axis[i] = Yc[i] - (a_axis[i]*Xc[i]);
		z[i]=sqrt(pow(d[i],2)+4.*pow(sigxy[i],2));
		length[i]=sqrt((sigx2[i] + sigy2[i] + z[i])/2.);
		if(sigx2[i] + sigy2[i] - z[i] >= 0) {
			width[i]=sqrt((sigx2[i] + sigy2[i] - z[i])/2.);
		}
		else{
			width[i] = 0;
		}
		dist[i] = sqrt(pow(Xc[i],2)+pow(Yc[i],2));
		azwidth[i] = sqrt((pow(Xc[i],2)*Yc2[i]-2.*Xc[i]*Yc[i]*XYc[i]+Xc2[i]*pow(Yc[i],2))/pow(dist[i],2));
		if(z[i] != 0) {
			U[i]=(1. + (d[i]/z[i]));
			V[i]=2. - U[i];
			bm[i] = (0.5*(U[i]*pow(Xc[i],2)+V[i]*pow(Yc[i],2)))-(2.*sigxy[i]*Xc[i]*Yc[i]/z[i]);
		}
		else{
			U[i] = nan("");
			bm[i] = nan("");
			V[i] = nan("");
			miss[i] = nan("");
			alpha[i] = nan("");
		}
		if(bm[i] < 0) {
			bm[i]=nan("");
			miss[i]=nan("");
			alpha[i]=nan("");
		}
		else if(bm[i] >= 0) {
			miss[i]=sqrt(bm[i]);
			alpha[i]=asin(miss[i]/dist[i])*(180./Pi);
		}
	}
//gamma
	if(gam == 1) {
		hillas[0] = event_size;
		hillas[1] = 0.1206 * Xc[0];
		hillas[2] = 0.1206*Yc[0];
		hillas[3] = con2;
		hillas[4] = 0.1206*length[0];
		hillas[5] = 0.1206*width[0];
		hillas[6] = 0.1206*dist[0];
		hillas[7] = 0.1206*dist[1];
		hillas[8] = 0.1206*azwidth[1];
		hillas[9] =  0.1206*miss[1];
		hillas[10] = alpha[0];
		hillas[11] = alpha[1];
		hillas[12] = a_axis[0];
		hillas[13] = b_axis[0];
		hillas[14] = 0.1206*x_max_coord;
		hillas[15] = 0.1206*y_max_coord;
	}
//hadron
	else{
		hillas[0] = event_size;
		hillas[1] = 0.1206 * Xc[0];
		hillas[2] = 0.1206*Yc[0];
		hillas[3] = con2;
		hillas[4] = 0.1206*length[0];
		hillas[5] = 0.1206*width[0];
		hillas[6] = 0.1206*dist[0];
		hillas[7] = 0.1206*dist[0];
		hillas[8] = 0.1206*azwidth[0];
		hillas[9] =  0.1206*miss[0];
		hillas[10] = alpha[0];
		hillas[11] = alpha[0];
		hillas[12] = a_axis[0];
		hillas[13] = b_axis[0];
		hillas[14] = 0.1206*x_max_coord;
		hillas[15] = 0.1206*y_max_coord;
	}
	return hillas;
}

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
			iacts.push_back(iact);
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
	file1.close();
}

void read_neigbours(int iact_i){
	for (int count = 0; count < 23; count++)
	{
		for (int coun = 0; coun < 64; coun++)
		{
			for (int cou = 0; cou < 6; cou++)
			{
				kkk[cou][coun][count] = -1;
				pos[cou][coun][count] = -1;
			}
		}
	}
	char str_line[100];
	if(iact_i == 1) {
		sprintf(str_line, "./neighbour_clean_experiment_iact1_375.txt");
	}
	else{
		sprintf(str_line, "./neighbour_clean_experiment_iact2_375.txt");
	}
	cout << iact_i << "\t" << str_line << endl;
	ifstream file1(str_line);
	if (!file1.is_open()) {
		cout << "file neighbours is not found" << endl;
		exit(1);
	}
	int q = 0;
	while(!file1.eof()) {
		int k[6]= {0},sos[6]= {0},kk=0,Nsos=0, ii=0, kx = 0,ix = 0;
		double x=0,y=0;
		getline(file1, line);
		if(file1.eof())
		{
			cout << "file neighbours size: " << q << endl;
			break;
		}
		q++;
		istringstream ist(line);
		ist >> kk >> ii;
		ist >> x_pos[ii][kk-1] >> y_pos[ii][kk-1] >> Nsos;
		//cout << kk << "\t" << ii << "\t" << Nsos <<"\t\t";
		if(Nsos < 6) {
			edge_pix[0].push_back(kk);
			edge_pix[1].push_back(ii);
		}
		//int zi = 0;
		/*for(int i = 0; i < 28; i++){
		   //cout << kk <<" " << exclud_clust[i] <<"\t" << ii << " " << exclud_numb[i] << endl;
		   if(kk == exclud_clust[i] && ii == exclud_numb[i]){
		    for(int j = 0; j < 6; j++){
		      kx = 0;
		      ix = 0;
		      ist >> kx >> ix;
		      kkk[j][ii][kk] = -1;
		      pos[j][ii][kk] = -1;
		    }
		    zi = 1;
		    break;
		   }
		   }*/
		//if((kk == exclud_clust[0] && ii == exclud_numb[0]) || (kk == exclud_clust[1] && ii == exclud_numb[1]) || (kk == exclud_clust[2] && ii == exclud_numb[2])){
		//if(zi == 0){
		for(int j = 0; j < Nsos; j++) {
			kx = 0;
			ix = 0;
			ist >> kx >> ix;
			//if((kx == exclud_clust[0] &&  ix == exclud_numb[0]) || (kx == exclud_clust[1] &&  ix == exclud_numb[1]) || (kx == exclud_clust[2] &&  ix == exclud_numb[2])){
			/*for(int i = 0; i < 28; i++){
			   if(kk == exclud_clust[i] && ii == exclud_numb[i]){
			    kkk[j][ii][kk] = -1;
			    pos[j][ii][kk] = -1;
			    zi = 0;
			    break;
			   }
			   }*/
			//if(zi == 0){
			kkk[j][ii][kk-1] = kx - 1;
			pos[j][ii][kk-1] = ix;
			//}
			//cout << kkk[j][ii][kk-1] << " " << pos[j][ii][kk-1] << " ";
		}
		//  }
		//cout << endl;
	}
	file1.close();
}

vector<vector<double> > cleaning(double bmp[64][23]){
	for(f = 0; f < 22; f++)
	{
		jj = 0;
		jjj = 0;
		for (int sc = 0; sc < 64; sc++)
		{
			if(bmp[sc][f] > 0) {
				jj++;
			}
			if(bmp[sc][f] <= 0) {
				bmp[sc][f] = 0;
				jjj++;
			}
		}
		if( jj > 0) {
			for (int sc = 0; sc < 64; sc = sc + 2)
			{
				//cout << sc << "\t" << f << "\t" << bmp[sc][f] << "\t" << pos[1][sc][f] << "\t" << kkk[1][sc][f] << "\t" << bmp[pos[1][sc][f]][kkk[1][sc][f]] << endl;
				if(bmp[sc][f]>edge1*sig[sc][f])
				{
					//cout << pos[0][sc][f] << "\t" << kkk[0][sc][f] << endl;
					if((bmp[pos[0][sc][f]][kkk[0][sc][f]]>edge2*sig[pos[0][sc][f]][kkk[0][sc][f]] && bmp[pos[0][sc][f]][kkk[0][sc][f]]>0) ||
					   (bmp[pos[1][sc][f]][kkk[1][sc][f]]>edge2*sig[pos[1][sc][f]][kkk[1][sc][f]] && bmp[pos[1][sc][f]][kkk[1][sc][f]]>0) ||
					   (bmp[pos[2][sc][f]][kkk[2][sc][f]]>edge2*sig[pos[2][sc][f]][kkk[2][sc][f]] && bmp[pos[2][sc][f]][kkk[2][sc][f]]>0) ||
					   (bmp[pos[3][sc][f]][kkk[3][sc][f]]>edge2*sig[pos[3][sc][f]][kkk[3][sc][f]] && bmp[pos[3][sc][f]][kkk[3][sc][f]]>0) ||
					   (bmp[pos[4][sc][f]][kkk[4][sc][f]]>edge2*sig[pos[4][sc][f]][kkk[4][sc][f]] && bmp[pos[4][sc][f]][kkk[4][sc][f]]>0) ||
					   (bmp[pos[5][sc][f]][kkk[5][sc][f]]>edge2*sig[pos[5][sc][f]][kkk[5][sc][f]] && bmp[pos[5][sc][f]][kkk[5][sc][f]]>0))
					{
						//  cout << "11111" << "\t" << i << "\t" << f << endl;
						// cout << bmp[i][f] << endl;
					}
					else
					{
						bmp[sc][f] = 0;
					}
				}
				else if(bmp[sc][f]>edge2*sig[sc][f])
				{
					if((bmp[pos[0][sc][f]][kkk[0][sc][f]]>edge1*sig[pos[0][sc][f]][kkk[0][sc][f]] && bmp[pos[0][sc][f]][kkk[0][sc][f]]>0) ||
					   (bmp[pos[1][sc][f]][kkk[1][sc][f]]>edge1*sig[pos[1][sc][f]][kkk[1][sc][f]] && bmp[pos[1][sc][f]][kkk[1][sc][f]]>0) ||
					   (bmp[pos[2][sc][f]][kkk[2][sc][f]]>edge1*sig[pos[2][sc][f]][kkk[2][sc][f]] && bmp[pos[2][sc][f]][kkk[2][sc][f]]>0) ||
					   (bmp[pos[3][sc][f]][kkk[3][sc][f]]>edge1*sig[pos[3][sc][f]][kkk[3][sc][f]] && bmp[pos[3][sc][f]][kkk[3][sc][f]]>0) ||
					   (bmp[pos[4][sc][f]][kkk[4][sc][f]]>edge1*sig[pos[4][sc][f]][kkk[4][sc][f]] && bmp[pos[4][sc][f]][kkk[4][sc][f]]>0) ||
					   (bmp[pos[5][sc][f]][kkk[5][sc][f]]>edge1*sig[pos[5][sc][f]][kkk[5][sc][f]] && bmp[pos[5][sc][f]][kkk[5][sc][f]]>0))
					{
						// cout << "22222" << "\t" << i << "\t" << f << endl;
						//cout << bmp[sc][f] << endl;
					}
					else
					{
						bmp[sc][f] = 0;
					}
				}
				else
				{
					bmp[sc][f] = 0;
				}
			}
			j = 0;
			for (int sc = 0; sc < 64; sc++)
			{
				if(bmp[sc][f] == 0)
				{
					j++;
				}
			}
		}
	}
	vector<vector<double> > vector_pixel( 5, vector<double> (0));
	for (int count = 0; count < 23; count++)
	{
		for (int coun = 0; coun < 64; coun+=2)
		{
			if(bmp[coun][count] != 0) {
				number_of_pixels++;
				vector_pixel[0].push_back(count+1);
				vector_pixel[1].push_back(coun);
				vector_pixel[4].push_back(bmp[coun][count]);
				vector_pixel[2].push_back(x_pos[coun][count]);
				vector_pixel[3].push_back(y_pos[coun][count]);
			}
		}
	}
	return vector_pixel;
}

int main(){
	read_param();
	for(int n = 0; n < iacts.size(); n++) {
		numb_clean = 0;
		read_neigbours(iacts[n]);
		char str_line[200];
		sprintf(str_line, "%s%s_iact%d.txt", out_path.c_str(), input_file_name.c_str(), iacts[n]);
		ifstream file(str_line);
		sprintf(str_line, "%s%s_clean_iact%d_%s.txt", out_path.c_str(), input_file_name.c_str(), iacts[n], out_add.c_str());
		ofstream fout(str_line);
		sprintf(str_line, "%s%s_hillas_iact%d_%s.txt", out_path.c_str(), input_file_name.c_str(), iacts[n], out_add.c_str());
		ofstream fout1(str_line);
		if (!file.is_open()) {
			cout << "data файл не найден" << endl;
		}
		else {
			for (int count = 0; count < 23; count++)
			{
				for (int coun = 0; coun < 64; coun++)
				{
					sig[coun][count] = 1;
				}
			}
			while(!file.eof()) {
				getline(file, line);
				//cout << line << endl;
				if(!file.eof()) {
					int numb = 0, i_scat = 0, number_of_pixels = 0, edge = 0;
					double energy = 0, distance = 0, angle = 0, x_cam = 0, y_cam = 0, x_ground = 0, y_ground = 0, xmax = 0;
					for (int count = 0; count < 23; count++)
					{
						for (int coun = 0; coun < 64; coun++)
						{
							y_pos[coun][count]=0;
							x_pos[coun][count]=0;
							bmp[coun][count]=0;
						}
					}
					stringstream ist(line);
					//ist >> numb >> i_scat >> clust >> number_of_pixels >> energy >> distance >> angle;
					ist >> numb >> i_scat >> number_of_pixels >> energy >> x_ground >> y_ground >> tet_center >> fi_center >> tet_source >> fi_source >> xmax;

					double xToSource=1*sin(tet_source)*cos(fi_source);
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
					double YToSourcecamatan2_cm = (YToSourcecamatan2/(180./(475.*Pi)));

					for(int i = 0; i < number_of_pixels; i++) {
						int channel = 0, cluster = 0;
						getline(file, line);
						stringstream ist(line);
						ist >> cluster >> channel;
						ist >> x_pos[channel][cluster-1] >> y_pos[channel][cluster-1] >> bmp[channel][cluster-1];
					}
					vector<vector<double> > vector_pixel( 5, vector<double> (0));
					vector_pixel = cleaning(bmp);
					if(vector_pixel[0].size() > 3) {
						double *hillas;
						//cout << vector_pixel[0].size() << endl;
						hillas = get_hillas(vector_pixel, XToSourcecamatan2_cm, YToSourcecamatan2_cm, gam);
						fout << setw(6);
						fout1 << setw(0);
						fout << numb_clean << setw(6) << numb << setw(6) << i_scat << setw(6) << vector_pixel[0].size() << setw(13) << energy << setw(13) << XToSourcecamatan2_cm << setw(13) << YToSourcecamatan2_cm << setw(13) << x_ground << setw(13) << y_ground << setw(13) << hillas[12] << setw(13) << hillas[13] << endl;
						for(int i = 0; i < vector_pixel[0].size(); i++) {
							for(int kk = 0; kk < edge_pix[0].size(); kk++) {
								if(vector_pixel[0][i] == edge_pix[0][kk]) {
									if(vector_pixel[1][i] == edge_pix[1][kk]) {
										edge = 1;
									}
								}
							}
							fout << setw(0);
							for(int ii = 0; ii < 5; ii++) {
								fout << vector_pixel[ii][i] << setw(10);
							}
							fout << endl;
						}
						fout1 << numb_clean << setw(7) << numb << setw(7) << i_scat << setw(4) << vector_pixel[0].size() << setw(9) << energy << setw(9) << setw(13);
						for(int i = 0; i < 16; i++) {
							fout1 <<  hillas[i] << setw(13);
						}
						fout1 << edge << setw(13) << XToSourcecamatan2_cm << setw(13) << YToSourcecamatan2_cm << setw(13) << x_ground << setw(13) << y_ground << setw(13) << tet_center << setw(13) << fi_center << setw(13) << tet_source << setw(13) << fi_source << setw(13) << xmax << endl;
						numb_clean++;
					}
				}
			}
		}
		file.close();
		cout << edge_pix[0].size() << endl;
		edge_pix[0].clear();
		edge_pix[1].clear();
	}
	return 0;
}
