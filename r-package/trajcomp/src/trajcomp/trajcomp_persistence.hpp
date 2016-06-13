
#ifndef TRAJCOMP_PERSISTENCE_HPP
#define TRAJCOMP_PERSISTENCE_HPP

namespace trajcomp{

namespace tools{
	
		double toRadians(double val) {
			return val * (M_PI/180);
		}

		double toDegree(double val) {
			return val * (180/M_PI);
		}

		double calc_orientation(double lon1, double lon2, double lat1, double lat2)
		{
			lon1 = trajcomp::tools::toRadians(lon1);
			lon2 = trajcomp::tools::toRadians(lon2);
			lat1 = trajcomp::tools::toRadians(lat1);
			lat2 = trajcomp::tools::toRadians(lat2);
			//http://www.sunearthtools.com/tools/distance.php#txtDist_3
			double delta_lon = lon2-lon1;

			if(abs(delta_lon) > M_PI) {
				if(delta_lon > 0.0) {
					delta_lon = -(2.0 * M_PI - delta_lon);
				} else {
					delta_lon = (2.0 * M_PI -delta_lon);
				}
			}

			double delta_phi = log(tan((lat2/2) + (M_PI/4)) / tan((lat1/2) + (M_PI/4)));

			double orientation = atan2(delta_lon , delta_phi);
			//orientation = fmod(trajcomp::tools::toDegree(orientation)+360.0, 360.0);
			return trajcomp::tools::toDegree(orientation);
		}

};


	
	// Implementation provided by Jan Eicken (2015)

template<class TrajectoryType, class CurvType, class ThresholdType=double>
class persistence_segmentation_impl
{

public:
	/*
	*Diese Funktion berechnet den Krümmungsradius an jedem Punkt.
	* Der Abstand h sollte möglich klein gewählt werden
	*
	*/
	void calc_curvature(TrajectoryType &traj, CurvType &curvs)
	{
		curvs.clear();

		int first_orientation = (int)  fmod(trajcomp::tools::calc_orientation(traj[0][1], traj[1][1], traj[0][0], traj[1][0])+360, 360);
		curvs.push_back({first_orientation,0,0});

		for(size_t i = 1; i < traj.size()-1; i++) {

			int last_orientation = (int) fmod(trajcomp::tools::calc_orientation(traj[i-1][1], traj[i][1], traj[i-1][0], traj[i][0])+360, 360);
			int cur_orientation = (int) fmod(trajcomp::tools::calc_orientation(traj[i][1], traj[i+1][1], traj[i][0], traj[i+1][0])+360, 360);
			//cout << "last orientation: " << last_orientation << endl;
			//cout << "next orientation: " << cur_orientation << endl;
			int new_curv = 0;
			int opposite_orientation = (last_orientation+180) % 360;
			if(last_orientation > 180) {
				if(cur_orientation > opposite_orientation && cur_orientation < last_orientation) {
					new_curv = last_orientation - cur_orientation;
				} else {
					new_curv = -((cur_orientation - last_orientation + 360)%360);
				}
			} else {
				if(cur_orientation < opposite_orientation && cur_orientation > last_orientation) {
					new_curv = -(cur_orientation - last_orientation);
				} else {
					new_curv = (last_orientation - cur_orientation + 360) % 360;
				}
			}
			//cout << "New curv: " << new_curv << endl;
			curvs.push_back({new_curv, 0, 0});
		}

		curvs.push_back({0,0,0});
		//cout << "PRINT CURVS: " << curvs.size() << endl;
		for(size_t j=0; j < curvs.size(); j++) {
			//cout << j << ": " << curvs[j][0] << endl;
		}
	}
	/*
	* Berechnet lokale Minima und Maxima der gegebenen Trajektorie.
	* Die Eingabe muss aus Krümmungsradien bestehen, minima und maxima werden in übergebenen vectoren gescpeichert
	* In global_minima und global_maxima wird der Index der jeweiligen globalen extrema gespeichert
	*
	*/
	void calc_extrema(CurvType &curv, vector<size_t> &minima)
	{
		bool isClimbing = true;
		minima.clear();
		//Prüfe erstes Element der Trajektorie auf Minimum oder Maximum
		if(curv[0][0] < curv[1][0]) {
			minima.push_back(0);
			curv[0][2] = 1;
		} else {

			//curv[0][1] = 1;
			isClimbing = false;
		}
		//Prüfe innere Elemente auf Minima und Maxima
		for(size_t i = 1; i < curv.size(); i++) {
			//Last element
			if(i == curv.size()-1) {
				//cout << "last element: " << curv[i][0] << " climbing: " << isClimbing <<endl;
				if(isClimbing) {
					if(curv[i][0] > curv[i-1][0]) {
						curv[i][1] = 1;
					} else {
						minima.push_back(i);
						curv[i][2] = 1;
					}
				} else {
					if(curv[i][0] < curv[i-1][0]) {
						minima.push_back(i);
						curv[i][2] = 1;
					} else {
						curv[i][1] = 1;
					}
				}
			} else if(isClimbing) {
				if(curv[i][0] > curv[i+1][0]) {
					//is maxima
					//cout << "found maxima: " << curv[i][0] << endl;
					isClimbing = false;
					curv[i][1] = 1;
				}
			} else {
				if(curv[i][0] < curv[i+1][0]) {
					//is minima
					//cout << "found minima: " << curv[i][0] << endl;
					isClimbing = true;
					minima.push_back(i);
					curv[i][2] = 1;
				}
			}

		}

		//Sort minima array
		std::sort(std::begin(minima),
	                std::end(minima),
	                [&](int i1, int i2) { return curv[i1] < curv[i2]; } );
		/*for(size_t k = 0; k < minima.size(); k++) {
			cout << "min" << curv[minima[k]][0] << " ismin: " << curv[minima[k]][2] << " ismax: " << curv[minima[k]][1]<< endl;
		}
		for(size_t b = 0; b < curv.size(); b++) {
			cout << "point " << curv[b][0] << " ismin: " << curv[b][2] << " ismax: " << curv[b][1]<< endl;
		}
		*/
	}

	/* Component besteht aus:
	* - res[0] = Index des Minimums, dient gleichzeitig als Identifikator
	* - res[1] = Index des Anfangspunktes
	* - res[2] = Index des Endpunktes
	* - res[3] = Index des Maximums
	*/
	void calc_component(CurvType &curv, vector<size_t> &res){

		//Keep adding smallest neighbours until maximum is found
		//cout << "Suche nach punkten in Minimum component: " << curv[res[0]][0] << endl;
		while(true)
		{

			if((res[1] >= 1) && (res[2] <= curv.size()-2)) {


				if(curv[res[1]-1][0] < curv[res[2]+1][0]) {
					res[1]--;
					//Maximum found left
					if(curv[res[1]][1] == 1) {
						res[3] = res[1];
						//cout <<  " Maximum left1: " << curv[res[1]][0] << endl;

						return;
					}
				} else {
					res[2]++;
					//Maximum found right
					if(curv[res[2]][1] == 1) {
						res[3] = res[2];
						//cout <<  " Maximum right1: " << curv[res[2]][0] << endl;
						return;
					}
				}
			} else if (res[1] >= 1) {
				res[1]--;
				//Maximum found left
				if((curv[res[1]][1] == 1) && (res[1] > 0)) {
					res[3] = res[1];
					//cout <<  " Maximum left: " << curv[res[1]][0] << endl;

					return;
				}


			} else if (res[2] <= curv.size()-2) {
				res[2]++;
				//Maximum found right
				if((curv[res[2]][1] == 1) && (res[2]) < curv.size()-1) {

					res[3] = res[2];
					//cout <<  " Maximum right: " << curv[res[2]][0] << " of minimum: " << curv[res[0]][0] << endl;

					return;
				}

			} else {
				return;
			}

		}
	}

	/*
	* Berchnet die für den online algorithmus gebrauchten bars
	* Problem zur Zeit noch: Jedes Mal wenn ein Maximum in zwei Komponenten gefunden wurde, setze ich das Maximum-Flag im
	* curvature vector auf 0 füge das Minimum wieder in den minima vector ein. d.h. aber auch auch, dass die Komponente wieder von
	* Anfang an aufgebaut wird.
	*/
	void calc_bars(CurvType &curv, vector<size_t> &minima, vector<vector<size_t> > &bars)
	{
	
//		cout << "Entering calc_bars" << endl;
		vector<vector<size_t>> components;
		components.clear();
		bars.clear();
		//Gehe durch alle minima von geringstem bis höchstem
		for(size_t i = 0; i < minima.size(); i++){
	//		cout << "Sweep line from minimum: " << curv[minima[i]][0] << endl << i << "/" << minima.size() << endl << endl;
			//Erstelle neue Komponente mit dem gefundenen Mimimumswert als Standardwert
			//@WRONG: martin added hack to avoid ever growing minima chased by i == minima.size() -1
			if (i == minima.size() -1) 
			   break;
			vector<size_t> new_comp(4);
			new_comp.clear();
			new_comp.push_back(minima[i]);
			new_comp.push_back(minima[i]);
			new_comp.push_back(minima[i]);
			new_comp.push_back(minima[i]);
			bool foundOld = false;
			bool foundMaxima = false;

			//Suche nach Minimum in bestehenden Komponenten. Setze derzeitige Komponente auf bestehende falls gefunden.
			for(size_t n = 0; n < components.size(); n++) {
				//cout << "Component: " << curv[components[n][0]][0] << " minimum: " << curv[minima[i]][0] << endl;
				if(components[n][0] == minima[i]) {
					//cout << "Found component with minimum" <<  curv[minima[i]][0] << endl;
					new_comp.assign(components[n].begin(), components[n].end());
					foundOld = true;
					components.erase(components.begin() + n);
				}
			}
		//	cout << "Marker before calc_component" << endl;
			calc_component(curv, new_comp);
		//	cout << "Marker after calc_component" << endl;
			//cout << endl;
			//denotates if a component with same maximum was found
			//cout << "searching for components with same maximum..." << endl;
			for(size_t j = 0; j < components.size(); j++) {
				//cout << "New comp left: " << new_comp[1] << " right: " << new_comp[2] << endl;
				//cout << "Old comp left: " << components[j][1] << " right: " << components[j][2] << endl;
				if(new_comp[3] == components[j][3]) {
					//cout << "found same maximum: "<< curv[new_comp[3]][0] << " for " << curv[new_comp[0]][0] << " and " << curv[components[j][0]][0] << endl;


					if(curv[new_comp[0]][0] < curv[components[j][0]][0]) {
						vector<size_t> new_bar = {components[j][0], components[j][3]};
						bars.push_back(new_bar);

						//Set left part of component
						if(components[j][0] > new_comp[0]) {
							components[j][1] = new_comp[1];
						//Set right part of component
						} else {
							components[j][2] = new_comp[2];
						}
						components[j][0] = new_comp[0];
						//minima.push_back(new_comp[0]);


					} else {
						vector<size_t> new_bar = {new_comp[0], new_comp[3]};
						bars.push_back(new_bar);

						//Set left part of component
						if(components[j][0] > new_comp[0]) {
							components[j][1] = new_comp[1];
						//Set right part of component
						} else {
							components[j][2] = new_comp[2];
						}
						//minima.push_back(components[j][0]);
					}
					minima.insert(minima.begin()+i+1, components[j][0]);

					//cout << "New Component is from " << curv[components[j][1]][0] << " to " << curv[components[j][2]][0] << endl;
					foundMaxima = true;
					//curv[components[j][0]][3] == 0;
					//the found components were merged, minimum has to be inserted into minima vector again

				}
			}
			if(!foundMaxima) {
				//cout << "No component found, add new one: min " << curv[new_comp[0]][0] << " max " << curv[new_comp[1]][0] << endl;
				components.push_back(new_comp);
			}
			if((i == minima.size()-1) && (components.size() > 1)) {
				//cout << "last minimum reached, components left" << endl;
				//Set possible maximum Minimum to minimal Minimum
				size_t maxMin = minima[0];
				for(size_t n = 0; n < components.size(); n++) {
					//cout << "component with minimum: " << curv[components[n][0]][0] << endl;
					if(curv[maxMin][0] < curv[components[n][0]][0]) {
						maxMin = components[n][0];
					}
				}
				//cout << "insert min: " << curv[maxMin][0] << endl;
				minima.insert(minima.begin()+i+1, maxMin);
			}
			//cout << endl;

		}
		//Add last component as bar
		//bars.push_back({components[components.size()-1][0], components[components.size()-1][3]});
		/*
		for(size_t k = 0; k < components.size(); k++) {

			cout << "component #" << k << " start: " << curv[components[k][0]][0] << " end: " << curv[components[k][3]][0] << endl;
		}

		for (size_t f = 0; f < bars.size(); f++) {
			cout << "bar: "<< f << " left: " << curv[bars[f][0]][0] << " right: " << curv[bars[f][1]][0] << endl;
		}
		*/

	}
	
	#define TRAJCOMP_PERSISTENCE_DEBUG

	TrajectoryType operator()(TrajectoryType &traj, ThresholdType beta)
	{
		//Schritt 1: Ableitung einer Kurvenfunktion


		//TODO: delete, only for testing
		CurvType curv(3);
		curv.clear();
		calc_curvature(traj, curv);
		cout << "curv length: " << curv.size() << endl;
		//initialiseTestVector(curv);

		//Schritt 2: Berechne maxima und minima
		vector<size_t> minima;
		calc_extrema(curv, minima);
		//cout << "done with extrema" << endl;
		for(size_t b = 0; b < curv.size(); b++) {
			//cout << "Point: " << curv[b][0] << " min: " << curv[b][2] << " max: " << curv[b][1] << endl;
		}
		//cout << "Extrema complete" << endl;
		
		vector<vector<size_t> > bars;
		calc_bars(curv, minima, bars);

		//cout << "Bars are complete" << endl;

		bool startPoint = false;
		bool endPoint = false;
		vector<double> res(1);
		res.clear();
		//beta-persistent simplification
		//cout << "bars size: " << bars.size() << endl;
		for(size_t i = 0; i < bars.size(); i++) {
			if(abs(curv[bars[i][0]][0] - curv[bars[i][1]][0]) < beta) {
				//cout << "Beta bigger: " << curv[bars[i][0]][0]+ curv[bars[i][1]][0] << endl;
				
				bars.erase(bars.begin() + i);
				i--;
			} else {
				//cout << "Beta smaller" << endl;
				//Check if end points already are in simplification
				if(bars[i][0] == 0 || bars[i][1] == 0) {
					//cout << "start point found" << endl;
					startPoint = true;
				}
				if(bars[i][1] == curv.size()-1 || bars[i][0] == curv.size()-1) {
					//cout << "end point found" << endl;
					endPoint == true;
				}
				//cout << "Push to result set: " << curv[bars[i][0]][0] << " + " << curv[bars[i][1]][0] << endl;
				res.push_back(bars[i][0]);
				res.push_back(bars[i][1]);
			}
		}
		sort(res.begin(), res.end());
		//cout << endl;
		//cout << "size: " << res.size() << endl;
		//cout << "Bars before start point check" <<endl;
		for(size_t j = 0; j < res.size(); j++) {
			//cout << curv[res[j]][0] << endl;
		}

		//Special case for start and end point
		if(!startPoint) {
			//Check points after start point and insert
			size_t toInsert = 1;
			bool insert = false;
			if(curv[res[0]][1] == 1){
				//first inserted extremum is maximum
				for(size_t k = 1; k < res[0]; k++) {
					if((curv[k][0] < curv[toInsert][0]) && (curv[k][2] == 1)) {
						//Inserted point between start and first inserted has to be minimum
						toInsert = k;
						insert = true;
					}
				}
			} else {
				//first insterted extremum is minimum
				for(size_t j = 1; j < res[0]; j++) {
					if((curv[j][0] > curv[toInsert][0]) && (curv[j][1] == 1)) {
						//Inserted point between start and first inserted has to be maximum
						toInsert = j;
						insert = true;
					}
				}
			}
			if(toInsert > -1) {
				res.insert(res.begin(), toInsert);
			}
			res.insert(res.begin(), 0);

		}
		if(!endPoint) {
			//cout << "End point not in simplification" << endl;
			//Check points before an end point and insert
			size_t toInsert = curv.size()-2;
			bool insert = false;
			if(curv[res[res.size()-1]][1] == 1){
				//last inserted extremum is maximum
				for(size_t k = curv.size()-1; k > res[res.size()-1]; k--) {
					//cout << "check last points" << endl;
					if((curv[k][0] < curv[toInsert][0]) && (curv[k][1] == 1)) {
						toInsert = k;
						insert = true;
					}
				}
			} else {
				//last insterted extremum is minimum
				for(size_t j = curv.size()-1; j > res[res.size()-1]; j--) {
					if((curv[j][0] > curv[toInsert][0]) && (curv[j][2] == 1)) {
						toInsert = j;
						insert = true;
					}
				}
			}
			if(insert){
				//cout << "found point before last point: " << toInsert << endl;
				res.insert(res.end(), toInsert);
			}
			res.insert(res.end(), curv.size()-1);
		}
		TrajectoryType resTraj;
		resTraj.clear();
		for(size_t n = 0; n < res.size(); n++) {

			//cout << traj[res[n]][0] << " " << traj[res[n]][1] << endl;
			resTraj.push_back({traj[res[n]][0], traj[res[n]][1]});
		}

		return resTraj;
	}
};


template<class TrajectoryType, class ThresholdType=double>
TrajectoryType persistence_segmentation(TrajectoryType &traj, ThresholdType beta)
{
	 TrajectoryType res(2);
	persistence_segmentation_impl<TrajectoryType, vector< vector<int> >, ThresholdType> ps;
	res = ps(traj, beta);
	return res;
}


};//namespace



#endif
