

/*template <class TrajectoryType, class DistanceType=double, class PtsType=int>
class DBSCAN_segmentation_impl
{
	public:
		std::string exec(const char* cmd)
		{
			FILE* pipe = popen(cmd, "r");
		    if (!pipe) return "ERROR";
		    char buffer[128];
		    std::string result = "";
		    while(!feof(pipe)) {
		    	if(fgets(buffer, 128, pipe) != NULL)
		    		result += buffer;
		    }
		    pclose(pipe);
		    return result;
		}
		double get_mean_in_cluster(vector<vector<double>> &coords)
		{
			double mean = 0.0;
			for(size_t i = 0; i < coords.size(); i++) {
				mean += coords[i][1];
			}
			mean = mean/coords.size();
			size_t meanID = coords.size()-1;
			for(size_t j = 0; j < coords.size()-1; j++) {
				//std::cout << meanID << " ID; Current diff to mean: " << fabs(coords[meanID][1] - mean);
				//std::cout << " " << coords[j][0] << " ID; Next diff to mean: " << fabs(coords[j][1] - mean) << endl;
				if(fabs(coords[j][1] - mean) < fabs(coords[meanID][1] - mean)) {
				//	std::cout << "FOUND SMALLER " << coords[j][0] << endl << endl;
					meanID =  j;
				}
			}
			//std::cout << "Found meanest: " << coords[meanID][0] << endl;
			return coords[meanID][0];
		}

		bool sortByID(const vector<double> &a, const vector<double> &b){return a[0] < b[0];}

		TrajectoryType operator()(TrajectoryType &traj, DistanceType epsilon,
			PtsType minPts, const std::string &dataName, const std::string &pathToData,
			const std::string path_to_elki, bool isDBSCAN)
		{
			//Schritt 1: sample data (Evtl auch davor in der main aufrufen)
			std::cout << "start dbscan..." << std::endl;
			//Schritt 2: Rufe DBSCAN auf übergebenen Datenpfad auf
			std::stringstream cmd;
			cout <<  "data name: " << dataName << endl;
			cout << "path to data: " << pathToData << endl;
			if(isDBSCAN) {
				cmd << "java -jar " << path_to_elki << " KDDCLIApplication -dbc.in " << pathToData <<
				" -algorithm clustering.DBSCAN -dbscan.epsilon " << epsilon << " -dbscan.minpts " <<
				minPts << " -resulthandler ResultWriter -out clustering/" << dataName;
			} else {
				cmd << "java -jar " << path_to_elki << " KDDCLIApplication -dbc.in " << pathToData <<
				" -algorithm clustering.optics.OPTICS " << " -optics.minpts " <<
				minPts << " -resulthandler ResultWriter -out clustering/" << dataName;
				 TrajectoryType myres(2);
				return myres;
			}

			const std::string& tmp = cmd.str();
			const char* p = tmp.c_str();
			string result = exec(p);

			cout << "done clustering" << endl;

			//TODO Überprüfen, ob bereits in trajektorie vorhanden
			//Erster Punkt von Trajektorie muss im Ergebnis sein
			vector<size_t> resIDs;
			TrajectoryType res(2);
			res.clear();

			//{ID, lat+lon}
			vector<vector<double> > coords(2);

			try {
				//Gehe gefundene Cluster in Ordner durch
				int i = 0;
				while(true)
				{
					std::stringstream path_name;
					path_name << "clustering/" << dataName << "/cluster_" << i << ".txt";
					std::ifstream file(path_name.str().c_str());
					if(!file) {
						cout << "No more clusters found... " << endl;
						break;
					} else {
						cout << "found cluster: " << path_name.str() << endl;
					}
				//cout << "Found file " << endl;
					coords.clear();
					string line;
					while(getline(file, line)) {
						size_t found = line.find("#");

						if(found != string::npos){
							continue;
						}
						std::istringstream iss(line);
						string IDString;
						double lat, lon;
						if(iss >> IDString >> lat >> lon) {
						} else {
							throw(std::runtime_error("Error: Can't convert values in cluster."));
						}
						IDString = IDString.substr(IDString.find("=")+1, string::npos);
						iss.str("");
	 					iss.clear();
						iss.str(IDString);
						double ID;
						iss >> ID;
						//std::cout << "CLUSTERID: " << ID << std::endl;
						//std::cout << lat << " : " << lon << std::endl;
						coords.push_back({ID,lat+lon});
					}
					//Schritt 3: Suche repräsentierbare Punkte in Cluster
					size_t meanestID = get_mean_in_cluster(coords);
					cout << "push id back: " << meanestID << endl;
					resIDs.push_back(meanestID);
			    	file.close();
					i++;
				}

				//Letzter Punkt in Trajektorie muss im Ergebnis sein
				//if(std::find(std::begin(resIDs), std::end(resIDs), traj.size()-1) != std::end(resIDs)) {
				//	resIDs.push_back(traj.size()-1);
				//}
				if(resIDs.size() > 0) {
					std::sort(resIDs.begin(), resIDs.end());
					for(int i = 0; i < resIDs.size(); i++) {
						res.push_back(traj[resIDs[i]]);
					}
				} else {
					cout << "no clusters found" << endl;
				}
			}catch (const std::bad_alloc&) {
			  return res;
			}

			cout << "return clustering results" << endl;
			return res;

		}

};

template <class TrajectoryType>
TrajectoryType DBSCAN_segmentation(TrajectoryType &traj, double epsilon,
	int minPts, string dataName, string pathToData, string path_to_elki, bool isDBSCAN)
{
	DBSCAN_segmentation_impl<TrajectoryType, double, int> dbs;
	//boost::filesystem::remove_all("clustering/");
	system("exec rm -r clustering/*");
	TrajectoryType ret = dbs(traj,epsilon, minPts, dataName, pathToData, path_to_elki, isDBSCAN);

	//TODO: What to do with noise points?

	return ret;
}
* 
* */

//--------------------START THRESHOLD SAMPLING---------------------------


template <class TrajectoryType, class TimesType, class ThresholdType=double>
class threshold_sampling_impl
{
	public:

		TrajectoryType operator()(TrajectoryType &traj, TimesType &times,
			ThresholdType velocity_thresh, ThresholdType orientation_thresh)
		{
			std::cout << std::setprecision(10);
			TrajectoryType result(2);
			result.clear();
			TimesType resultTimes(1);
			resultTimes.clear();
			result.push_back({traj[0][0], traj[0][1]});

			resultTimes.push_back(times[0]);
			result.push_back({traj[1][0], traj[1][1]});
			resultTimes.push_back(times[1]);

			//calculate safe area
			for(size_t i = 2; i<traj.size(); i++) {
				size_t resultSize = result.size();

				//START Calculate upper and lower angle and radius for trajectory and sample based last two points.

				//For velocity vector of trajectory, take last 2 elements in original trajectory
				vector<double> velocity_vector_traj = trajcomp::tools::calc_velocity_vector(traj[i-2][1], traj[i-1][1],
													traj[i-2][0], traj[i-1][0], times[i-2], times[i-1]);
				//For velocity vector of sample, take last 2 elements in result trajectory
				vector<double> velocity_vector_sample = trajcomp::tools::calc_velocity_vector(result[resultSize-2][1], result[resultSize-1][1], result[resultSize-2][0],
												result[resultSize-1][0], resultTimes[resultSize-2], resultTimes[resultSize-1]);

				//Calculate the radius for trajectory and sample based using the calculated velocity vectors
				double upper_radius_traj = ((times[i] - times[i-1]) / 1000.0) *
												(velocity_vector_traj[0]*(1 + velocity_thresh));
				double lower_radius_traj = ((times[i] - times[i-1]) / 1000.0) *
												(velocity_vector_traj[0]*(1 - velocity_thresh));

				double upper_radius_sample = ((times[i] - times[i-1]) / 1000.0) *
												(velocity_vector_sample[0]*(1 + velocity_thresh));
				double lower_radius_sample = ((times[i] - times[i-1]) / 1000.0) *
												(velocity_vector_sample[0]*(1 - velocity_thresh));

				cout << "Upper radius traj: " << upper_radius_traj << " lower: " << lower_radius_traj << endl;
				cout << "Uppe radius sample: " << upper_radius_sample << " lower: " << lower_radius_sample << endl;
				//Calculate orientation angle for both from given slopes
				double upper_slope_traj = velocity_vector_traj[1] + orientation_thresh;
				double lower_slope_traj = velocity_vector_traj[1] - orientation_thresh;


				double upper_slope_sample = velocity_vector_sample[1] + orientation_thresh;
				double lower_slope_sample = velocity_vector_sample[1] - orientation_thresh;
				cout << "Upper slope traj: " << upper_slope_traj << " lower: " << lower_slope_traj << endl;
				cout << "Uppe slope sample: " << upper_slope_sample << " lower: " << lower_slope_sample << endl;
				//START Calculate angle and radius for point to insert
				vector<double> velocity_vector_insert = trajcomp::tools::calc_velocity_vector(traj[i-1][1], traj[i][1],
													traj[i-1][0], traj[i][0], times[i-1], times[i]);

				double radius_insert = ((times[i] - times[i-1]) / 1000.0) *
												velocity_vector_insert[0];
				cout << "Radius insert: " << radius_insert <<endl;
				double angle_insert = velocity_vector_insert[1];
				cout << "angle insert: " << angle_insert << endl;

				//Check if new points angle and radius is inside of the calculated
				bool isInRadius = false;
				bool isInAngle = false;
				if((radius_insert >= lower_radius_traj) && (radius_insert >= lower_radius_sample)
						&& (radius_insert <= upper_radius_traj) && radius_insert <= upper_radius_sample) {
					cout << "Trajectory: " << i << " is in Radius: " << endl;
					isInRadius = true;
				}

				if((angle_insert <= upper_slope_traj) && (angle_insert <= upper_slope_sample)
						&& (angle_insert >= lower_slope_sample) && (angle_insert >= lower_slope_traj)) {
					isInAngle = true;
				}
				if(!isInAngle || !isInRadius) {
					result.push_back({traj[i][0], traj[i][1]});
					resultTimes.push_back(times[i]);
				}
			}

			cout << "Original size: " << traj.size() << " threshold sampling result size: " << result.size() << endl;

			return result;
		}
};


template<class TrajectoryType, class TimesType>
TrajectoryType threshold_sampling(TrajectoryType &traj, TimesType &times, double velocity_thresh, double orientation_thresh)
{
	threshold_sampling_impl<TrajectoryType, TimesType, double> ts;

	TrajectoryType ret = ts(traj, times, velocity_thresh, orientation_thresh);
	return ret;
}
