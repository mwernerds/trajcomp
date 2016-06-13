#ifndef PUBLIC_API_H_INC
#define PUBLIC_API_H_INC
#include<vector>
#include<string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace pt = boost::property_tree;


/* Pre-compiled Settings */
struct Settings
{
	pt::ptree tree;
  std::vector<std::string> distances;
  std::vector<int> lengths;
  std::vector<std::vector<std::string>> features;
};


// Some forward decls
std::string feature_dispatchbyname(std::vector<std::vector<double>> &q, std::string name);
double string_distance_dispatchbyname(std::string &t, std::string &q, std::string name, int k);
double distance_dispatchbyname(std::vector<std::vector<double>> &t, std::vector<std::vector<double>> &q, std::string name);
std::string zip_features(std::vector<std::vector<double>> &q, std::vector<std::string> featureNames);

Settings &getSettings(int i);


#endif
