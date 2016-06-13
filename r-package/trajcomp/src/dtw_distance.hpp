#ifndef DTW_DISTANCE_HPP
#define DTW_DISTANCE_HPP


double get_distance(std::vector<double> &a, std::vector<double> &b)
{
  double xdiff = a[0] - b[0];
  double ydiff = a[1] - b[1];
  return xdiff*xdiff + ydiff*ydiff;
}

namespace dtw_distance {

  // dtw implementation corresponding to "Fast Similarity Search of Multidimensional Time Series", X.Gong et al., 2015
  double dtw(std::vector<std::vector<double>> &t, std::vector<std::vector<double>> &q, double wc = 0.1){

    int n = t.size();
    // warping constraint
    int c = wc*n;

    double dm[n*n]; // distance "matrix"
    std::fill(&dm[0], &dm[0]+(n*n), std::numeric_limits<double>::infinity());


    double tmpDistance, accumulatedDistance = 0.0;


    dm[0] = get_distance(t[0], q[0]);

    // fill in the distance matrix using dynamic programming
    for (int i = 1; i < n; i++){
      for (int j = std::max(1,i-c); j <= std::min(n-1,i+c); j++){
        tmpDistance = std::min(std::min(dm[(i-1)*n + j], dm[(i-1)*n + j-1]), dm[i*n + j-1]);
        dm[i*n + j] = tmpDistance + get_distance(t[j], q[j]);
      }
    }
    accumulatedDistance = dm[(n - 1)*n + n - 1];

    return sqrt(accumulatedDistance);
  }
}
#endif
