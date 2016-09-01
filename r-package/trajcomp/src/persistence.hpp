#include <Rcpp.h>
#include <cassert>
using namespace Rcpp;


namespace persistence{

//basic vector math
typedef std::vector<double> vec2;
typedef std::vector<vec2> Trajectory;

double length_sqr(const vec2 & v){
  return v[0] * v[0] + v[1] * v[1];
}
double length(const vec2 & v){
  return sqrt(length_sqr(v));
}

double dot(const vec2 & v, const vec2 & w){
  return v[0] * w[0] + v[1] * w[1];
}

//positive if v is left of w negative else
double left_right(const vec2 & v, const vec2 & w){
  return v[0] * w[1] - v[1] * w[0];
}

//1.0 if v is left of w -1 else
double lr_sgn(const vec2 & v, const vec2 & w){
  return (left_right(v, w) < 0) ? -1.0: 1.0;
}

double to_deg(double rad){
  return (rad * 180.0) / M_PI;
}

vec2 mult(const vec2 & v, double x){
  return {v[0] * x, v[1] * x};
}

void normalize(vec2 & v){
  auto len = length(v);
  if(len != 0){
    v[0] *= 1/len;
    v[1] *= 1/len;
  }
}

vec2 diff(const vec2 & v, const vec2 & w){
  return {v[0] - w[0], v[1] - w[1]};
}
vec2 sum(const vec2 & v, const vec2 & w){
  return {v[0] + w[0], v[1] + w[1]};
}

double distance(const vec2 & v, const vec2 & w){
  return length(diff(v,w));
}

double minimum_distance(const vec2 &v, const vec2 &w, const vec2 &p) {
  // Return minimum distance between line segment vw and point p
  const double l2 = length_sqr(diff(v, w));  // i.e. |w-v|^2 -  avoid a sqrt
  if (l2 == 0.0) return distance(p, v);   // v == w case
  // Consider the line extending the segment, parameterized as v + t (w - v).
  // We find projection of point p onto the line. 
  // It falls where t = [(p-v) . (w-v)] / |w-v|^2
  // We clamp t from [0,1] to handle points outside the segment vw.
  const double t = std::max(0.0, std::min(1.0, dot(diff(p , v), diff(w , v)) / l2));
  const vec2 projection = sum(v ,mult( diff(w , v), t));  // Projection falls on the segment
  return distance(p, projection);
}

struct CurveElement{
  double val;
  int index;
  vec2 vertex;
};

typedef std::vector<CurveElement> Curve;

void print_curve(Curve &c){
  for(int i = 0; i < c.size(); i++){
    if(c[i].vertex.size() == 2){
      std::cout << "v["<<i<<"].("<<c[i].index<<"): "<< c[i].vertex[0] <<", "<<c[i].vertex[1]
                << " = " << c[i].val
                << std::endl;
    }
    else{
      std::cout << "v["<<i<<"]: EMPTY!"<< c[i].vertex.size()<<" - " << c[i].index << std::endl;
    }
  }
  
}
Curve traj_to_curve(const std::vector<vec2> &t){
  auto res = Curve(t.size());
  
  for(auto i = 0; i < t.size(); ++i){
    //edge cases
    if(i == 0 || i == t.size()-1){
      res[i] = {0, i, t[i]};
    }else{
      auto v = t[i];
      auto last = t[i-1];
      auto v1 = diff(last, v);
      auto next = t[i+1];
      auto v2 = diff(v, next);
      
      auto len1 = length(v1);
      v1[0] *= 1/len1;
      v1[1] *= 1/len1;
      
      auto len2 = length(v2);
      v2[0] *= 1/len2;
      v2[1] *= 1/len2;
      
      auto value = dot(v1, v2);
      auto deg = to_deg( acos(value) * lr_sgn(v1, v2) );
      if(len1 == 0 || len2 == 0){
        deg = 0;
      }
      res[i] = {deg, i, t[i]};
      
      //std::cout << i<<": "<<v1[0]<<", "<<v1[1] <<" | "<<v2[0]<<", "<<v2[1] <<" | "<<deg<<std::endl;
    }
  }
  //print_curve(res);
  return res;
};

struct Extrema{
  std::vector<int> min, max;
};

Extrema extrema(const Curve &curve){
  
  const int fl_un = 0, fl_min = 1, fl_max = 2;
  Extrema res;
  
  if(curve.size() < 3){
    std::cout <<"CURVE IS EMPTY"<<std::endl;
    return res;
  }
  
  auto d_min = curve[0], d_max = curve[0];
  auto search = fl_un;
  auto i_un = 0;
  
  for(auto i = 1; i < curve.size(); ++i){
    auto x = curve[i];
    if(search == fl_un){
      if(x.val > d_min.val){
        search = fl_max;
        d_max = x;
        res.min.push_back(i_un);
      }
      if(x.val < d_min.val){
        search = fl_min;
        d_min = x;
        res.max.push_back(i_un);
      }
    } else if(search == fl_min){
      if(x.val > d_min.val){
        search = fl_max;
        d_max = x;
        res.min.push_back(i-1);
      }else{
        d_min = x;
      }
    } else if(search == fl_max){
      if(x.val < d_max.val){
        search = fl_min;
        d_min = x;
        res.max.push_back(i-1);
      }else{
        d_max = x;
      }
    }
  }
  
  //if no min or max found
  if(res.min.size() == 0 || res.max.size() == 0){
    res.min.push_back(0);
    res.max.push_back(curve.size()-1);
    return res;
  }
  
  //add end, its either min or max
  if(res.min.back() > res.max.back()){
    res.max.push_back(curve.size()-1);
  }else{
    res.min.push_back(curve.size()-1);
  }
  return res;
}

struct Component{
  int left, right;
  int min, max;
  bool finished;
};

struct Bar{
  int min, max;
};

struct Result{
  std::vector<Bar> bars;
  std::vector<Component> comps;
  std::vector<int> used;
  Curve pruned;
};

Curve betaPruning(const std::vector<Bar> &bars, const Curve &curve, double beta = 0.01){
  //insertion sort structure: std::set<>
  std::set<CurveElement, std::function<bool (const CurveElement&, const CurveElement&)>> pruned([](const CurveElement& a, const CurveElement& b){
    return a.index < b.index;
  });
  
  auto end_p = curve.size()-1;
  pruned.insert(curve[0]);
  pruned.insert(curve[end_p]);
  
  for(auto b : bars){
    //prune by beta
    if(curve[b.max].val - curve[b.min].val > beta){
      //avoid duplicate start and end
      if(b.min != 0 && b.min != end_p ){
        pruned.insert(curve[b.min]);
      }
      if(b.max != 0 && b.max != end_p ){
        pruned.insert(curve[b.max]);
      }
    }
  }
  //copy out sorted curve
  Curve pruned_curve;
  std::copy(pruned.begin(), pruned.end(), std::back_inserter(pruned_curve));
  return pruned_curve;
}

Result PersistenceAlg(const Curve &curve, double beta = 0.01, int debug_iter = -1){
  
  auto ex = persistence::extrema(curve);
  auto min = ex.min;
  auto max = ex.max;
  
  if(curve.size() < 3 || min.size() == 0 || max.size() == 0){
    std::cout <<"EMPTY CURVE"<<std::endl;
    Result x;
    x.pruned = curve;
    return x;
  }
  auto used = std::vector<int>(curve.size());
  std::fill(used.begin(), used.end(), -1);
  
  auto max_set = std::set<int>(max.begin(), max.end());
  
  auto comps = std::vector<Component>();
  for(auto m:min){
    comps.push_back({m,m,m,-1, false});
  }
  auto active_comps = std::vector<int>();
  for(int i = 0; i < comps.size(); ++i){
    active_comps.push_back(i);
  }
  
  
  //setup loop
  bool finished = false;
  int iter = -1;
  if(debug_iter < 0){
    debug_iter = curve.size()*(int)comps.size();
  }
  
  int max_iter = std::min((int)curve.size()*(int)comps.size(), debug_iter);
  
  while (! finished && ++iter < max_iter && active_comps.size() > 0){
    
    //grow bars
    int i = iter%((int) active_comps.size());
    
    //std::cout<<"iter: "<<iter << " %% "<< i<<"/"<<active_comps.size() <<std::endl;
    
    auto current_i = active_comps[i];
    auto &c1 = comps[current_i];
    
    if(c1.left == 0 && c1.right == (int)curve.size()-1){
      
      finished = true;
      c1.finished = true;
      //REMOVE
      active_comps.erase(active_comps.begin()+i);
      
      if(max_set.count(c1.left) > 0){
        c1.max = c1.left;
      }
      else if(max_set.count(c1.right) > 0){
        c1.max = c1.right;
      }
      
      //std::cout<<"Finished!!!!!!!!!!!!!!!! at: "<< iter<<std::endl;
    }
    
    if(!c1.finished){
      auto l = c1.left - 1;
      auto r = c1.right + 1;
      
      //decide grow direction
      int x = 0;
      if((curve[l].val < curve[r].val || r == curve.size())  && l != -1  ){
        c1.left = l;
        x = l;
      }else{
        c1.right = r;
        x = r;
      }
      
      if(x== -1){
        x=0;
      }
      if(x == curve.size()){
        x = curve.size()-1;
      }
      
      //if is _max(x)
      if(max_set.count(x) > 0 ){
        
        c1.finished = true;
        c1.max = x;
        //REMOVE
        active_comps.erase(active_comps.begin()+i);
        
        //if max is used
        if(used[x] != -1){
          
          //merge
          auto & c2 = comps[used[x]];
          auto m1 = c1.min;
          auto m2 = c2.min;
          if (c2.min > c1.min){
            m1 = c2.min;
            m2 = c1.min;
          }
          
          c2.left = std::min(c1.left, c2.left);
          c2.right = std::max(c1.right, c2.right);
          c2.min = m2;
          c2.max = -1;
          c2.finished = false;
          //ADD
          active_comps.push_back(used[x]);
          
          c1.left = std::min(x, m1);
          c1.right = std::max(x, m1);
          c1.min = m1;
          //REDUNDANT see above
          c1.max = x;
          c1.finished = true;
          ////REMOVE
          //active_comps.erase(active_comps.begin()+i);
        }
        used[x] = current_i;
      }
    }
    //}
  }
  
  //make bars
  std::vector<Bar> bars;
  
  for(auto c : comps){
    if(c.max >=0){
      std::cout << " b: " << c.min <<", "<< c.max <<std::endl;
      bars.push_back({c.min,c.max});
    }
  }
  
  return {bars, comps, used, persistence::betaPruning(bars, curve, beta)};
};

Curve prune_curve_dist(const Curve &curve, double scale){
  if(curve.size() < 2){
    return curve;
  }
  Curve res;
  res.reserve(curve.size());
  res.push_back(curve[0]);
  
  for(auto i = 1; i < curve.size()-1;  ++i){
    auto it = curve[i], it2 = curve[i+1];
    //check distance
    if(distance(it.vertex, it2.vertex) < scale){
      //pruning
      if(it.val < it2.val){
        //avoid double insertion
        if(it.index != res.back().index)
          res.push_back(it);
      }else {
        res.push_back(it2);
      }
    }else{
      //no pruning
      res.push_back(it);
    }
  }
  res.push_back(curve.back());
  return res;
}

Curve prune_curve_dist_to_segment(const Curve &curve, double epsilon){
  if(curve.size() < 2){
    return curve;
  }
  Curve res;
  res.reserve(curve.size());
  res.push_back(curve[0]);
  vec2 last_pos = curve[0].vertex;
  for(auto i = 1; i < curve.size()-1;  ++i){
    auto it = curve[i], it2 = curve[i+1];
    //check distance
    if(minimum_distance(last_pos,it2.vertex, it.vertex) < epsilon){
      last_pos = it2.vertex;
      res.push_back(it2);
      ++i;
    }else{
      //no pruning
      last_pos = it.vertex;
      res.push_back(it);
    }
  }
  res.push_back(curve.back());
  return res;
}


void recalc_angles_inplace(Curve &curve){
  for(auto i = 0; i < curve.size(); ++i){
    //edge cases
    if(i == 0 || i == curve.size()-1){
      curve[i].val = 0;
    }else{
      auto v = curve[i].vertex;
      auto last = curve[i-1].vertex;
      auto v1 = diff(last, v);
      auto next = curve[i+1].vertex;
      auto v2 = diff(v, next);
      
      auto len1 = length(v1);
      v1[0] *= 1/len1;
      v1[1] *= 1/len1;
      
      auto len2 = length(v2);
      v2[0] *= 1/len2;
      v2[1] *= 1/len2;
      
      
      auto value = dot(v1, v2);
      auto deg = to_deg( acos(value)*lr_sgn(v1, v2) );
      if(len1 == 0 || len2 == 0){
        deg = 0;
      }
      curve[i].val = deg;
    }
  }
}

Curve persistenceMultiRes(const Curve &curve,  double beta, int levels){
  auto _curve = curve;
  std::cout<<"START"<<std::endl;
  print_curve(_curve);
  
  for(int i = 0; i < levels; ++i){
    if(i != 0){
      recalc_angles_inplace(_curve);
    }
    auto p_result = persistence::PersistenceAlg(_curve, beta, -1);
    std::cout<<"IT: "<<i<<" pruned: "<<std::endl;
    print_curve(p_result.pruned);
    _curve = persistence::prune_curve_dist(p_result.pruned, pow(2,i));
    std::cout<<"IT: "<<i<<std::endl;
    print_curve(_curve);
    
  }
  std::cout<<"END"<<std::endl;
  print_curve(_curve);
  return _curve;
}

Curve persistenceDist(const Curve &curve,  double beta, double epsilon, int iterations){
  auto _curve = curve;
  //print_curve(_curve);
  
  auto p_result = persistence::PersistenceAlg(_curve, beta, -1);
  _curve = p_result.pruned;
  for(int i = 0; i < iterations; ++i){
    // if(i != 0){
    //   recalc_angles_inplace(_curve);
    // }
    _curve = persistence::prune_curve_dist_to_segment(_curve, epsilon);
  }
  //print_curve(_curve);
  return _curve;
}


}//*** END OF NAMESPACE PERSISTENCE ***//



// [[Rcpp::export]]
NumericVector persistence_curve(NumericMatrix T) {
  std::vector<std::vector<double>> trajectory,query;

  for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
    trajectory.push_back({T(i,0),T(i,1)});

  auto curve = persistence::traj_to_curve(trajectory);
  auto res = std::vector<double>();
  for(auto &c: curve){
    res.push_back(c.val);
  }
  return wrap(res);
}

// [[Rcpp::export]]
NumericVector persistence_extrema(NumericMatrix T) {
  std::vector<std::vector<double>> trajectory,query;
  
  for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
    trajectory.push_back({T(i,0),T(i,1)});
  
  auto curve = persistence::traj_to_curve(trajectory);
  auto ext = persistence::extrema(curve);
  auto res = std::vector<double>();
  for(auto &c: ext.min){
    res.push_back(c);
  }
  for(auto &c: ext.max){
    res.push_back(c);
  }
  return wrap(res);
}

// [[Rcpp::export]]
NumericVector persistence_bars(NumericMatrix T, NumericVector it = -1) {
  persistence::Trajectory trajectory;

  for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
   trajectory.push_back({T(i,0),T(i,1)});

  auto curve = persistence::traj_to_curve(trajectory);
  auto p_result = persistence::PersistenceAlg(curve, 0, it(0));

  std::vector<int> bars;
  for(auto b : p_result.bars){
    bars.push_back( b.min);
    bars.push_back( b.max);
  }

  return wrap(bars);
}

// [[Rcpp::export]]
NumericMatrix persistence_pruned(NumericMatrix T, NumericVector Beta = 0,  NumericVector it = -1) {
  persistence::Trajectory trajectory;

  for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
    trajectory.push_back({T(i,0),T(i,1)});

  auto curve = persistence::traj_to_curve(trajectory);
  auto p_result = persistence::PersistenceAlg(curve, Beta(0), it(0));

  NumericMatrix resultMatrix = NumericMatrix(p_result.pruned.size(), 2) ;
  
  for(int i = 0; i < p_result.pruned.size(); i++){
    /*
     if(p_result[i].vertex.size() == 2){
     std::cout << "v["<<i<<"]: "<< p_result[i].vertex[0] <<", "<<p_result[i].vertex[1] << std::endl;
     }
     else{
     std::cout << "v["<<i<<"]: EMPTY!"<< p_result[i].vertex.size()<<" - " << p_result[i].index << std::endl;
     }
     */
    NumericVector temp  =  wrap(p_result.pruned[i].vertex);
    resultMatrix(i,_) = temp;
  }
  return resultMatrix;
}

// [[Rcpp::export]]
NumericVector persistence_comps(NumericVector T, NumericVector Beta = 0,  NumericVector it = -1) {
  persistence::Curve curve;

  for (size_t i=0; i < T.size(); i++)		//@todo: remove copy by an adapter class @Martin
    curve.push_back({T[i], 0, {0,0}});

  auto p_result = persistence::PersistenceAlg(curve, Beta(0), it(0));

  std::vector<int> comps;
  for(auto b : p_result.comps){
    comps.push_back( b.left);
    comps.push_back( b.min);
    comps.push_back( b.right);
  }

  return wrap(comps);
}

// [[Rcpp::export]]
NumericVector persistence_multires_index(NumericMatrix T, NumericVector Beta = 0, NumericVector Levels = 5) {
  persistence::Trajectory trajectory;

  for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
    trajectory.push_back({T(i,0),T(i,1)});

  auto curve = persistence::traj_to_curve(trajectory);
  auto p_result = persistence::persistenceMultiRes(curve, Beta(0), Levels(0));

  std::vector<int> res;
  for(auto r : p_result){
    res.push_back( r.index);
  }

  return wrap(res);
}

// [[Rcpp::export]]
NumericMatrix persistence_multires(NumericMatrix T, NumericVector Beta = 0, NumericVector Levels = 5) {
  persistence::Trajectory trajectory;
  try{
    for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
      trajectory.push_back({T(i,0),T(i,1)});
    
    auto curve = persistence::traj_to_curve(trajectory);
    
    auto p_result = persistence::persistenceMultiRes(curve, Beta(0), Levels(0));
    
    // std::vector<std::vector<double>> res;
    // for(auto r : p_result){
    //   res.push_back( r.vertex);
    // }
    // 
    // return wrap(res);
    // 
    
    NumericMatrix resultMatrix = NumericMatrix(p_result.size(), 2) ;
    for(int i = 0; i < p_result.size(); i++){
      /*
       if(p_result[i].vertex.size() == 2){
       std::cout << "v["<<i<<"]: "<< p_result[i].vertex[0] <<", "<<p_result[i].vertex[1] << std::endl;
       }
       else{
       std::cout << "v["<<i<<"]: EMPTY!"<< p_result[i].vertex.size()<<" - " << p_result[i].index << std::endl;
       }
       */
      NumericVector temp  =  wrap(p_result[i].vertex);
      resultMatrix(i,_) = temp;
    }
    return resultMatrix;
  }catch(int e){
    std::cout<<"ERROR: "<< e <<std::endl;
  }
  return NumericMatrix();
}
// [[Rcpp::export]]
NumericMatrix persistence_dist(NumericMatrix T, NumericVector Beta = 0, NumericVector Epsilon = 5, NumericVector Iterations = 5) {
  persistence::Trajectory trajectory;
  
    for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
      trajectory.push_back({T(i,0),T(i,1)});
    
    auto curve = persistence::traj_to_curve(trajectory);
    
    auto p_result = persistence::persistenceDist(curve, Beta(0), Epsilon(0), Iterations(0));
    
    NumericMatrix resultMatrix = NumericMatrix(p_result.size(), 2) ;
    
    for(int i = 0; i < p_result.size(); i++){
      /*
       if(p_result[i].vertex.size() == 2){
       std::cout << "v["<<i<<"]: "<< p_result[i].vertex[0] <<", "<<p_result[i].vertex[1] << std::endl;
       }
       else{
       std::cout << "v["<<i<<"]: EMPTY!"<< p_result[i].vertex.size()<<" - " << p_result[i].index << std::endl;
       }
       */
      NumericVector temp  =  wrap(p_result[i].vertex);
      resultMatrix(i,_) = temp;
    }
    return resultMatrix;
}

// [[Rcpp::export]]
NumericVector persistence_test_bars(NumericVector T, NumericVector Beta = 0, NumericVector it = -1) {
  persistence::Curve curve;
  
  for (size_t i=0; i < T.size(); i++)		//@todo: remove copy by an adapter class @Martin
    curve.push_back({T[i], 0, {0,0}});
  
  auto p_result = persistence::PersistenceAlg(curve, Beta(0), it(0));
  
  std::vector<int> bars;
  for(auto b : p_result.bars){
    bars.push_back( b.min);
    bars.push_back( b.max);
  }
  
  return wrap(bars);
}

// [[Rcpp::export]]
NumericVector persistence_test_comps(NumericVector T, NumericVector Beta = 0,  NumericVector it = -1) {
  persistence::Curve curve;
  
  for (size_t i=0; i < T.size(); i++)		//@todo: remove copy by an adapter class @Martin
    curve.push_back({T[i], 0, {0,0}});
  
  auto p_result = persistence::PersistenceAlg(curve, Beta(0), it(0));
  
  std::vector<int> comps;
  for(auto b : p_result.comps){
    comps.push_back( b.left);
    comps.push_back( b.min);
    comps.push_back( b.right);
  }
  
  return wrap(comps);
}